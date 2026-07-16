# ==================================================
# SWOR BOOSTED REGRESSION TREE (BRT) MODELLING PIPELINE
# ==================================================
#
# DESCRIPTION
# This script fits, evaluates, and validates Boosted Regression Tree (BRT)
# models for SWOR presence/absence data using environmental predictors
# derived from MOM6 ocean model outputs.
#
# The workflow includes:
#
# 1. Data preparation
#    - Load SWOR environmental dataset
#    - Format dates and response variable (PresAbs)
#
# 2. Model fitting (BRT using gbm.step)
#    - Single model fit
#    - 10-model ensemble to quantify uncertainty
#
# 3. Model interpretation
#    - Variable importance
#    - Partial dependence plots
#
# 4. Model evaluation
#    - Deviance explained
#    - AUC, TSS, Sensitivity, Specificity
#
# 5. Validation schemes
#    - 75/25 train-test split
#    - Leave-One-Year-Out (LOO) cross-validation
#    - Full-data evaluation (100%)
#
# 6. Output saving
#    - Model objects (.rds)
#    - Evaluation metrics (.csv)
#    - Figures (response curves, importance plots)
#
# NOTES
# - Designed for ecological niche modelling / fisheries applications
# - Uses gbm + dismo implementation of BRT (Elith et al.)
# - Response variable must be binary (PresAbs = 0/1)
#
# ==================================================

rm(list = ls())

# --------------------------------------------------
# LIBRARIES
# --------------------------------------------------

library(dplyr)
library(gbm)
library(dismo)
library(caret)
library(suncalc)
library(PresenceAbsence)

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

swor_monthly <- read.csv("FINAL/Gridded/1_Data/swor_env_MOM6_gridded_daily.csv")

data <- swor_monthly

data$date <- as.Date(data$date, format = "%Y-%m-%d")

# --------------------------------------------------
# OPTIONAL: FEATURE ENGINEERING (COMMENTED)
# --------------------------------------------------
# data$lunar <- getMoonIllumination(data$date)$fraction * 100

# --------------------------------------------------
# DEFINE PREDICTOR VARIABLES
# --------------------------------------------------

gbm.x <- c(
  "bbv", "ild", "sos", "ssh",
  "ssu_rotate", "ssv_rotate", "tos",
  "deptho", "tos_sd", "ssh_sd", "deptho_sd",
  "eke", "moon_phase"
)

# --------------------------------------------------
# TRAIN-TEST SPLIT (75/25)
# --------------------------------------------------

set.seed(123)

trainIndex <- createDataPartition(data$PresAbs, p = 0.75, list = FALSE)

train <- data[trainIndex, ]
test  <- data[-trainIndex, ]

# --------------------------------------------------
# BRT FITTING FUNCTIONS
# --------------------------------------------------

fit.brt <- function(data, gbm.x, gbm.y, lr) {
  
  gbm.step(
    data = data,
    gbm.x = gbm.x,
    gbm.y = gbm.y,
    family = "bernoulli",
    tree.complexity = 3,
    learning.rate = lr,
    bag.fraction = 0.60
  )
}

fit.brt.n10 <- function(data, gbm.x, gbm.y, lr, iterations = 10) {
  
  models <- vector("list", iterations)
  
  for (i in 1:iterations) {
    models[[i]] <- fit.brt(data, gbm.x, gbm.y, lr)
  }
  
  models
}

# --------------------------------------------------
# DEVIANCE EXPLAINED
# --------------------------------------------------

dev_eval <- function(model_object) {
  
  null <- model_object$self.statistics$mean.null
  res  <- model_object$self.statistics$mean.resid
  
  ((null - res) / null) * 100
}

# --------------------------------------------------
# SINGLE MODEL FIT
# --------------------------------------------------

set.seed(123)

model_single <- fit.brt(train, gbm.x, "PresAbs", lr = 0.05)

summary_res <- summary(model_single)

write.csv(summary_res,
          "summary_results_gridded_0.05.csv",
          row.names = FALSE)

dev_single <- dev_eval(model_single)

write.csv(dev_single,
          "summary_dev_single_gridded_0.05.csv",
          row.names = FALSE)

saveRDS(model_single,
        "SWOR_brt_single_gridded_0.05.rds")

# --------------------------------------------------
# 10-MODEL ENSEMBLE
# --------------------------------------------------

models_10 <- fit.brt.n10(train, gbm.x, "PresAbs", lr = 0.05)

saveRDS(models_10,
        "SWOR_brt_10models_gridded_0.05.rds")

devs <- sapply(models_10, dev_eval)

write.csv(
  data.frame(
    mean_dev = mean(devs),
    sd_dev   = sd(devs)
  ),
  "summary_dev_10models_gridded_0.05.csv",
  row.names = FALSE
)

# --------------------------------------------------
# MODEL VISUALISATION
# --------------------------------------------------

gbm.plot(model_single,
         smooth = TRUE,
         plot.layout = c(4, 4),
         write.title = TRUE)

dev.print(png,
          file = "partial_curves_gridded_0.05.png",
          width = 12, height = 8, units = "in", res = 300)

ggInfluence(model_single)

dev.print(png,
          file = "importance_noPh.png",
          width = 12, height = 8, units = "in", res = 300)

# --------------------------------------------------
# VALIDATION FUNCTIONS
# --------------------------------------------------

LOO_eval <- function(DataInput, gbm.x, gbm.y, lr = 0.01, tc = 3) {
  
  DataInput$Year <- as.numeric(format(as.Date(DataInput$date), "%Y"))
  
  years <- sort(unique(DataInput$Year))
  
  out <- data.frame(
    Year = years,
    Deviance = NA,
    AUC = NA,
    TSS = NA
  )
  
  for (i in seq_along(years)) {
    
    y <- years[i]
    
    train <- DataInput %>% filter(Year != y)
    test  <- DataInput %>% filter(Year == y)
    
    model <- gbm.step(
      data = train,
      gbm.x = gbm.x,
      gbm.y = gbm.y,
      family = "bernoulli",
      tree.complexity = tc,
      learning.rate = lr,
      bag.fraction = 0.6
    )
    
    preds <- predict.gbm(
      model,
      test,
      n.trees = model$gbm.call$best.trees,
      type = "response"
    )
    
    dev <- calc.deviance(test$PresAbs, preds, calc.mean = TRUE)
    
    e <- evaluate(p = preds[test$PresAbs == 1],
                  a = preds[test$PresAbs == 0])
    
    out$Deviance[i] <- dev
    out$AUC[i] <- e@auc
    out$TSS[i] <- max(e@TPR + e@TNR - 1)
  }
  
  out
}

eval_7525 <- function(DataInput, gbm.x, gbm.y, lr = 0.01, tc = 3) {
  
  idx <- sample(nrow(DataInput), 0.75 * nrow(DataInput))
  
  train <- DataInput[idx, ]
  test  <- DataInput[-idx, ]
  
  model <- gbm.step(train, gbm.x, gbm.y,
                    family = "bernoulli",
                    tree.complexity = tc,
                    learning.rate = lr,
                    bag.fraction = 0.6)
  
  preds <- predict.gbm(model, test,
                       n.trees = model$gbm.call$best.trees,
                       type = "response")
  
  e <- evaluate(p = preds[test$PresAbs == 1],
                a = preds[test$PresAbs == 0])
  
  data.frame(
    AUC = e@auc,
    TSS = max(e@TPR + e@TNR - 1)
  )
}

eval_100_percent <- function(DataInput, gbm.x, gbm.y, lr = 0.01, tc = 3) {
  
  model <- gbm.step(DataInput, gbm.x, gbm.y,
                    family = "bernoulli",
                    tree.complexity = tc,
                    learning.rate = lr,
                    bag.fraction = 0.6)
  
  preds <- predict.gbm(model, DataInput,
                       n.trees = model$gbm.call$best.trees,
                       type = "response")
  
  e <- evaluate(p = preds[DataInput$PresAbs == 1],
                a = preds[DataInput$PresAbs == 0])
  
  data.frame(
    AUC = e@auc,
    TSS = max(e@TPR + e@TNR - 1)
  )
}

# --------------------------------------------------
# RUN VALIDATION
# --------------------------------------------------

loo_eval  <- LOO_eval(data, gbm.x, "PresAbs", lr = 0.01, tc = 3)
eval_75   <- eval_7525(data, gbm.x, "PresAbs", lr = 0.05, tc = 3)
eval_100  <- eval_100_percent(data, gbm.x, "PresAbs", lr = 0.01, tc = 3)

# --------------------------------------------------
# SAVE VALIDATION RESULTS
# --------------------------------------------------

write.csv(loo_eval,  "SDR_LOO_eval_gridded.csv", row.names = FALSE)
write.csv(eval_75,   "SDR_7525_eval_gridded.csv", row.names = FALSE)
write.csv(eval_100,  "SDR_100_eval_gridded.csv", row.names = FALSE)

saveRDS(loo_eval,  "SDR_LOO_eval.rds")
saveRDS(eval_75,   "SDR_7525_eval.rds")
saveRDS(eval_100,  "SDR_100_eval.rds")