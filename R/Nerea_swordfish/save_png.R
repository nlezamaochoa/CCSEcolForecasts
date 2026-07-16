library(ncdf4)
# library(maps)
# library(mapdata)



### generalized function for operationalization
make_png_operationalization=function(r,spname, get_date,outDir,type){ ### does what it says
  
  png(paste0(outDir,spname,"/",spname,"_", get_date,"_",type,".png"), width=5, height=7, units="in", res=400)
  par(ps=10) #settings before layout
  layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
  #layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
  par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
  
  par(mar=c(4,4,1,1)) # I usually set my margins before each plot
  #pal <- colorRampPalette(c("blue", "grey", "red"))
  #pal <- colorRampPalette(c("white","blue", "cyan", "yellow", "red"))
   pal <- colorRampPalette(c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c")) #Revised palette
  # ncolors <- 100
  #pal <- colorRampPalette(c("purple4", "white", "blue"))
  ncolors <- 100
  #min<- minValue(meanPredR)
  #max <- maxValue(meanPredR)
 # breaks <- seq(0.02,0.3,,ncolors+1)
  #breaks <- round(breaks, digits = 6)
  breaks2 <- seq(0,1,,ncolors+1)
 # breaks <- c(-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.45,-0.4,-0.36,seq(-0.08573341,0.2963,,83),0.36,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1)
  image(r,col=pal(ncolors), breaks2=breaks2,ylab="", xlab="")
  maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
 # contour(r, add=TRUE, col="black",levels=c(.5,.75))
  box()
  
  par(mar=c(4,4,0,1)) # I usually set my margins before each plot
  levs <- breaks2[-1] - diff(breaks2)/2
  image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks2, ylab="", xlab="", yaxt="n")
  mtext(paste0("Probability of ", spname," presence ",get_date," (",type,")",sep=" "), side=1, line=2.5)
  
  box()
  
  dev.off() # closes device
}


make_png_operationalization_SE=function(r,spname, get_date,outDir,type){ ### does what it says
  
  png(paste0(outDir,spname,"/",spname,"_", get_date,"_",type,".png"), width=5, height=7, units="in", res=400)
  par(ps=10) #settings before layout
  layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), heights=c(4,1), widths=7)
  #layout.show(2) # run to see layout; comment out to prevent plotting during .pdf
  par(cex=1) # layout has the tendency change par()$cex, so this step is important for control
  
  par(mar=c(4,4,1,1)) # I usually set my margins before each plot
  pal <- colorRampPalette(c("coral3","coral","cadetblue3","cadetblue1","antiquewhite","cadetblue1","cadetblue3","coral","coral3"))
  ncolors <- 100
  breaks <- seq(-0.2,0.2,,ncolors+1)
  image(r, col=pal(ncolors), breaks=breaks, ylab="", xlab="", xlim=c(-130,-115.5),ylim=c(30,47))
  maps::map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
  contour(r, add=TRUE, col="black",levels=c(.5,.75))
  box()
  
  par(mar=c(4,4,0,1)) # I usually set my margins before each plot
  levs <- breaks[-1] - diff(breaks)/2
  image(x=levs, y=1, z=as.matrix(levs), col=pal(ncolors), breaks=breaks, ylab="", xlab="", yaxt="n")
  mtext(paste0("Probability of ", spname," presence ",get_date," (",type,")",sep=" "), side=1, line=2.5)
  
  box()
  
  dev.off() # closes device
}

