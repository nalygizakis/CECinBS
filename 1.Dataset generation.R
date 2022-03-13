#Clean workspace
rm(list=ls(all=TRUE))

#Load dependencies
library(gstat)
library(mgcv)
library(maps)
library(mapdata)
library(RColorBrewer)
library(akima)
library(maptools)
library(classInt)
library(scales)
library(plot3D)
library(magick)

#Load the dataset
load("1.Dataset generation_workspace.RData")

#Selection of sampling year
areas <- as.data.frame(seawater_survey2016)
selected_year <- "2016"

#Define functions
mefilled.contour <- function(x = seq(0, 1, length.out = nrow(z)), 
                             y = seq(0, 1, length.out = ncol(z)), 
                             z, 
                             xlim = range(x, finite = TRUE), 
                             ylim = range(y, finite = TRUE), 
                             zlim = range(z, finite = TRUE), 
                             levels = pretty(zlim, nlevels),
                             nlevels = 20, 
                             color.palette = cm.colors, 
                             col = color.palette(length(levels) - 1),
                             plot.title, 
                             plot.axes, 
                             key.title, 
                             key.axes, 
                             asp = NA, 
                             xaxs = "i", 
                             yaxs = "i",
                             las = 1, 
                             axes = TRUE, 
                             frame.plot = axes, ...){
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54

  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1

  plot.new()

  if (missing(key.axes)) {
    if (axes) 
      print(1)

  }
  else key.axes

  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)

  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

plotMapPoints<-function(area=area, 
                        xlim=c(19.5,29.5), 
                        ylim=c(34.8,41.8), 
                        col="grey60", 
                        fill=T,
                        type ="p", 
                        col_points="red",cex=0.5,
                        plotclr=c("pink1","red", "red4","black"), 
                        symbol.size=1.1){
  library(gstat)
  library(mgcv)
  library(maps)
  library(mapdata)
  library(RColorBrewer)
  library(akima)
  library(maptools)
  library(classInt)
  library(scales)

  min(area$Abundance)
  max(area$Abundance)
  min(log(area$Abundance))
  max(log(area$Abundance))
  

  plotvar <- area$Abundance

  plotclr<-plotclr
  nclr <- length(plotclr)
  
  class <- classIntervals(plotvar, n=nclr, style="quantile")

  class

  colcode <- findColours(class, plotclr)
  
  map(database = "world", xlim=xlim, ylim = ylim, resolution = 0, col=col, fill=fill,bg="white",
      xlab="Longitude",ylab="Latitude")
  map.axes()
  points(area$LON, area$LAT, pch=16, type="p",col=colcode,cex=symbol.size,xlab="Longitude",ylab="Latitude")

  legend<-names(attr(colcode, "table"))

  
  text(34.5,41,"Turkey", cex = 1)
  text(29.7,46.2,"Ukraine", cex = 1)
  text(40.5,45.0,"Russia", cex = 1)
  text(42.4,42.4,"Georgia", cex = 1)
  text(27.7,43.6,"Bulgaria", cex = 1)
  text(27.8,44.6,"Romania", cex = 1)  
  text(29,45.18,"Danube", cex = 2.0,font=2, col="#ADD8E6")  
}

ylim=c(44,47)
xlim=c(27,34)

pred.points<-rbind(
  cbind(x=seq(from=29,to=30,length.out=length(seq(from=44.7, to=45.5,by=0.01))),seq(from=44.7, to=45.5,by=0.01)),
  cbind(x=seq(from=29,to=30,length.out=length(rev(seq(from=44.7, to=45.5,by=0.01)))),rev(seq(from=44.7, to=45.5,by=0.01))),
  
  cbind(x=seq(from=30,to=31,length.out=length(seq(from=44, to=44.5,by=0.01))),seq(from=44, to=44.5,by=0.01)),
  cbind(x=seq(from=30,to=31,length.out=length(rev(seq(from=44, to=44.5,by=0.01)))),rev(seq(from=44, to=44.5,by=0.01))),
  
  cbind(x=seq(from=31,to=32,length.out=length(seq(from=44, to=47,by=0.01))),seq(from=44, to=47,by=0.01)),
  cbind(x=seq(from=31,to=32,length.out=length(rev(seq(from=44, to=47,by=0.01)))),rev(seq(from=44, to=47,by=0.01))),
  
  cbind(x=seq(from=32,to=33,length.out=length(seq(from=44, to=45,by=0.01))),seq(from=44, to=45,by=0.01)),
  cbind(x=seq(from=32,to=33,length.out=length(rev(seq(from=44, to=45,by=0.01)))),rev(seq(from=44, to=45,by=0.01))),
  
  cbind(x=seq(from=33,to=34,length.out=length(seq(from=42.5, to=43,by=0.01))),seq(from=42.5, to=43,by=0.01)),
  cbind(x=seq(from=33,to=34,length.out=length(rev(seq(from=42.5, to=43,by=0.01)))),rev(seq(from=42.5, to=43,by=0.01))),
  
  cbind(x=seq(from=34,to=35,length.out=length(seq(from=42.5, to=43,by=0.01))),seq(from=42.5, to=43,by=0.01)),
  cbind(x=seq(from=34,to=35,length.out=length(rev(seq(from=42.5, to=43,by=0.01)))),rev(seq(from=42.5, to=43,by=0.01))),
  
  cbind(x=seq(from=35,to=36,length.out=length(seq(from=42.5, to=43,by=0.01))),seq(from=42.5, to=43,by=0.01)),
  cbind(x=seq(from=35,to=36,length.out=length(rev(seq(from=42.5, to=43,by=0.01)))),rev(seq(from=42.5, to=43,by=0.01)))
)

pred.points<-as.data.frame(pred.points)
names(pred.points)<-c("LON","LAT")

#####1st figure (Point plot)-----
i<-6
dir.create("point_plots")
pdf("point_plots/pointplot.pdf", width = 7, height = 5, paper="special", fillOddEven=T)
input<-t(areas)
for(i in 6:ncol(areas)){
  
  area<-data.frame(LAT=as.numeric(areas[,2]),
                   LON=as.numeric(areas[,3]),
                   Abundance=as.numeric(areas[,c(i)]),
                   mzrt=names(areas)[c(i)])
  
  plotMapPoints(area,ylim=c(44,47), xlim=c(27,34), col="white", fill=TRUE,
                type ="p", col_points="red", cex=1.2, 
                plotclr=c("#f2f0f7","#cbc9e2", "#9e9ac8","#756bb1","#54278f"),
                symbol.size=1.8)
  
  print(i)
}
dev.off()


#####2nd plot (spatial distribution)-----
dir.create("spatial_distribution")
pdf("spatial_distribution/spatial_distribution.pdf", width = 7, height = 5,
    paper="special", fillOddEven=T)
i<-6
for(i in 6:ncol(areas)){
  
  input<-data.frame(year=rep(selected_year,length(areas[,1])),
                    LAT=as.numeric(areas[,2]),
                    LON=as.numeric(areas[,3]),
                    Abundance=areas[,c(i)],
                    mzrt=names(areas)[c(i)],
                    names=areas$Sample)
  input$year<-as.factor(input$year)
  
  area<-input
  area$Abundance[area$Abundance==0]<-1
  m1 <- gam(Abundance~LON+LAT,family=Gamma(link=log),data=area)

  model.res <- resid(m1, type="deviance")
  res.grid <- interp(area$LON, area$LAT, model.res, duplicate="median")
  mypalette <- palette(gray(seq(.9,.4,len=6)))
  
  mydata2 <- data.frame(model.res, area$LON, area$LAT)
  coordinates(mydata2) <- c("area.LON", "area.LAT")

  pred <- predict(m1, pred.points, type="link")

  min(pred);max(pred)
  pred <- rescale(pred, c(1,100))
  mygrid  <- interp(x=pred.points$LON, y=pred.points$LAT, z=pred, duplicate="mean",linear=FALSE)
  mydataare <- cbind(pred.points$LON,pred.points$LAT,pred)

  mefilled.contour(x=mygrid, 
                   xlim=xlim, 
                   ylim=ylim, 
                   zlim=c(0,100),
                   col=rev(heat.colors(10)), nlev=10,
                   xlab="", ylab="",
                   axes=FALSE,
                   frame.plot=FALSE,
                   main="",
                   plot.axes=(map(database="world2", 
                                  add=T, 
                                  resolution=0, 
                                  fill=TRUE,
                                  border=NA,
                                  col="white",bg="light blue")))
  print(i)
 }
dev.off()


#Open the pdf with Adobe reader and export as png images----
files_for_processing<-list.files("./spatial_distribution", pattern=".png", full.names = TRUE)
files_for_processing_short<-list.files("./spatial_distribution", pattern=".png", full.names = FALSE)
dir.create("./spatial_distribution/png_cropped")
new_path<-'./spatial_distribution/png_cropped'
i<-1
for(i in 1:length(files_for_processing)){
  temp_image <- image_read(files_for_processing[i])
  temp_image <-image_crop(temp_image, "870x810+500+150") 
  temp_image <-image_trim(temp_image)
  image_write(temp_image, 
              path = paste0(new_path,"/",files_for_processing_short[i]), 
              format = "png")
  print(i)
}