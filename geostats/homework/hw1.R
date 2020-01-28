library(maps)
library(akima)
source("~/spring2020/geostats/image.legend.r")

# Swiss Contour
swiss <- read.table("~/spring2020/geostats/datasets/swiss.txt",header=T)
swiss.conc <- interp(swiss$long,swiss$lat,swiss$rainfall,nx=150,ny=150)
#swissdf = data.frame(swiss)
image(swiss.conc,xlab="long",ylab="lat",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Greyscale Map of Rainfall",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
map("world","switzerland")
image.legend(0,0,zlim=range(swiss$rainfall),       # Puts a legend on the map with upper
             col=rev(heat.colors(24)))                     #   left corner at (3600,25).
points(swiss$long,swiss$lat,pch=16,cex=1.0)                      # Overlays concentration locations.
contour(swiss$rainfall, xlab="x",labcex=1.5,ylab="y", # Creates a contour plot with the
        cex.axis=1.5,cex.lab=1.6,cex.main=1.2,        #   interpolated logconc values.
        main="Rainfall")

# Swiss histogram
#brk <- seq(0,max(swiss$rainfall),20)
hist(swiss$rainfall, xlab="Depth(cm)",ylab="Frequency",
     cex.lab=1.6,axes=F,       #   with no axes drawn (axes=F), axis
     main="Histogram of Rainfall Concentrations", #   labels, and a title using the
     cex.main=1.6,mgp=c(2.7,1,0),density=24,    #   break points in "breaks".
     angle=45)
axis(1,pos=0,cex.axis=1.5)                   # Appends an x-axis at position y=0
axis(2,pos=0,cex.axis=1.5)



