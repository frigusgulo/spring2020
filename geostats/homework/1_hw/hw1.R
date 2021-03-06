library(maps)
library(akima)
source("~/spring2020/geostats/scripts/image.legend.r")
source("~/spring2020/geostats/scripts/hscatter.r")
source("~/spring2020/geostats/scripts/corrplot.r")
source("~/spring2020/geostats/scripts/movewin.r")

#==========================================
# Swiss Contour
swiss <- read.table("~/spring2020/geostats/datasets/swiss.txt",header=T)
swiss.conc <- interp(swiss$long,swiss$lat,swiss$rainfall,nx=150,ny=150)

srfmean <- mean(swiss$rainfall)
srfvar <- var(swiss$rainfall)
srfrange <- range(swiss$rainfall)
#swissdf = data.frame(swiss)
image(swiss.conc,xlab="long",ylab="lat",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Greyscale Map of Rainfall",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
map("world","switzerland",lwd=4,add=T)
image.legend(6,47.8,zlim=range(swiss$rainfall),       # Puts a legend on the map with upper
             col=rev(heat.colors(24)))                     #   left corner at (3600,25).                  # Overlays concentration locations.

# Swiss histogram
#brk <- seq(0,max(swiss$rainfall),20)
hist(swiss$rainfall, xlab="Depth(cm)",ylab="Frequency",
     cex.lab=1.6,axes=F,       #   with no axes drawn (axes=F), axis
     main="Histogram of Rainfall Concentrations", #   labels, and a title using the
     cex.main=1.6,mgp=c(2.7,1,0),density=24,    #   break points in "breaks".
     angle=45)
axis(1,pos=0,cex.axis=1.5)                   # Appends an x-axis at position y=0
axis(2,pos=0,cex.axis=1.5)

par(mfrow=c(3,3),oma=c(0,0,0,0),            # Creates a 3x3 graphics window with no outer margin
    mar=c(0,0,3,0))                           #   (oma), and only a top margin (3) on each plot.
thresh <- round(quantile(swiss$rainfall,   # Computes indicator threshholds at each decile of
                         probs=seq(.1,.9,.1)),2)                   #   of the coal ash percentages (10%, 20%, ..., 90%)
for (i in 1:9){
  swissdf <- data.frame(swiss$long,swiss$lat,
                        swiss$rainfall>thresh[i])
  thresh.int <- interp(swiss$long,swiss$lat,swissdf[,3],nx=100,ny=100)
  image(thresh.int,xlim=c(6.0,10.5),ylim=c(45.8,47.8),col=c(0,1),
        axes=F, main=paste("Rainfall % >",thresh[i], sep=""))
}


#========================================================================================================
# Wolf Camp Data
library(geoR)
data(wolfcamp)
wcmean <- mean(wolfcamp$data)
wcvar <- var(wolfcamp$data)
wcrange <- range(wolfcamp$data)
coords = wolfcamp$coords
wolfcampwin <- movewin(coords[,1],coords[,2],wolfcamp$data,70,70,10,10)

wolfcamptrans <- log(wolfcamp$data)
wolfcamp.conc <- interp(coords[,1], coords[,2],wolfcamptrans,nx=150,ny=150)
par(mfrow=c(1,1),mar=c(5,4,4,2))
image(wolfcamp.conc,xlab="Easting",ylab="Northing",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Piezometric-Head Heights (meters above sea level)",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
image.legend(-220,130,zlim=range(wolfcamptrans),
             col=rev(heat.colors(24)))
hist(wolfcamptrans, xlab="Meters Above Sea Level",ylab="Frequency",
     cex.lab=1.6,axes=F,       #   with no axes drawn (axes=F), axis
     main="Histogram of Piezometric Head Heights", #   labels, and a title using the
     cex.main=1.6,mgp=c(2.7,1,0),density=24,    #   break points in "breaks".
     angle=45)
axis(1,pos=0,cex.axis=1.5)                   # Appends an x-axis at position y=0
axis(2,pos=0,cex.axis=1.5)
#corrplot <- function(x,y,u,v,h,numlags,dtol=NA,atol=NA){
corrplot(x=coords[,1],y=coords[,2],wolfcamptrans,wolfcamptrans,
         h=c(20,0),5,dtol=10,atol=15)
corrplot(x=coords[,1],y=coords[,2],wolfcamptrans,wolfcamptrans,
         h=c(0,20),5,dtol=10,atol=15)

#========================================================================================================
# Phytophthora Data
phyto <- read.table("~/spring2020/geostats/datasets/phytoph.txt",header=T)
phy_perc_diz <-100*(sum(phyto[,3]==1)/sum(phyto[,3]==0))
phy_moise_mean <- mean(phyto[,4])
phy_moise_var <- var(phyto[,4])
phy_moist_range <- range(phyto[,4])
phyto.conc <- interp(phyto[,1],phyto[,2],phyto[,4],nx=150,ny=150)
image(phyto.conc,xlab="X",ylab="Y",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Soil Moisture Percentage",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(phyto[,1],phyto[,2],cex=1.5)
points(phyto[,1][phyto[,3]==1],phyto[,2][phyto[,3]==1],pch=16,cex=1.5)
corrplot(phyto[,1],phyto[,2],phyto[,4],phyto[,3],h=c(1,0),5,dtol=1,atol=15)
corrplot(phyto[,1],phyto[,2],phyto[,4],phyto[,3],h=c(0,1),5,dtol=1,atol=15)
#corrplot <- function(x,y,u,v,h,numlags,dtol=NA,atol=NA){
#================================
# Chorley Cancer
library(spatstat)
data("chorley")
incin <- chorley.extra$incin
larynx_perc <- 100*(sum(chorley$marks == "larynx")/(sum(chorley$marks=="larynx") + sum(chorley$marks=="lung")))
lung_perc <- 100 - larynx_perc
quadrat.test(chorley[chorley$marks == "larynx"], nx=5, ny=5,
             alternative=c( "clustered"),
             method=c("MonteCarlo"),
             conditional=TRUE, CR=1,
             nsim=1999)
quadrat.test(chorley[chorley$marks == "lung"], nx=5, ny=5,
             alternative=c( "clustered"),
             method=c("MonteCarlo"),
             conditional=TRUE, CR=1,
             nsim=1999)
mean_diff <- mean(sqrt((chorley$x - incin$x)^2 + (chorley$y - incin$y)^2)) # calculate average distance to incinerator
window <- as.owin(c(346.6, 364.1, 412.6, 430.3))
samples <- 100
csr <- rpoispp(lambda=1,win=window,nsim=samples)
means <- integer(samples)
for (i in 1:samples){
  ind = paste("Simulation",i,sep=" ")
  val <- mean(sqrt((csr[[ind]][["x"]] - incin$x)^2 + (csr[[ind]][["y"]] - incin$y)^2))
  means[i] <- val
}
probs_close = pnorm(mean_diff, mean=mean(means),sd=sd(means),lower.tail=TRUE)

chorley.extra$plotit()

