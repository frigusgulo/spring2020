library(tripack)                # Needed for "tri.mesh" triangulation
library(sp)                     # Needed for "point.in.polygon" function
library(splancs)                # Needed for "areapl" function
library(spatstat)               # Needed for "ppp" and "delaunay" functions
library(maps)
library(akima)
library(pracma)
library(e1071)   
source("~/spring2020/geostats/scripts/triangfuncs.r") # Loads triangulation functions for use with polydec.
source("~/spring2020/geostats/scripts/polydec.r")# Loads R-script with all of my functions
source("~/spring2020/geostats/scripts/celldec.r")
source("~/spring2020/geostats/scripts/image.legend.r")
source("~/spring2020/geostats/scripts/allfunctions.r")


snow <- read.table("~/spring2020/geostats/homework/2_hw/snowdep.txt",header=T,sep = ",")

# 3 a - box plot, histogram, heatmap, correlogram (both axes)
snow.conc <- interp(snow$x,snow$y,snow$depth,nx=1500,ny=1500)

image(snow.conc,xlab="long",ylab="lat",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Greyscale Map of Snowfall",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
image.legend(6,47.8,zlim=range(snow$depth),       # Puts a legend on the map with upper
             col=rev(heat.colors(24)))

boxplot(snow$depth)


hist(snow$depth, xlab="Depth(cm)",ylab="Frequency",
     cex.lab=1.6,axes=F,       #   with no axes drawn (axes=F), axis
     main="Histogram of snowfall measurements", #   labels, and a title using the
     cex.main=1.6,mgp=c(2.7,1,0),density=24,    #   break points in "breaks".
     angle=45)
axis(1,pos=0,cex.axis=1.5)                   # Appends an x-axis at position y=0
axis(2,pos=0,cex.axis=1.5)


x=snow$x
y=snow$y
v = snow$depth
u = snow$depth
dtol = 0.25
atol = 25
h <- c(0,1,2,3,4,5,6)                              # Sets a vector of h-values.
out <- matrix(nrow=length(h),ncol=3)               # Defines an hx3 blank matrix.
for (i in 0:6) out[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,u,v,c(0,i),dtol=dtol,atol=atol))                        #   calculates cross-functions for each.
cch <- out[,1]
cph <- out[,2 ]
cgh <- out[,3]     # Defines the vectors of cross-covariances,
#   cross-corr's, & cross-semivariograms.
par(mfrow=c(2,2))                                  # Sets up a 2x2 graphics window.
plot(h,cch,type="n",axes=F,xlab="",ylab="")        # Creates a completely blank plot.
plot(h,cch,type="l",xlab="|h|",ylab="Cross - C(h)",# Plots the cross-correlations vs. h
     main="Cross-Covariogram for U,V",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h,cph,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Cross-Correlogram for U,V",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h[-1],cgh[-1],type="l",xlab="|h|",cex.lab=1.5,# Plots the cross-semivariogram vs. h
     ylab="Cross - Gamma(h)",cex.axis=1.3,            #   with sizes controlled by "cex".
     main="Cross-Semivariogram for U,V",cex.main=1.5)

out_y <- matrix(nrow=length(h),ncol=3)               # Defines an hx3 blank matrix.
for (i in 0:6) out_y[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(i,0),dtol=dtol,atol=atol))                        #   calculates cross-functions for each.
cch_y <- out_y[,1]
cph_y <- out_y[,2]
cgh_y <- out_y[,3]     # Defines the vectors of cross-covariances,
#   cross-corr's, & cross-semivariograms.
par(mfrow=c(2,2))                                  # Sets up a 2x2 graphics window.
plot(h,cch_y,type="n",axes=F,xlab="",ylab="")        # Creates a completely blank plot.
plot(h,cch_y,type="l",xlab="|h|",ylab="Cross - C(h)",# Plots the cross-correlations vs. h
     main="Cross-Covariogram for U,V",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h,cph_y,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Cross-Correlogram for U,V",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h[-1],cgh_y[-1],type="l",xlab="|h|",cex.lab=1.5,# Plots the cross-semivariogram vs. h
     ylab="Cross - Gamma(h)",cex.axis=1.3,            #   with sizes controlled by "cex".
     main="Cross-Semivariogram for U,V",cex.main=1.5)



# 3 b - sample mean, discuss biased (it is, samples were clustered)

snowmean = mean(snow$depth)
snowskew = skewness(snow$depth)

# 3 c - polygonal declustering

pd = polydec(snow$x,snow$y,peel=7)
pd = pd / sum(pd)
pd_mean = sum(snow$depth*pd)


par(mfrow=c(1,1))
# 3 d - cell declustering
cellMeans <- matrix(nrow=30,ncol=30)  
for (i in 1:30){
  for (j in 1:30){
  cellMeans[i,j] = celldec(snow$x,snow$y,snow$depth,i,j)

  }
}
rotate_one_eighty <- function(x) { rot90(x, 2) }
cellMeans = rotate_one_eighty(cellMeans)
contour(x=seq(1,30),y=seq(1,30),z=cellMeans)
