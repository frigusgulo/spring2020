library(rgdal)
library(sp)
library(splancs)
gravity = read.table("~/spring2020/geostats/homework/2_hw/gravity255.txt",header=T,sep=",")
source("~/spring2020/geostats/scripts/point.crossval.r")
source("~/spring2020/geostats/scripts/LongLatToUTM.r")
source("~/spring2020/geostats/scripts/triangfuncs.r") # Loads triangulation functions for use with polydec.
source("~/spring2020/geostats/scripts/polydec.r")# Loads R-script with all of my functions
source("~/spring2020/geostats/scripts/celldec.r")


xy <- LongLatToUTM(gravity$long, gravity$lat,11)
p <- 2
r <- c(20000,40000,60000,80000,100000)

x = xy$X
y = xy$Y
z = gravity$gravity
for (i in 1:5){
  point.crossval(xy$X,xy$Y,gravity$gravity,2,r[i])
}

p_ = c(0.2,0.5,1.0,2.0,5.0,10.0)
for (i in 1:length(p_)) {
  point.crossval(x,y,z,p_[i],20000)
}


#QQ plots
out <- point.crossval(xy$X,xy$Y,gravity$gravity,1.0,40000)
qqplot( gravity$gravity,out$pred.pd)
qqplot( gravity$gravity,out$pred.tr)
qqplot( gravity$gravity,out$pred.ls)
qqplot( gravity$gravity,out$pred.iw)

#scatters
plot(gravity$gravity,out$pred.pd,main="Prediction Polygonal Declustering")
lines(lowess( gravity$gravity,out$pred.pd))

plot( gravity$gravity,out$pred.tr,main="Prediction Triangulation")
lines(lowess( gravity$gravity,out$pred.tr))

plot( gravity$gravity,out$pred.ls,main="Prediction Local Sample Mean")
lines(lowess(gravity$gravity,out$pred.ls))

plot(gravity$gravity,out$pred.iw,main="Prediction Inverse Distance")
lines(lowess(gravity$gravity,out$pred.iw))



