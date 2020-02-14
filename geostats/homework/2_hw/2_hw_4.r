library(rgdal)
library(sp)

gravity = read.table("~/spring2020/geostats/homework/2_hw/gravity255.txt",header=T,sep=",")
source("~/spring2020/geostats/scripts/point.crossval.r")
source("~/spring2020/geostats/scripts/LongLatToUTM.r")
source("~/spring2020/geostats/scripts/triangfuncs.r") # Loads triangulation functions for use with polydec.
source("~/spring2020/geostats/scripts/polydec.r")# Loads R-script with all of my functions
source("~/spring2020/geostats/scripts/celldec.r")


xy <- LongLatToUTM(gravity$long, gravity$lat,11)
p <- 2
r <- c(20000,40000,60000,80000,100000)

for (i in 1:5){
  point.crossval(xy$X,xy$Y,gravity$gravity,2,r[i])
}
