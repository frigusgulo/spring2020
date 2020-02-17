library(rgdal)
library(sp)

gravity = read.table("~/spring2020/geostats/homework/2_hw/gravity255.txt",header=T,sep=",")
source("~/spring2020/geostats/scripts/point.crossval.r")
source("~/spring2020/geostats/scripts/LongLatToUTM.r")
source("~/spring2020/geostats/scripts/triangfuncs.r") # Loads triangulation functions for use with polydec.
source("~/spring2020/geostats/scripts/polydec.r")# Loads R-script with all of my functions
source("~/spring2020/geostats/scripts/celldec.r")


xy = LongLatToUTM(gravity$long, gravity$lat,11)
x = xy$X; y = xy$Y
z = gravity$gravity
p = 2
r = c(20000,40000,60000,80000,100000)

for (i in 1:5){
  point.crossval(x,y,z,p,r[i])
}

p_ = c(0.2,0.5,1.0,2.0,5.0,10.0)
for (i in 1:length(p_)) {
  point.crossval(x,y,z,p_[i],20000)
}


#QQ plots
out <- point.crossval(x,y,z,0.5, 30000)

plot(sort(out$pred.pd), sort(gravity$gravity))
plot(sort(out$pred.tr), sort(gravity$gravity))
plot(sort(out$pred.ls), sort(gravity$gravity))
plot(sort(out$pred.is), sort(gravity$gravity))



# f
# polydec
pd = polydec(snow$x,snow$y,peel=7)
pd = pd / sum(pd)
pd_mean = sum(snow$depth*pd)

#celldec
celldec(snow$x,snow$y,snow$depth,i,j)

#idw

#triangulation
