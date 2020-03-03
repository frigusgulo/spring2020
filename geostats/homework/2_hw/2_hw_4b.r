library(sp)
library(rgdal)
gravity = read.csv("~/spring2020/geostats/homework/2_hw/gravity255.txt")
xy = LongLatToUTM(gravity$long,gravity$lat,11)
x = xy$x
y = xy$y
out = point.crossval(x,y,gravity$gravity, 0.5,20000)