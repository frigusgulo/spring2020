library(gstat)
library(lattice) 
library(sp)
library(geoR)
library(map)
library(usmap)
library(TSdist)
source("~/spring2020/geostats/scripts/allfunctions.r")
source("~/spring2020/geostats/scripts/getCD.r")

x = c(130,150,170,156,151,161,137,160,139,145)
y = c(231,228,234,228,218,241,240,218,220,234)
v = c(530.3,402.3,173.1,458.4,594.3,141.9,112.5,580.4,535.9,333.2)

par(mfrow=c(1, 1))
plot(x, y, type="n",
     xlab='x', ylab='y', main='Walker Lake V')
text(x, y, format(round(v, 1)),
     cex=0.7)
########################

#PART B
# determine azi of maximum continuity
angles = seq(0,165,15)
df = data.frame(x=x,y=y,data=v)
coordinates(df) = ~x+y
walk.vario = variogram(data~1,data=df,width=1,cutoff=20,alpha=angles,tol.hor = 35 )
plot(walk.vario,main='Directional Variograms For Walker Lake Data' )

azi = 90
nug = 0
sill = 2000
range = 30
CD = getCD(x,y,155,235,nug,sill,range,azi,0,"spher")
C = CD$C
D = CD$D
ones = rep(1,length(x)+1)
Cinv = solve(C)
A = (Cinv%*%ones%*%t(ones)%*%Cinv%*%D)*(1/(t(ones)%*%Cinv%*%ones))[1]
B = Cinv%*%ones*( (1/t(ones)%*%Cinv%*%ones)[1])
Z = Cinv%*%D
w1 = Z - A + B

ones = rep(1,length(x))
Cinv = solve(CD$C[1:10,1:10])
MGLS = t(ones)%*%Cinv%*%v / t(ones)%*%Cinv%*%ones
Var = sum(v - MGLS)
STD = sqrt(Var)

#########################################
# Part C

euclideanDistance = function(vec1,vec2){
  dist = rep(0,NROW(vec1))
  z = NROW(vec1)
  for (i in 1:z){
    dist[i] = sqrt( (vec1[i,1] -vec2[1,1])^2 + (vec1[i,2] -vec2[1,2])^2 )
  }
  return(dist)
}

dist = euclideanDistance(cbind(x,y),cbind(155,235))
IDW = 1 / dist^2
pred = IDW%*%v

################################
# Part C
locx = 155
locy = 235
range = c(20,40,50)