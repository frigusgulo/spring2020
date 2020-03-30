library(gstat)
library(lattice) 
library(sp)
library(geoR)
library(maps)
library(usmap)
#library(spatial)
library(fields)
library(akima)
source("~/spring2020/geostats/scripts/allfunctions.r")
source("~/spring2020/geostats/scripts/rotateaxis.r")
source("~/spring2020/geostats/scripts/panel.gamma0.r")



scallops = read.table("~/spring2020/geostats/datasets/scallops.txt",header=T)
X = scallops$long
Y = scallops$lat
catch_t = scallops$tcatch
catch_t_log = log(catch_t+1)
summary(catch_t)
catch_t.conc = interp(X,Y,catch_t_log,nx=150,ny=150)
par(mfrow=c(1,3))
image(catch_t.conc,xlab="Longitude",ylab="Latitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(X,Y)
hist(catch_t,main="Total Catch")
hist(catch_t_log,main="Log Transformed Total Catch")

##########################

#  PART B

catch_t.conc = interp(X,Y,catch_t_log,nx=75,ny=75)
par(mfrow=c(1,1))
map("usa",xlim=c(-75,-71),ylim=c(38.2,41.5),fill=T,col=2)
points(scallops$long,scallops$lat,pch=16,cex=.75)
contour(catch_t.conc,add=T)
############################

# PART C

xy = cbind(X,Y)
coords_trans = rotateaxis(xy,52)
catch_t_loggamma = log(catch_t + 1)
x = coords_trans$newx
y = coords_trans$newy
w = y-min(y)
nlfit = nls(catch_t_loggamma ~ (delta*x + w^(alpha-1)*
  exp(-w/beta + gamma)) , start=list(delta=0,alpha=5.96,
  beta=0.07766,gamma=12.8755))
'''
solved for parameters with python:

import numpy as np                                                                
from numpy.linalg import solve                                                    
A = np.array([[np.log(0.25),-0.25,1],[np.log(0.5),-0.5,1],[np.log(1),-1,1]])      
b = np.array([np.log(1),np.log(5),np.log(1)])                                     
b = np.Transpose(b) 
x = solve(A,b)
x[0] -= 1
x[1] = 1/x[1] 

In [11]: x                                                                                
Out[11]: array([ 5.96578428,  0.07766687, 12.8755033 ])


'''
summary(nlfit)

# Could refer to parameters as:
# =============================
delta = summary(nlfit)$par[1]
alpha = summary(nlfit)$par[2]
beta = summary(nlfit)$par[3]
gamma = summary(nlfit)$par[4]
# =============================

zfunc = function(x,y){
delta = -0.08790   
alpha = 9.40226 
beta  = 0.05821 
gamma = 15.80569 
w = y-min(y)
z = (delta*x + w^(alpha-1)*exp(-w/beta + gamma))
return (z)
}
  
catch_t_log_trend = zfunc(x,y)
catch_t_log_resids = catch_t_log - catch_t_log_trend


catch_t_log.conc = interp(X,Y,catch_t_log,nx=150,ny=150)
catch_t_log_trend.conc = interp(x,y,catch_t_log_trend,nx=150,ny=150)
catch_t_log_resids.conc = interp(x[catch_t_log_resids!=-Inf],y[catch_t_log_resids!=-Inf],catch_t_log_resids[catch_t_log_resids!=-Inf],nx=150,ny=150)

par(mfrow=c(1,3))
image(catch_t_log.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(x,y)

image(catch_t_log_trend.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Inferred Trend of Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(x,y)

image(catch_t_log_resids.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Inferred Trend Residuals of Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(x,y)

#########################################

# PART D

par(mfrow=c(1,2))
plot(catch_t_log,catch_t_log_resids,ylab='Residuals',xlab='Trend',main='Log Transformed Total Catch Fitted Trend vs Residuals')
abline(0,0)
qqnorm(catch_t_log_resids[catch_t_log_resids!=-Inf],main="Normal Q-Q Plot",xlab= "Theoretical Quantiles")
qqline(catch_t_log_resids,distribution = qnorm, probs = c(0.25, 0.75) )

##########################################

# PART E
par(mfrow=c(1,2))
angles = seq(0,165,15)
df <- data.frame(x,y,data=catch_t_log_resids)
coordinates(df) = ~x+y
catch_t_log_resids.vario_dir = variogram(data~1,data=df,width=0.1,cutoff=1,alpha=angles,tol.hor = 30 )
plot(catch_t_log_resids.vario_dir,main='Directional Variograms For Log Transformed Catch Residuals' )
rose(catch_t_log_resids.vario_dir,2.3)


#############################################

# PART F

xyplot(gamma ~ dist|as.factor(dir.hor),
  data=catch_t_log_resids.vario_dir,layout=c(4,3), panel =
  panel.gamma0, gamma0=2.3, main=
  "Interpolated Directional Ranges")

coords = cbind(x,y)
anisopars = cbind(90,2.55)
coords_corrected = coords.aniso(coords,anisopars)

newx = coords_corrected[,1]
newy = coords_corrected[,2]
df_c = data.frame(x=newx,y=newy,data=catch_t_log_resids)
coordinates(df_c) = ~x+y
iso_var = variogram(data~1,data=df_c,width=0.1,cutoff=1)
par(mfrow=c(1,2))
plot(iso_var,main="Isotropic Variogram For Corrected Points")
###########################################

# PART G
model1 = fit.variogram(iso_var,vgm("Sph"))
plot(iso_var, model1 ,main="Empirical Variogram and Fitted Values")
summary(model1)

#############################################

#TODO: H

# How do you plot a fitted variogram along different azimuths?
plot(catch_t_log_resids.vario_dir,model1,main="Anisotropy-Corrected Variograms" )

#############################################

# PART I

poly = chull(x,y)
minx = min(x)
maxx = max(x)
miny = min(y)
maxy = max(y)
seqx = seq(minx,maxx,0.2)
seqy = seq(miny,maxy,0.2)
poly.in = polygrid(seqx,seqy, cbind(x,y)[poly,])
points(poly.in)
scall.rot <- data.frame(x,y,resid=catch_t_log_resids)
coordinates(scall.rot) = ~x+y
krige.out = krige(resid ~ 1, scall.rot ,poly.in, model=model1)

#####################################

# part J

pred_catch = exp(catch_t_log_trend + krige.out)

par(c(1,2))

















