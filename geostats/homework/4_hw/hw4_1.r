library(gstat)
library(lattice) 
library(sp)
library(geoR)
library(map)
library(usmap)
source("~/spring2020/geostats/scripts/allfunctions.r")
source("~/spring2020/geostats/scripts/rotateaxis.r")
source("~/spring2020/geostats/scripts/panel.gamma0.r")



scallops = read.table("~/spring2020/geostats/datasets/scallops.txt",header=T)
X = scallops$long
Y = scallops$lat
catch_t = scallops$tcatch
catch_t_log = log(catch_t)
summary(catch_t)
catch_t.conc = interp(X,Y,catch_t,nx=150,ny=150)
par(c(1,3))
image(catch_t.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(X,Y)
hist(catch_t,label="Total Catch")
hist(catch_t,label="Log Transformed Total Catch")

##########################

par(c(1,1))
usmap::plot_usmap(include = .new_england)
contour(x=X,y=Y,z=catch_t_log,main="Log Transformed Total Catch Contour")
points(X,Y)

############################

xy = cbin(X,Y)
coords_trans = rotateaxis(xy,52)
catch_t_loggamma = log(catch_t + 1)
newx = coords_trans$x
newy = coords_trans$y
nlfit = nls(catch_t_loggamma ~ (newx*delta + ((newy - min(newy))**(alpha-1))*exp(-(y - min(y))/beta + gamma) ) ,start=list(delta=0,alpha=2.2,beta=10,gamma=10))
summary(nlfit)
catch_t_log_trend = nlfit(newx,newy)
catch_t_log_resids = catch_t_log - catch_t_log_trend


catch_t_log.conc = interp(X,Y,catch_t_log,nx=150,ny=150)
catch_t_log_trend.conc = interp(newx,newy,catch_t_log_trend,nx=150,ny=150)
catch_t_log_resids.conc = interp(newx,newy,catch_t_log_resids,nx=150,ny=150)

par(c(1,3))
image(catch_t_log.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(X,Y)

image(catch_t_log_trend.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Inferred Trend of Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(newx,newy)

image(catch_t_log_resids.conc,xlab="Longitude",ylab="Lattitude",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="Inferred Trend Residuals of Log Transformed Total Catch",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(newx,newy)

#########################################

par(c(1,2))
plot(catch_t_log_trend,catch_t_log_resids,ylab='Residuals',xlab='Trend',main='Log Transformed Total Catch Fitted Trend vs Residuals')
qqnorm(catch_t_log_resids,main="Normal Q-Q Plot",xlab= "Theoretical Quantiles")
qqline(catch_t_log_resids,distribution = qnorm, probs = c(0.25, 0.75) )

##########################################
angles = seq(0,180,15)
df = data.frame(x=newx,y=newy,data=catch_t_log_resids)
catch_t_log.vario_dir = variogram(df,width=10,cutoff=100,alpha=angles,tol.hor = 15)
rose(catch_t_log.vario_dir,5)

catch_t_log.vario_omni = variogram(df,width=10,cutoff=100)
rose(catch_t_log.vario_omni,5)

#############################################
coords = cbin(newx,newy)
#TODO: F

# How do I use the xyplot function with panel.gamma0 to determine the direction of maximum spatial continuity?
# 
anisopairs = cbin(rho,theta)
coord_corrected = coords.aniso(coords,anisopairs)
df = data.frame(x=coords_corrected[,1],y=coords_corrected[,2],data=catch_t_log_resids)

par(c(1,2))
iso_var = variogram(df,width=10,cutoff=100)
rose(iso_car,5)
###########################################
fit = fit.variogram(catch_t_log.vario_omni,vgm("Sph"))
plot(catch_t_log.vario_omni, fit,main="Empirical Variogram and Fitted Values")
summary(fit)

#############################################

#TODO: H

# How do you plot a fitted variogram along different azimuths?

#############################################

poly = chull(newx,newy)
minx = min(newx)
maxx = max(newx)
miny = min(newy)
maxy = max(newy)
seqx = seq(minx,maxx,0.25)
seqy = seq(miny,maxy,0.25)
poly.in = polygrid(seqx,seqy, cbin(newx,newy)[poly,])
coordinates(poly.in) = ~x+y
krige.out = krige(catch_t_log_resids, poly.in, model=model1.out)

#####################################
















