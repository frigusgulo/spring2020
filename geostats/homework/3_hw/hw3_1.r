library(gstat)
library(lattice) 
library(sp)
source("~/spring2020/geostats/scripts/allfunctions.r")
source("~/spring2020/geostats/scripts/movewin.r")
phytopthora = read.table("~/spring2020/geostats/datasets/phytoph.txt",header=T)
phy.mat = matrix(phytopthora$moist,nrow=20,byrow=T)
phy.mat[5,4] = NA
phy.mp = medpolish(phymat,na.rm=T)
phy.trend = phy.mat - phy.mp$residuals
par(mfrow=c(1,3)) 
all <- c(phy.mat,phy.trend)
zmin <- min(all)
zmax <- max(all[all!=max(all)])
image(phy.mat,zlim=range(all,na.rm=T),cex.axis=1.5,
      xlab="Columns",ylab="Rows",cex.lab=1.6,col=rev(gray.colors(24)))
image.legend(9.5,3,zlim=range(phy.mat,na.rm=T),col=rev(gray.colors(24))) 
title("Original Moisture %’s",cex.main=1.5)
image(phy.trend,zlim=range(phy.trend,na.rm=T),cex.axis=1.5,xlab="Columns",ylab="Rows",cex.lab=1.6,col=rev(gray.colors(24)))
image.legend(9.5,3,zlim=range(phy.trend,na.rm=T),col=rev(gray.colors(24)))
title("Median Polish Trend",cex.main=1.5)
image(phy.mp$resid,zlim=range(phy.mp$resid, na.rm=T),xlab="Columns",ylab="Rows",cex.lab=1.6,col=rev(gray.colors(24)))
image.legend(9.5,3,zlim=range(phy.mp$resid,na.rm=T),col=rev(gray.colors(24)))
title("Median Polish Residuals",cex.main=1.5)

#boxplot for row and column trends in the residual

phy.movewin = movewin(x=1:20,y=1:20,v=phy.mp$residuals,wx=3,wy=3,wox=1,woy=1)

x=1:20
y=1:20
v = phy.mp$resid

vario = variogram()

x=1:20
y=1:20
v = phy.mp$resid
u = phy.mp$resid
h <- c(0,1,2,3,4,5,6)                              # Sets a vector of h-values.
out <- matrix(nrow=length(h),ncol=3)               # Defines an hx3 blank matrix.
for (i in 0:6) out[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(i,0)))                        #   calculates cross-functions for each.
cch <- out[,1]; cph <- out[,2]; cgh <- out[,3]     # Defines the vectors of cross-covariances,
#   cross-corr's, & cross-semivariograms.
par(mfrow=c(2,2))                                  # Sets up a 2x2 graphics window.
plot(h,cch,type="n",axes=F,xlab="",ylab="")        # Creates a completely blank plot.
plot(h,cch,type="l",xlab="|h|",ylab="Cross - C(h)",# Plots the cross-correlations vs. h
     main="Residual Covariogram X axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h,cph,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Residual Correlogram X axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)                        
plot(h[-1],cgh[-1],type="l",xlab="|h|",cex.lab=1.5,# Plots the cross-semivariogram vs. h
     ylab="Cross - Gamma(h)",cex.axis=1.3,            #   with sizes controlled by "cex".
     main="Residual Semivariogram X axis",cex.main=1.5)