library(geoR)
library(spatial)
library(fields)
source("~/spring2020/geostats/scripts/corrplot.r")
source("~/spring2020/geostats/scripts/hscatter.r")
landsat <- read.table("~/spring2020/geostats/homework/2_hw/landsateng.txt",header=T)
landsat_mat = as.matrix(landsat)
bandCovariance = cov(landsat_mat[,3:9])

for (i in 3:9){
  image.plot(interp(landsat_mat[,1],landsat_mat[,2],landsat_mat[,i],nx=100,ny=100))
}

###############################################
corrplot(x=landsat_mat[,1],y=landsat_mat[,2],landsat_mat[,4],landsat_mat[,4],
         h=c(1,0),5,dtol=0,atol=0)
corrplot(x=landsat_mat[,1],y=landsat_mat[,2],landsat_mat[,4],landsat_mat[,4],
         h=c(0,1),5,dtol=0,atol=0)
corrplot(x=landsat_mat[,1],y=landsat_mat[,2],landsat_mat[,6],landsat_mat[,6],
         h=c(1,0),5,dtol=0,atol=0)
corrplot(x=landsat_mat[,1],y=landsat_mat[,2],landsat_mat[,6],landsat_mat[,6],
         h=c(0,1),5,dtol=0,atol=0)

#################################################

# Plots the 3 cross-functions on page 38 of the class notes
# =========================================================]
x=landsat_mat[,1]
y=landsat_mat[,2]
v = landsat_mat[,4]
u = landsat_mat[,6]
h <- c(0,1,2,3,4,5,6)                              # Sets a vector of h-values.
out <- matrix(nrow=length(h),ncol=3)               # Defines an hx3 blank matrix.
for (i in 0:6) out[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(0,i)))                        #   calculates cross-functions for each.
cch <- out[,1]; cph <- out[,2]; cgh <- out[,3]     # Defines the vectors of cross-covariances,
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
  hscatter(x,y,v,u,c(i,0)))                        #   calculates cross-functions for each.
cch_y <- out_y[,1]; cph_y <- out_y[,2]; cgh_y <- out_y[,3]     # Defines the vectors of cross-covariances,
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