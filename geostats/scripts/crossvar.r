# Plots the cross h-scatterplots for the Walker Lake (n=100) V&U-data, and
# plots the cross-covariogram, cross-correlogram, and cross-semivariogram
# functions for the U,V data.  The hscatter function was used to generate
# the cross h-scatterplots.
# ========================================================================
walk100 <- read.table("Data/walk100.txt",header=T) # Reads in walk100 Walker Lake data.
x <- walk100$x                         # x is set to the x-values of "walk100".
y <- 11 - walk100$y                    # y is set to 11 - the y-values of "walk100".
v <- walk100$v                         # v is set to the V-values of "walk100".
u <- walk100$u                         # u is set to the U-values of "walk100".

# Creates the 5 cross-h scatterplots on page 35 of the class notes
# ================================================================
par(mfrow=c(3,2))                      # Sets up a 3x2 graphics window.
source("Scripts/hscatter.r")           # Loads and compiles hscatter function.
hscatter(x,y,v,u,c(0,0))               # Produces a scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,1))               # Produces a cross h=(0,1)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,2))               # Produces a cross h=(0,2)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(0,3))               # Produces a cross h=(0,3)-scatterplot of u vs. v.
hscatter(x,y,v,u,c(1,0))               # Produces a cross h=(1,0)-scatterplot of u vs. v.

# Plots the 3 cross-functions on page 38 of the class notes
# =========================================================
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
