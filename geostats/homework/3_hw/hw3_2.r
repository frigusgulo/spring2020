library(geoR)
data("soja98")

MO = soja98$MO
X = soja98$X
Y = soja98$Y

# EDA: Box plot, corellograms, heatmap, normal distribution, skewness
MO.conc <- interp(X,Y,MO,nx=150,ny=150)
par(c(1,2))
image(MO.conc,xlab="X",ylab="Y",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="MO concentrations",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
points(X,Y)
image.legend(120,80,zlim=range(MO),       # Puts a legend on the map with upper
             col=rev(heat.colors(24)))  

boxplot(MO,main="Percentage MO",ylab="Percentage")
hist(MO, xlab = "Percentage MO")

x = X
y = Y
u = MO
v = MO

h <- c(0,1,2,3,4)                              # Sets a vector of h-values.
outx <- matrix(nrow=length(h),ncol=3)
outy <- matrix(nrow=length(h),ncol=3)  
for (i in 0:6) outx[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(i,0)))                        #   calculates cross-functions for each.
cphx <- outx[,2];   # Defines the vectors of cross-covariances,
#   cross-corr's, & cross-semivariograms.
for (i in 0:6) outy[i+1,] <- as.numeric(            # Loops through distances of 0 to 6 and
  hscatter(x,y,v,u,c(0,i)))                        #   calculates cross-functions for each.
cphy <- outy[,2];   # Defines the vectors of cross-covariances,
#   cross-corr's,

par(mfrow=c(1,2))                                  # Sets up a 2x2 graphics window.
plot(h,cphx,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Residual Correlogram X axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)  

plot(h,cphy,type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Residual Correlogram Y axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)  