library(geoR)
source("Scripts/corrplot.r")
source("~/spring2020/geostats/scripts/rose.r")

data("soja98")

MO = soja98$MO
X = soja98$X
Y = soja98$Y

# EDA: Box plot, corellograms, heatmap, normal distribution, skewness
MO.conc <- interp(X,Y,MO,nx=150,ny=150)
par(c(1,3))
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

h <- c(10,20,30,40,50)

outx <- corrplot(x,y,v,u,c(10,0),5,dtol=5,atol=15)
outy <- corrplot(x,y,v,u,c(0,10),5,dtol=5,atol=15)
##########################################################################################

par(mfrow=c(1,3))                                  # Sets up a 2x2 graphics window.
plot(h,outx[,3],type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Residual Correlogram X axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5)  

plot(h,outy[,3],type="l",xlab="|h|",ylab="Cross - p(h)",# Plots the cross-correlogram vs. h
     main="Residual Correlogram Y axis ",cex.lab=1.5,    #   with sizes controlled by "cex".
     cex.axis=1.3,cex.main=1.5) 


MO.mat <- matrix(MO,nrow=5,byrow=T) # Converts phytophthora data to a matrix.
MO.mp <- medpolish(MO.mat, na.rm=T)          # Performs median polish on the moistures.
MO.res <- MO.mp$residuals
image(MO.res,xlab="X",ylab="Y",cex.lab=1.6,   # Creates greyscale map of interpolated
      main="MO Residuals",   #   log concentrations using the colors
      cex.axis=1.5,col=rev(heat.colors(24)),        #   in "heat.colors" in reverse, with
      cex.main=1.2) #   axis labels and a title.
##########################################################################################

par(mfrow=c(1,2))  
MO.df = data.frame(x=x,y=y,data=MO)
MO.df_res = data.frame(x=x,y=y,data=MO.res)
angles = seq(0,180,15)
MO.vario = variogram(MO.df,width=10,cutoff=100,alpha=angles,tol.hor = 15)
rose(MO.vario,5)
MO.vario_res = variogram(MO.df_res,width=10,cutoff=100,alpha=angles,tol.hor = 15)
rose(MO.vario_res,5)