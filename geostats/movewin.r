movewin <- function(x=x,y=y,v=v,wx=wx,wy=wy,wox=wox,woy=woy){

# This function computes moving window means and variances for any type of data.
# Arguments to this function are:
#      x = a vector of the x-coordinates of the data
#      y = a vector of the y-coordinates of the data
#      v = a vector of the response variable values
#     wx = width of the window in the x-direction
#     wy = width of the window in the y-direction
#    wox = width of the x-overlap between adjacent windows
#    woy = width of the y-overlap between adjacent windows
# The function returns a data frame with new (x,y) coordinates corresponding to the
# midpoints of the windows, the number of points in each window (numvals), the
# moving window means (means), and the moving window standard deviations (sdevs).

  xedge <- min(x[x!=min(x)]-min(x))/2
  yedge <- min(y[y!=min(y)]-min(y))/2
  dimx <- max(x) - min(x) + 2*xedge
  dimy <- max(y) - min(y) + 2*yedge
  numx <- trunc((dimx-wx)/(wx-wox) + 1)
  numy <- trunc((dimy-wy)/(wy-woy) + 1)
  xmid <- matrix(rep(min(x)-xedge+(0:(numx-1))*(wx-wox),numy),nrow=numx)+(wx/2)
  ymid <- matrix(rep(min(y)-yedge+(0:(numy-1))*(wy-woy),numx),nrow=numx,byrow=T)+(wy/2)
  isin <- array(dim=c(numx,numy,length(x)))
  for (i in 1:length(x)) isin[,,i] <- (abs(x[i]-xmid)<=(wx/2))+(abs(y[i]-ymid)<=(wy/2))==2
  numvals <- apply(isin,c(1,2),sum)
  means <- matrix(nrow=numx,ncol=numy)
  sdevs <- matrix(nrow=numx,ncol=numy)
  for (i in 1:numx){
    for (j in 1:numy){
      if (numvals[i,j]>0){
        means[i,j] <- mean(v[isin[i,j,]==T])
        sdevs[i,j] <- sd(v[isin[i,j,]==T])
      }
    }
  }
  data.frame(x=c(xmid),y=c(ymid),numvals=c(numvals),means=c(means),sdevs=c(sdevs))
}
