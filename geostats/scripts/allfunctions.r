# ========================================================================
# IMAGE.LEGEND
# This file contains all of my R-functions to date.  You can include it at
# the top of your code using the "source" function before each session, to
# avoid having to remember which functions you need to load prior to use.
# Some of these functions use specific libraries as well, so these will
# also need to be loaded separately.
# ========================================================================

image.legend <- 
    function(x,y, zlim, at.z = NULL, col = heat.colors(12), legnd=NULL, 
             lwd = max(3,32/length(col)), bg = NA, bty = "", ...) 
  ## * kein y.i -- Benutzer soll rein ueber lwd steuern; sollte reichen. 
  ## * legnd koennte interessant sein, falls Text geschrieben werden soll 
  ## (weiss mal wieder nicht, wie man aus legnd legend als option 
  ## macht) 
  ## * lwd wird per default in Abh. von col gewaehlt. 
{ 
    ## Purpose: 
    ## Authors: Martin Maechler, 9 Jul 2001 
    ## Martin Schlather, 24 Jul 2001 

  if (!is.null(legnd) && is.null(at.z)) 
      stop("at.z must be given if legnd is") ## falls legnd darf at.z 
    ## nicht automatisch gewaehlt werden 

    if(!is.numeric(zlim) || zlim[1] > zlim[2]) 
        stop("`zlim' must be numeric; zlim[1] <= zlim[2]") 
    if(is.null(at.z)) { 
        ## hier ein Versuch in Abhaengigkeit von n 
        ## die Anzahl der labels zu bestimmen: 
        n <- min(5, max(1,length(col)/10)) 
        at.z <- pretty(zlim,n=n,min.n=max(n %/% 3,1)) 

        ## es sieht nicht schoen aus, wenn pretty die letzte oder 
        ## erste zahl weit ausserhalb des zlim legt. 
        ## heuristisch nur 25% (oder so) ueberschreitung bzgl 
        ## intervalllaenge zulassen: 
        tol <- diff(at.z)[1] / 4 
        at.z <- at.z[(at.z>=zlim[1]-tol) & (at.z<=zlim[2]+tol)] 
      } 
    if(!is.numeric(at.z) || is.unsorted(at.z)) 
        stop("`at.z' must be numeric non-decreasing") 
    n.at <- length(at.z) 
    nc <- length(col) 
    if(n.at >= nc) 
        stop("length(at.z) must be (much) smaller than length(col)") 
    dz <- diff(zlim) 
    ## The colors must run equidistantly from zlim[1] to zlim[2]; 
    ## col[i] is for z-interval zlim[1] + [i-1, i) * dz/nc ; i = 1:nc 
    ## i.e., an at.z[] value z0 is color i0 = floor(nc * (z0 - zlim[1])/dz) 
    at.i <- floor(nc * (at.z - zlim[1])/dz ) 
    ## Possibly extend colors by `background' to the left and right 
    bgC <- if(is.null(bg)) NA else bg 
    if((xtra.l <- 1 - at.i[1]) > 0) { 
        at.i <- at.i + xtra.l 
        col <- c(rep(bgC, xtra.l), col) 
    } 
    if((xtra.r <- at.i[n.at] - nc) > 0) 
        col <- c(col, rep(bgC, xtra.r)) 
    lgd <- character(length(col)) 

    ## folgende if-Anweisung ist neu: 
    if (is.null(legnd)) lgd[at.i] <-format(at.z, dig = 3) 
    else { 
      if (length(legnd)!=length(at.z)) 
        stop("at.z and legnd must have the same length") 
      lgd[at.i] <- legnd 
    } 
    if((V <- R.version)$major <= 1 && V$minor <= 3.0 && V$status == "") 
{ 
        ## stop-gap fix around the bug that "NA" is not a valid color: 
        if(is.na(bgC)) { 
            lgd <- lgd[!is.na(col)] 
            col <- col[!is.na(col)] 
        } 
    } 
    legend(x,y, legend = rev(lgd), col = rev(col), 
           y.i = lwd/16, bty = bty, lwd = lwd, bg = bg, ...) 
} 

# ========================================================================
# PPNORM
# Defines a function to create normal PP-plots which takes four arguments:
#        x = the vector of data values,
#   points = the values at which normal quantiles are calculated,
#   xticks = the x-axis tick mark locations desired,
#   xlabel = the x-axis label.
# ========================================================================
ppnorm <- function(x,points,xticks,xlabel){
  repx <- matrix(rep(x,length(points)),nrow=length(points),byrow=T)
  reppt <- matrix(rep(points,length(x)),nrow=length(points),byrow=F)
  diff <- reppt-repx
  cumfreq <- apply(diff,1,function(a) mean(a>=0))
  zq <- qnorm(cumfreq,0,1)
  minx <- 2*points[1]-points[2]
  plot(points,zq,xlab="",ylab="",xlim=c(minx,
    2*points[length(points)]-points[length(points)-1]),ylim=c(-2.5,2.5),
    pch="+",axes=F,cex=1.5,type="b")
  axis(1,at=xticks,tck=.02,pos=-2.5,mgp=c(3,.5,0),cex.axis=1.5)
  axis(2,at=qnorm(c(1,2,5,10,20,50,80,90,95,98,99)/100,0,1),tck=.02,
    pos=minx,mgp=c(3,.5,0),labels=c("1","2","5","10","20","50","80","90",
    "95","98","99"),las=2,cex.axis=1.3)
  mtext(side=1,line=2.3,xlabel,cex=1.5)
  mtext(side=2,line=2.0,"Cumulative Frequency (%)",cex=1.5)
  i25 <- sum(cumfreq<.25)
  end25 <- points[i25] + (points[i25+1]-points[i25])*((.25-cumfreq[i25])/
    (cumfreq[i25+1]-cumfreq[i25]))
  lines(c(minx,end25),c(qnorm(.25,0,1),qnorm(.25,0,1)),lty=8)
  lines(c(end25,end25),c(-2.5,qnorm(.25,0,1)),lty=8)
  text(end25,-2.62,"Q1",cex=1.2)
  i75 <- sum(cumfreq<.75)
  end75 <- points[i75] + (points[i75+1]-points[i75])*((.75-cumfreq[i75])/
    (cumfreq[i75+1]-cumfreq[i75]))
  lines(c(minx,end75),c(qnorm(.75,0,1),qnorm(.75,0,1)),lty=8)
  lines(c(end75,end75),c(-2.5,qnorm(.75,0,1)),lty=8)
  text(end75,-2.62,"Q3",cex=1.2)
}

# =================================================================================
# MOVEWIN
# This function computes moving window means and variances for any type of data.
# Arguments to this function are:
#      x = a vector of the x-coordinates of the data
#      y = a vector of the y-coordinates of the data
#      v = a vector of the response variable values
#     wx = width of the window in the x-direction
#     wy = width of the window in the y-direction
#     wo = width of the overlap between adjacent windows
# The function returns a data frame with new (x,y) coordinates corresponding to the
# midpoints of the windows, the number of points in each window (numvals), the
# moving window means (means), and the moving window standard deviations (sdevs).
# =================================================================================
movewin <- function(x=x,y=y,v=v,wx=wx,wy=wy,wo=wo){
  edge <- min(x[x!=min(x)]-min(x))/2
  dimx <- max(x) - min(x) + 2*edge
  dimy <- max(y) - min(y) + 2*edge
  numx <- trunc((dimx-wx)/(wx-wo) + 1)
  numy <- trunc((dimy-wy)/(wy-wo) + 1)
  xmid <- matrix(rep(min(x)-edge+(0:(numx-1))*(wx-wo),numy),nrow=numx)+(wx/2)
  ymid <- matrix(rep(min(y)-edge+(0:(numy-1))*(wy-wo),numx),nrow=numx,byrow=T)+(wy/2)
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

# ===============================================================================
# HSCATTER
# This function generates h-scatterplots or cross h-scatterplots in any direction
# desired, for use with lattice and geostatistical data.  Arguments to this
# function are:
#     x = a vector of the x-coordinates of the grid
#     y = a vector of the y-coordinates of the grid
#     u = a vector of the response variable U values
#     v = a vector of the response variable V values
#     h = a 2x1 vector indicating the direction h
#  dtol = a distance tolerance (used with non-lattice data)
#  atol = an angle tolerance (used with non-lattice data)
# The function creates the h-scatterplot and returns the covariance, correlation,
# and semivariogram of the paired values.  If an h-scatterplot is desired, u & v
# should be the same variable, and if a cross h-scatterplot is desired, u & v
# should be different.
# ================================================================================
hscatter <- function(x,y,u,v,h,dtol=NA,atol=NA){
 if (missing(dtol) && missing(atol)){
  umat <- tapply(u,list(y,x),function(z)z)
  vmat <- tapply(v,list(y,x),function(z)z)
  if (h[1] > 0) vmat <- vmat[,-((max(x)+1-h[1]):max(x))]
  if (h[1] > 0) umat <- umat[,-(1:h[1])]
  if (h[2] > 0) vmat <- vmat[-((max(y)+1-h[2]):max(y)),]
  if (h[2] > 0) umat <- umat[-(1:h[2]),]
  if (h[1] < 0) vmat <- vmat[,-(1:-h[1])]
  if (h[1] < 0) umat <- umat[,-((max(x)+1+h[1]):max(x))]
  if (h[2] < 0) vmat <- vmat[-(1:-h[2]),]
  if (h[2] < 0) umat <- umat[-((max(y)+1+h[2]):max(y)),]
  umat1 <- tapply(u,list(y,x),function(z)z)
  vmat1 <- tapply(v,list(y,x),function(z)z)
  if (h[1] > 0) vmat1 <- vmat1[,-(1:h[1])]
  if (h[1] > 0) umat1 <- umat1[,-((max(x)+1-h[1]):max(x))]
  if (h[2] > 0) vmat1 <- vmat1[-(1:h[2]),]
  if (h[2] > 0) umat1 <- umat1[-((max(y)+1-h[2]):max(y)),]
  if (h[1] < 0) vmat1 <- vmat1[,-((max(x)+1+h[1]):max(x))]
  if (h[1] < 0) umat1 <- umat1[,-(1:-h[1])]
  if (h[2] < 0) vmat1 <- vmat1[-((max(y)+1+h[2]):max(y)),]
  if (h[2] < 0) umat1 <- umat1[-(1:-h[2]),]
  semivar <- (1/(2*length(c(vmat))))*sum((c(umat) -
             c(umat1))*(c(vmat1) - c(vmat)))
 }
 else{
  n <- length(u)
  dmat <- matrix(0,nrow=n,ncol=n)
  xy <- cbind(x,y)
  dmat[lower.tri(dmat)] <- dist(xy)
  dmat <- dmat + t(dmat)
  xdiff <- matrix(rep(x,n),nrow=n,byrow=T) - matrix(rep(x,n),nrow=n)
  ydiff <- matrix(rep(y,n),nrow=n,byrow=T) - matrix(rep(y,n),nrow=n)
  angle <- atan(ydiff/xdiff)*180/pi
  disth <- sqrt(sum(h*h))
  angleh <- atan(h[2]/h[1])*180/pi
  index <- seq(1,n*n)
  dmat <- c(dmat); angle <- c(angle)
  iloc <- index[(abs(dmat-disth)<dtol)&(abs(angle-angleh)<atol)]
  irow <- trunc(((iloc-1)/n)+1)
  icol <- iloc - trunc((iloc-1)/n)*n
  vmat <- v[icol]; umat <- u[irow]
  semivar <- (1/(2*length(vmat)))*sum((u[irow]-u[icol])*(v[irow]-v[icol]))
 }
 xlims <- c(min(u),max(u))
 ylims <- c(min(v),max(v))
 lab1 <- paste(deparse(substitute(u)),"(t)",sep="")
 lab2 <- paste(deparse(substitute(v)),"(t+h)",sep="")
 plot(umat,vmat,pch="+",xlab="",ylab="",xlim=xlims,ylim=ylims,axes=F,cex=1.2)
 axis(1,tck=.02,pos=min(u,v),mgp=c(3,.5,0),cex.axis=1.3)
 axis(2,tck=.02,pos=min(u,v),mgp=c(3,.5,0),cex.axis=1.3)
 mtext(side=1,line=2.0,lab1,cex=1.3)
 mtext(side=2,line=2.3,lab2,cex=1.3)
 mtext(side=3,line=0,paste("h = (",h[1],",",h[2],")",sep=""),cex=1.3)
 lines(c(min(u,v),max(u,v)),c(min(u,v),max(u,v)),lty=8)
 numpairs <- length(c(vmat))
 cat("\n Correlation:   p(h)=", cor(c(vmat),c(umat)))
 cat("\n Covariance:    C(h)=", cor(c(vmat),c(umat))*sqrt(var(c(umat))*
                              var(c(vmat)))*(numpairs-1)/numpairs)
 cat("\n Semivariogram: g(h)=", semivar)
 cat("\n")
 data.frame(covar=cor(c(vmat),c(umat))*sqrt(var(c(umat))*var(c(vmat)))*
   (numpairs-1)/numpairs,corr=cor(c(vmat),c(umat)),semivar=semivar)
}

corrplot <- function(x,y,u,v,h,numlags,dtol=NA,atol=NA){

# This function generates the covariogram, correlogram, and semivariogram plots in the
# direction specified by "h", and with "numlags" lags in that direction.  Additionally,
# the values of these functions are returned as a data frame.  Arguments to this
# function are:
#         x = a vector of the x-coordinates of the grid
#         y = a vector of the y-coordinates of the grid
#         u = a vector of the response variable U values
#         v = a vector of the response variable V values
#         h = a 2x1 vector indicating the direction h
#   numlags = the number of lags at which to determine the functions
#      dtol = a distance tolerance for non-lattice data (optional)
#      atol = an angle tolerance for non-lattice data (optional)
# If the cross-covariogram, cross-correlogram, and cross-semivariogram are desired,
# then u & v should be different.  Otherwise, u & v should be the same.

correlog <- rep(0,numlags)
covariog <- rep(0,numlags)
semivar <- rep(0,numlags)

for (i in 1:numlags){
  h2 <- h*i
  if (missing(dtol)) dtol <- max(abs(h2))*.0001
  if (missing(atol)) atol <- .0001
  n <- length(u)
  dmat <- matrix(0,nrow=n,ncol=n)
  xy <- cbind(x,y)
  dmat[lower.tri(dmat)] <- dist(xy)
  dmat <- dmat + t(dmat)
  xdiff <- matrix(rep(x,n),nrow=n,byrow=T) - matrix(rep(x,n),nrow=n)
  ydiff <- matrix(rep(y,n),nrow=n,byrow=T) - matrix(rep(y,n),nrow=n)
  angle <- atan(ydiff/xdiff)*180/pi
  disth <- sqrt(sum(h2*h2))
  angleh <- atan(h2[2]/h2[1])*180/pi
  index <- seq(1,n*n)
  dmat <- c(dmat); angle <- c(angle)
  iloc <- index[(abs(dmat-disth)<=dtol)&(abs(angle-angleh)<=atol)]
  irow <- trunc(((iloc-1)/n)+1)
  icol <- iloc - trunc((iloc-1)/n)*n
  if (h[2]==0){
    keep <- irow<icol
    irow <- irow[keep]; icol <- icol[keep]
  }
  vmat <- v[icol]; umat <- u[irow]
  nh <- sum(!is.na((u[irow]-u[icol])*(v[irow]-v[icol])))
  semivar[i] <- (1/(2*nh))*sum((u[irow]-u[icol])*(v[irow]-v[icol]),na.rm=T)
  numpairs <- min(sum(!is.na(vmat)),sum(!is.na(umat)))
  umatvec <- c(umat); vmatvec <- c(vmat)
  include <- !is.na(vmatvec)&!is.na(umatvec)
  correlog[i] <- cor(vmatvec[include],umatvec[include])
  covariog[i] <- cor(vmatvec[include],umatvec[include])*
    sqrt(var(umatvec,na.rm=T)*var(vmatvec,na.rm=T))*(numpairs-1)/numpairs
}
include <- !is.na(v)&!is.na(u)
covar0 <- cor(v[include],u[include])*sqrt(var(u,na.rm=T)*var(v,na.rm=T))*
  (numpairs-1)/numpairs
corr0 <- cor(v[include],u[include])
if (h[1]==0 && h[2]>0) azimuth <- 0
if (h[1]==0 && h[2]<0) azimuth <- 180
if (h[1]!=0) azimuth <- 90-atan(h[2]/h[1])*180/pi

par(mfrow=c(2,2))
xlims <- c(0,numlags); ylims <- c(min(covariog,covar0,0),max(covariog,covar0))
plot(0:numlags,c(covar0,covariog),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.0,"Covariogram C(h)",cex=1.0)
xlims <- c(0,numlags); ylims <- c(min(correlog,corr0,0),max(correlog,corr0))
plot(0:numlags,c(corr0,correlog),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.7,"Correlogram p(h)",cex=1.0)
xlims <- c(0,numlags); ylims <- c(0,max(semivar))
plot(0:numlags,c(0,semivar),type="b",xlab="",ylab="",xlim=xlims,ylim=ylims,
  axes=F,cex=1.0)
axis(1,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
axis(2,tck=.02,pos=0,mgp=c(3,.5,0),cex.axis=0.9)
mtext(side=1,line=2.0,"Lag (h)",cex=1.0)
mtext(side=2,line=2.7,"Semivariogram Gamma(h)",cex=1.0)
plot(c(0,1),c(0,1),xlab="",ylab="",axes=F,pch=" ")
text(.5,.6,paste("Azimuth = ",round(azimuth,1),sep=""),cex=1.5)

data.frame(lag=1:numlags,corr=correlog,covar=covariog,semivar=semivar)
}

# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle, and returns the 3 interior angles of the triangle.
# =========================================================================================
angles.tri <- function(x,y){
  pi <- 3.14159
  ord <- order(x)
  y <- y[ord]; x <- sort(x)
  ratdiff1 <- (y[2]-y[1])/(x[2]-x[1])
  ratdiff2 <- (y[3]-y[1])/(x[3]-x[1])
  ratdiff3 <- (y[3]-y[2])/(x[3]-x[2])
  theta <- rep(0,3)
  theta[1] <- abs(atan(ratdiff1)-atan(ratdiff2))*180/pi
  theta[3] <- abs(atan(ratdiff2)-atan(ratdiff3))*180/pi
  theta[2] <- 180-abs(atan(ratdiff1)-atan(ratdiff3))*180/pi
  theta[order(ord)]
}

# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle, and returns the unique point which is equidistant to all
# three points of the triangle.
# =========================================================================================

mid.tri <- function(x,y){
  y <- y[order(x)]; x <- sort(x)
  m1 <- (x[1]-x[2])/(y[2]-y[1])
  m2 <- (x[1]-x[3])/(y[3]-y[1])
  b1 <- ((y[2]^2-y[1]^2)-(x[1]^2-x[2]^2))/(2*(y[2]-y[1]))
  b2 <- ((y[3]^2-y[1]^2)-(x[1]^2-x[3]^2))/(2*(y[3]-y[1]))
  if (y[1]==y[2]) mid <- c((x[1]+x[2])/2,m2*(x[1]+x[2])/2+b2)
  if (y[1]==y[3]) mid <- c((x[1]+x[3])/2,m1*(x[1]+x[3])/2+b1)
  if ((y[1]!=y[2]) && (y[1]!=y[3]))
    mid <- solve(matrix(c(-m1,1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1,b2),nrow=2)
  mid
}

# =========================================================================================
# This function takes as arguments x- and y-coordinates as vectors of length 3 which define
# the coordinates of a triangle on the edge of the region, and returns the two intersection
# nodes from the equidistant point outside the region to the interior.
# =========================================================================================
edgenodes.tri <- function(x,y,mid,an){
  pi <- 3.14159
  obpt <- order(an)[3]
  acpts <- sort(order(an)[1:2])
  xmids <- (x[acpts]+x[obpt])/2
  ymids <- (y[acpts]+y[obpt])/2
  m1 <- (mid[2]-ymids)/(mid[1]-xmids)
  m2 <- (y[acpts[2]]-y[acpts[1]])/(x[acpts[2]]-x[acpts[1]])
  b1 <- mid[2] - m1*mid[1]
  b2 <- y[acpts[1]] - m2*x[acpts[1]]
  if (x[acpts[1]]==x[acpts[2]]){
     node1 <- c(x[acpts[1]],m1[1]*x[acpts[1]]+b1[1])
     node2 <- c(x[acpts[1]],m1[2]*x[acpts[1]]+b1[2])
  }
  if (xmids[1]==mid[1]) node1 <- c(mid[1],m2*mid[1]+b2)
  if (abs(xmids[2]-mid[1])<10^(-10)) node2 <- c(mid[1],m2*mid[1]+b2)
  if ((x[acpts[1]]!=x[acpts[2]]) && (xmids[1]!=mid[1])) node1 <-
    solve(matrix(c(-m1[1],1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1[1],b2),nrow=2)
  if ((x[acpts[1]]!=x[acpts[2]]) && (abs(xmids[2]-mid[1])>=10^(-10))) node2 <-
    solve(matrix(c(-m1[2],1,-m2,1),nrow=2,byrow=T))%*%matrix(c(b1[2],b2),nrow=2)
  matrix(c(node1,node2),nrow=2)
}

# =========================================================================================
# The "polydec" function takes a vectors of x- and y-coordinates and finds the influence
# polygons for each point, returning the areas of the polygons for use with the polygonal
# declustering global estimation method.  An outer edge is placed around the region at a
# distance roughly equal to half the minimum distance of a site to its neighbors, to handle
# edge effects.  To help place this boundary correctly, there is an optional "peels"
# argument that can be given.  The default is 1 meaning a single convex hull is taken.
# Generally, "peels" will not need to be any higher than 3 or 4.  This function calls the
# three functions in the "triangfuncs.r" file, so that this file must be compiled before
# executing "polydec".  Typing: "pd <- polydec(x,y,1)" will compute the influence polygons
# for all points in "x" and "y", returning the areas of these polygons in "pd".
# =========================================================================================
polydec <- function(x,y,peels=NA){
  if (missing(peels)) peels <- 1
  alldist <- dist(cbind(x,y))
  distmat <- matrix(0,nrow=attr(alldist,"Size"),ncol=attr(alldist,"Size"))
  distmat[lower.tri(distmat)] <- alldist
  distmat <- distmat + t(distmat)
  edgedist <- max(apply(distmat,1,function(a) a[order(a)[2]]))
  delaun <- triangles(tri.mesh(x,y))
  at <- t(delaun[,1:3])
  tx <- x; ty <- y
  ch <- rep(0,0)
  for (i in 1:peels){
    ch <- c(ch,chull(tx,ty))
    tx[ch] <- mean(tx)
    ty[ch] <- mean(ty)
  }
  atmat <- matrix(c(rep(at[1,],length(ch)),rep(at[2,],length(ch)),rep(at[3,],length(ch))),
    ncol=length(at[1,]),byrow=T)
  chmat <- matrix(rep(matrix(rep(ch,length(at[1,])),nrow=length(ch),byrow=F),3),ncol=length(at[1,]),
    byrow=F)
  equal <- atmat==chmat
  keep <- apply(equal,2,function(a) sum(a)>0)
  at <- at[,keep]
  pairs <- cbind(at[1:2,],at[c(1,3),],at[2:3,])
  extra.lab <- matrix(rep(pairs,6),nrow=2,byrow=F)
  pairdist <- sqrt((x[pairs[1,]]-x[pairs[2,]])^2+(y[pairs[1,]]-y[pairs[2,]])^2)
  pairs <- rbind(pairs,pairdist)
  end1.x <- (x[pairs[1,]]+x[pairs[2,]])/2 + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end2.x <- (x[pairs[1,]]+x[pairs[2,]])/2 - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end3.x <- x[pairs[1,]] + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end4.x <- x[pairs[2,]] - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end5.x <- x[pairs[1,]] - (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end6.x <- x[pairs[2,]] + (edgedist/pairs[3,])*(y[pairs[2,]]-y[pairs[1,]])/2
  end1.y <- (y[pairs[1,]]+y[pairs[2,]])/2 - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end2.y <- (y[pairs[1,]]+y[pairs[2,]])/2 + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end3.y <- y[pairs[1,]] - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end4.y <- y[pairs[2,]] + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end5.y <- y[pairs[1,]] + (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  end6.y <- y[pairs[2,]] - (edgedist/pairs[3,])*(x[pairs[2,]]-x[pairs[1,]])/2
  extra <- rbind(c(end1.x,end2.x,end3.x,end4.x,end5.x,end6.x),c(end1.y,end2.y,end3.y,end4.y,end5.y,end6.y))
  extra.out <- point.in.polygon(extra[1,],extra[2,],x[chull(x,y)],y[chull(x,y)])
  extra.x <- c(matrix(rep(extra[1,extra.out==F],2),nrow=2,byrow=F))
  extra.y <- c(matrix(rep(extra[2,extra.out==F],2),nrow=2,byrow=F))
  extra.lab <- c(t(extra.lab[,extra.out==F]))
  extra.nodes <- rbind(extra.x,extra.y,extra.lab)
  tr <- t(delaun[,1:3])
  an <- apply(tr,2,function(a,x,y) angles.tri(x[a],y[a]),x,y)

  maxan <- apply(an,2,max)
  an <- an[,maxan<180]
  tr <- tr[,maxan<180]
  acute <- apply(an,2,function(a) max(a)<90)
  tr.acute <- tr[,acute]
  mid.acute <- apply(tr.acute,2,function(a,x,y) mid.tri(x[a],y[a]),x,y)
  midseg1 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[1:2]]),
    mean(y[a[1:2]])),nrow=2,byrow=T),x,y)
  midseg2 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[c(1,3)]]),
    mean(y[a[c(1,3)]])),nrow=2,byrow=T),x,y)
  midseg3 <- apply(tr.acute,2,function(a,x,y) matrix(c(mean(x[a[2:3]]),
    mean(y[a[2:3]])),nrow=2,byrow=T),x,y)
  segs.ac <- cbind(rbind(x,y),matrix(rep(mid.acute,3),nrow=2),midseg1,
    midseg1,midseg2,midseg2,midseg3,midseg3)
  segs.ac <- rbind(segs.ac,c(1:length(x),c(t(tr.acute)),tr.acute[1,],tr.acute[2,],
    tr.acute[1,],tr.acute[3,],tr.acute[2,],tr.acute[3,]))

  an.ob <- an[,acute==F]
  tr.obtuse <- tr[,acute==F]
  mid.obtuse <- apply(tr.obtuse,2,function(a,x,y) mid.tri(x[a],y[a]),x,y)
  ext.nodes <- matrix(nrow=3,ncol=2*length(tr.obtuse[1,]))
  for (i in 1:length(tr.obtuse[1,])) ext.nodes[1:2,(2*i-1):(2*i)] <-
    edgenodes.tri(x[tr.obtuse[,i]],y[tr.obtuse[,i]],mid.obtuse[,i],an.ob[,i])
  ext.nodes <- cbind(ext.nodes,ext.nodes)
  obpoints0 <- apply(an.ob,2,function(a) order(a)[3])
  obpoints <- diag(apply(tr.obtuse,2,function(a,obpoints0) a[obpoints0],obpoints0))
  obpoints <- c(matrix(rep(obpoints,2),nrow=2,byrow=T))
  others1 <- rep(1,length(obpoints0))
  others1[obpoints0==1] <- 2
  others2 <- rep(3,length(obpoints0))
  others2[obpoints0==3] <- 2
  others1 <- diag(apply(tr.obtuse,2,function(a,others1) a[others1],others1))
  others2 <- diag(apply(tr.obtuse,2,function(a,others2) a[others2],others2))
  others <- c(matrix(c(others1,others2),nrow=2,byrow=T))
  ext.nodes[3,] <- c(obpoints,others)

  ob.in <- point.in.polygon(mid.obtuse[1,],mid.obtuse[2,],
    extra.x[chull(extra.x,extra.y)],extra.y[chull(extra.x,extra.y)])
  ob.x <- c(matrix(rep(mid.obtuse[1,ob.in==T],3),nrow=3,byrow=T))
  ob.y <- c(matrix(rep(mid.obtuse[2,ob.in==T],3),nrow=3,byrow=T))
  ob.lab <- c(tr.obtuse[,ob.in==T])
  ob.nodes <- rbind(ob.x,ob.y,ob.lab)

  mids.x <- c(matrix(rep((x[obpoints]+x[others])/2,2),nrow=2,byrow=T))
  mids.y <- c(matrix(rep((y[obpoints]+y[others])/2,2),nrow=2,byrow=T))
  mids.lab <- c(matrix(c(obpoints,others),nrow=2,byrow=T))
  mids.nodes <- rbind(mids.x,mids.y,mids.lab)

  allnodes <- cbind(segs.ac,ext.nodes,ob.nodes,mids.nodes,extra.nodes)
  dist.nodes <- sqrt((allnodes[1,]-x[allnodes[3,]])^2+(allnodes[2,]-y[allnodes[3,]])^2)
  min.dist <- apply(allnodes,2,function(a,x,y) min(sqrt((a[1]-x)^2 + (a[2]-y)^2)),x,y)
  minxydist <- min(dist(cbind(x,y)))
  omit.node <- abs(dist.nodes-min.dist)>(.00001*minxydist)
  allnodes <- allnodes[,omit.node==F]

  par(pty="s")
  xlims <- c(min(allnodes[1,]),max(allnodes[1,]))
  ylims <- c(min(allnodes[2,]),max(allnodes[2,]))
  plot(x,y,type="p",axes=T,xlab="",ylab="",xlim=xlims,ylim=ylims)
  polyarea <- rep(0,length(x))
  for (i in 1:length(x)){
    nodes <- allnodes[,allnodes[3,]==i]
    xord <- order(nodes[1,])
    polygon(nodes[1,chull(nodes[1,],nodes[2,])],nodes[2,chull(nodes[1,],nodes[2,])],
      density=0)
    polyarea[i] <- areapl(cbind(nodes[1,chull(nodes[1,],nodes[2,])],
      nodes[2,chull(nodes[1,],nodes[2,])]))
  }
  polyarea
}

# The function cell.dec performs cell declustering to estimate the
# global mean over some region.  There are 5 arguments to this
# function, as outlined below, and the function returns the global
# mean estimate.  The five arguments are:
#       x = the vector of x-coordinates of sample values,
#       y = the vector of y-coordinates of sample values,
#       v = the vector of response values at the sample locations,
#  xcells = the number of cells desired in the x-direction,
#  ycells = the number of cells desired in the y-direction.
#
# As an example, if we wanted to perform cell declustering on the
# Walker Lake data with 3 cells in the x-direction, and 5 cells in
# the y-direction, we would type:
#    > attach(walk470)
#    > globmean <- celldec(x,y,v,3,5)

celldec <- function(x,y,v,xcells,ycells){
  leftlim <- min(x)-(x[order(x)[2]]-min(x))/2
  rightlim <- max(x)+(max(x)-x[order(x)[length(x)-1]])/2
  downlim <- min(y)-(y[order(y)[2]]-min(y))/2
  uplim <- max(y)+(max(y)-y[order(y)[length(y)-1]])/2
  xgap <- (rightlim-leftlim)/xcells
  ygap <- (uplim-downlim)/ycells
  xmids <- (leftlim+rightlim)/2; ymids <- (downlim+uplim)/2
  if (xcells!=1) xmids <- seq(leftlim+xgap/2,rightlim-xgap/2,xgap)
  if (ycells!=1) ymids <- seq(downlim+ygap/2,uplim-ygap/2,ygap)
  xmat <- t(matrix(rep(x,length(xmids)),nrow=length(x)))
  ymat <- t(matrix(rep(y,length(ymids)),nrow=length(y)))
  xmids <- matrix(rep(xmids,length(x)),nrow=length(xmids))
  ymids <- matrix(rep(ymids,length(y)),nrow=length(ymids))
  distx <- matrix(abs(xmat-xmids),nrow=xcells)
  disty <- matrix(abs(ymat-ymids),nrow=ycells)
  xcelllab <- matrix(apply(distx,2,order),nrow=xcells)[1,]
  ycelllab <- matrix(apply(disty,2,order),nrow=ycells)[1,]
  celllab <- (ycelllab-1)*xcells + xcelllab
  d <- length(unique(celllab))
  numrep <- rep(0,length(x))
  for (i in 1:length(x)) numrep[i] <- sum(celllab==celllab[i])
  cdmean <- sum((1/d)*(v/numrep))
}

# The function point.crossval takes as arguments a vector of x-coordinates (x),
# a vector of y-coordinates (y), a vector of response values (z), a power for
# inverse distance weighting (p) and a radius of influence for both inverse
# distance weighting and the unweighted local sample mean (r).  It performs
# crossvalidation on the observed values using the point estimation methods:
# polygonal declustering, triangulation, unweighted local sample mean, and
# inverse distance weighting.  Tables of output for each of the four methods
# are produced with means, standard deviations, 5-number summaries, MSE, MAE,
# and the correlation between the true and estimated values.  Also returned
# with the function under the names pred.pd, pred.tr, pred.ls, pred.iw are the
# vectors of predicted values corresponding to the true values on which
# crossvalidation was performed.

# Libraries requried:  tripack, splancs

point.crossval <- function(x,y,z,p,r){
  n <- length(x); pred.pd <- rep(0,n); pred.iw <- rep(0,n)
  pred.tr <- rep(NA,n); pred.ls <- rep(0,n)
  for (i in 1:n){
    cat("i = ",i,"\n")
    dist0 <- sqrt((x-x[i])**2 + (y-y[i])**2)
    pred.pd[i] <- z[order(dist0)[2]]
    local <- seq(1,length(dist0))[-c(seq(1,length(dist0))[dist0>r],i)]
    if (length(local)<1) stop("Radius of Influence Too Small")
    pred.ls[i] <- mean(z[local])
    weights <- (1/dist0[local]**p)/sum(1/dist0[local]**p)
    pred.iw[i] <- sum(z[local]*weights)
    if (!is.element(i,chull(x,y))){
      delaun <- triangles(tri.mesh(x[-i],y[-i]))
      numtriang <- length(delaun[,1])
      area <- rep(0,numtriang)
      x2 <- x[-i]; y2 <- y[-i]; z2 <- z[-i]
      intri <- rep(F,numtriang)
      for (j in 1:numtriang)
        intri[j] <- in.convex.hull(tri.mesh(x2[delaun[j,1:3]],
                                   y2[delaun[j,1:3]]),x[i],y[i])
      tri <- rev(order(intri==T))[1]
      A1 <- areapl(cbind(c(x2[delaun[tri,1:2]],x[i]),c(y2[delaun[tri,1:2]],y[i])))
      A2 <- areapl(cbind(c(x2[delaun[tri,c(1,3)]],x[i]),c(y2[delaun[tri,c(1,3)]],y[i])))
      A3 <- areapl(cbind(c(x2[delaun[tri,2:3]],x[i]),c(y2[delaun[tri,2:3]],y[i])))
      A <- A1+A2+A3
      pred.tr[i] <- (A1*z[delaun[tri,3]] + A2*z[delaun[tri,2]] + A3*z[delaun[tri,1]])/A
    }
  }
  nt <- n - sum(is.na(pred.tr))
  nonna <- order(pred.tr)[1:nt]
  mse.pd <- (1/n)*sum((pred.pd-z)**2)
  mse.ls <- (1/n)*sum((pred.ls-z)**2)
  mse.iw <- (1/n)*sum((pred.iw-z)**2)
  mse.tr <- (1/nt)*sum((pred.tr[nonna]-z[nonna])**2)
  mae.pd <- (1/n)*sum(abs(pred.pd-z))
  mae.ls <- (1/n)*sum(abs(pred.ls-z))
  mae.iw <- (1/n)*sum(abs(pred.iw-z))
  mae.tr <- (1/nt)*sum(abs(pred.tr[nonna]-z[nonna]))
  corr.pd <- cor(pred.pd,z)
  corr.ls <- cor(pred.ls,z)
  corr.iw <- cor(pred.iw,z)
  if (nt>1) corr.tr <- cor(pred.tr[nonna],z[nonna])
  if (nt<=1) corr.tr <- NA
  numdec <- trunc(log10(mae.pd)-1)
  cat("Polygonal Declustering","\n======================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",round(mean(
    pred.pd),numdec),"   SD = ",round(sqrt(var(pred.pd)),numdec),"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.pd),
    numdec),")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.pd-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.pd-z),numdec),
    round(sqrt(var(pred.pd-z)),numdec),")\n")
  cat("MSE = ",round(mse.pd,numdec),"   MAE = ",round(mae.pd,numdec),
      "   Correlation = ",round(corr.pd,3),"\n\n")
  cat("Delaunay Triangulation","\n======================\n")
  cat("n = ",nt,"   True Mean = ",round(mean(z[nonna]),numdec),"   m = ",round(mean
    (pred.tr[nonna]),numdec),"   SD = ",round(sqrt(var(pred.tr[nonna])),
    numdec),"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.tr[nonna]),
    numdec),")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.tr[nonna]-
    z[nonna]),numdec),")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.tr[nonna]-z[nonna]),numdec),
    round(sqrt(var(pred.tr[nonna]-z[nonna])),numdec),")\n")
  cat("MSE = ",round(mse.tr,numdec),"   MAE = ",round(mae.tr,numdec),
      "   Correlation = ",round(corr.tr,3),"\n\n")
  cat("Local Sample Mean (Unweighted)","\n==============================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",
    round(mean(pred.ls),numdec),"   SD = ",round(sqrt(var(pred.ls)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.ls),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.ls-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.ls-z),numdec),
    round(sqrt(var(pred.ls-z)),numdec),")\n")
  cat("MSE = ",round(mse.ls,numdec),"   MAE = ",round(mae.ls,numdec),
      "   Correlation = ",round(corr.ls,3),"\n\n")
  cat("Inverse Distance Weighting","\n==========================\n")
  cat("n = ",n,"   True Mean = ",round(mean(z),numdec),"   m = ",
    round(mean(pred.iw),numdec),"   SD = ",round(sqrt(var(pred.iw)),numdec),
    "   p = ",p,"   r = ",r,"\n")
  cat("5-Number Summary (True Values): (",round(quantile(z),numdec),")\n")
  cat("5-Number Summary (Estimated Values): (",round(quantile(pred.iw),numdec),
    ")\n")
  cat("5-Number Summary (Residuals): (",round(quantile(pred.iw-z),numdec),
    ")\n")
  cat("Residual (Mean,SD) = (",round(mean(pred.iw-z),numdec),
    round(sqrt(var(pred.iw-z)),numdec),")\n")
  cat("MSE = ",round(mse.iw,(numdec+1)),"   MAE = ",round(mae.iw,(numdec+1)),
      "   Correlation = ",round(corr.iw,3),"\n\n")
  list(pred.pd=pred.pd,pred.tr=pred.tr,pred.ls=pred.ls,pred.iw=pred.iw)
}

# The rose function plots a rose diagram with the lengths of the rays given at
# the end of each ray.  This function takes two arguments: a variogram object
# and the value of the semivariogram (gamma) at which the rose plot is desired.
# =============================================================================
rose <- function(dat.var,critgam){
  numcases <- length(dat.var$dist)
  numdirec <- length(unique(dat.var$dir.hor))
  rose.len <- rep(0,numdirec)
  rose.azi <- rep(0,numdirec)
  j <- 1
  for (i in 1:(numcases-1)){
    prod <- (dat.var$gamma[i] - critgam)*(dat.var$gamma[i+1] - critgam)
    if (dat.var$dir.hor[i] != dat.var$dir.hor[i+1]) prod <- 0
    if (j > 1){
      if (rose.azi[j-1] == dat.var$dir.hor[i]) prod <- 0
    }
    if (prod < 0){
        rose.len[j] <- dat.var$dist[i] + (critgam - dat.var$gamma[i])*
                       (dat.var$dist[i+1] - dat.var$dist[i])/
                       (dat.var$gamma[i+1] - dat.var$gamma[i])
        rose.azi[j] <- as.numeric(unique(dat.var$dir.hor)[j])
        j <- j + 1
    }
  }
  rose.azi <- rose.azi*pi/180
  dscale <- max(rose.len)
  rose.len <- rose.len/dscale
  plot(0,0,type="n",axes=F,xlim=c(-1.1,1.1),ylim=c(-1.1,1.1),
    xlab="E-W Direction",ylab="")
  for (j in 1:(length(rose.len))){
    segments(0,0,(rose.len[j]*sin(rose.azi[j])),(rose.len[j]*cos(rose.azi[j])))
    text((rose.len[j]*sin(rose.azi[j]) + (.1*sin(rose.azi[j]))),
         (rose.len[j]*cos(rose.azi[j]) + (.1*cos(rose.azi[j]))),
         (format(round(rose.len[j]*dscale,2))))
    segments(0,0,-(rose.len[j]*sin(rose.azi[j])),
                 -(rose.len[j]*cos(rose.azi[j])))
  }
  maxx <- max(rose.len*sin(rose.azi))
  text(-(maxx+.25),0,"N-S Direction",srt=90,crt=90)
  title("Rose Diagram",srt=0,crt=0)
  print("Rose Lengths")
  rose.len*dscale
}

# This function takes an nx2 matrix of (x,y) coordinates and rotates the
# coordinates by the angle theta.  The function returns the new set of
# rotated coordinates under the name assigned to the function.  So, if
# you type: "newxy <- rotateaxis(xy,30)", newxy will be a nx2 matrix of the
# rotated (x,y) coordinates (rotated by 30 degrees).

rotateaxis <- function(xy, theta)
{
        pimult <- (theta * 2 * pi)/360
        newx <- c(cos(pimult), sin(pimult))
        newy <- c( - sin(pimult), cos(pimult))
        XY <- as.matrix(xy) %*% cbind(newx, newy)
        as.data.frame(XY)
}

panel.gamma0 <- function(x,y,gamma0=gamma0,span=2/3, ...){
  lofit <- loess.smooth(x,y,span=span)
  panel.xyplot(x,y,...)
  panel.xyplot(lofit$x,lofit$y,type="l")
  nump <- round(length(lofit$x)/2)
  dist0 <- approx(lofit$y[1:nump],lofit$x[1:nump],xout=gamma0)$y
  panel.segments(0,gamma0,dist0,gamma0)
  panel.segments(dist0,0,dist0,gamma0)
  parusr <- par()$usr
  panel.text(max(x),min(y),paste("d0=",format(round(dist0,4)),sep=""),adj=1)
}

getCD <- function(x,y,x0,y0,nugget,sill,range,azim,ratio,model){

# This function finds the matrix of covariances C between all sample locations
# given in the location vectors x,y, and the matrix of covariances D between
# the vector of sample locations (x,y) and the vector of prediction sites
# (x0,y0).  Arguments to this function are:
#     x      = a vector of the x-coordinates of the locations
#     y      = a vector of the y-coordinates of the locations
#     x0     = a vector of the x-coordinates of the k prediction points
#     y0     = a vector of the y-coordinates of the k prediction points
#     nugget = the nugget for the variogram model
#     sill   = the sill for the variogram model (sigma^2)
#     range  = the range for the variogram model
#     azim   = the direction (from 90) of maximum spatial continuity for
#              geometric anisotropy.
#     ratio  = the range ratio associated with the azimuth for geometric
#              anisotropy.
#     model  = the type of variogram model.  Choices are:
#                 "spher" = Spherical
#                 "expon" = Exponential
#                 "gauss" = Gaussian
#
# If the variogram model is isotropic, choose the range ratio to be 0.
#
# The function returns a list with the (n+1)x(n+1) matrix of covariances
# among the sample values labeled C, and the (n+1)xk matrix of covariances
# between the sample values and the prediction points labeled D.  Both
# matrices have the supplemented ones columns supplemented for use in the
# kriging equations.  Each column of values in the D matrix corresponds to
# a different prediction location.

  n <- length(x)
  k <- length(x0)
  xy <- matrix(c(x,x0,y,y0),nrow=(n+k),byrow=F)
  if (ratio > 0){
    az <- (90-azim)*pi/180
    rot <- matrix(c(cos(az),sin(az),-sin(az),cos(az)),nrow=2,byrow=T)
    newxy <- rot%*%t(xy)
    disx <- dist(newxy[1,],method="euclidean")
    disy <- dist(newxy[2,],method="euclidean")
    disxy <- sqrt((disx/range)**2 + (disy*ratio/range)**2)
  }
  else{
    disxy <- dist(xy,method="euclidean")/range
  }
  dmat <- matrix(0,(n+k),(n+k))
  dmat[lower.tri(dmat)] <- disxy
  dmat <- dmat + t(dmat)
  if (model == "spher"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*((1.5*dmat) - (0.5*(dmat)**3))*
            (dmat < 1) + (sill-nugget)*(dmat >= 1)
  }
  if (model == "expon"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*(1 - exp(-3*dmat))
  }
  if (model == "gauss"){
    vmat <- (nugget*(dmat > 0)) + (sill-nugget)*(1 - exp(-3*(dmat**2)))
  }
  
  vmat <- sill - vmat
  C <- rbind(cbind(vmat[1:n,1:n],rep(1,n)),t(rep(1,(n+1))))
  C[n+1,n+1] <- 0
  D <- rbind(matrix(vmat[1:n,(n+1):(n+k)],n,k),t(rep(1,k)))
  list(C=C, D=D)
}

ordkrige <- function(x,y,v,x0,y0,nugget,sill,range,azim,ratio,model){
  covs <- getCD(x,y,x0,y0,nugget,sill,range,azim,ratio,model)
  w <- solve(covs$C)%*%covs$D
  pred <- sum(w[-(length(v)+1)]*v)
  se.pred <- sqrt(sill - sum(w*covs$D))
  list(w=w,pred=pred,se.pred=se.pred)
}
