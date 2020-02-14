# ============================================================================ #
# Function to convert (Long,Lat) to UTM Coordinates: The function requires the #
# coordinates (x = Longitude, y = Latitude) and the longitude zone, which can  #
# be looked up for a given set of coordinates.                                 #
# ============================================================================ #
LongLatToUTM<-function(x,y,zone){
 xy <- data.frame(ID = 1:length(x), X = x, Y = y)
 coordinates(xy) <- c("X", "Y")
 proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
 res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
 return(as.data.frame(res))
}