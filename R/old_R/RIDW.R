#' Residual IDW interpolation
#'
#' @author Bastian Manz
#'
#' @description  This function calculates the residual (additive bias) at every 
#' pixel-point pair and interpolates them to the satellite field using inverse-
#' distance weighting (IDW) at every time step 
#'
#' @param grd is a SpatialPointsDataFrame of gridded rainfall data.
#' @param obs is a SpatialPointsDataFrame of observed rainfall data
#' @param only_obs return only columns for which there are observations
#'
#' @return A SpatialPointsDataFrame of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs.
#' @examples
#' # RIDW_out <- RIDW(sat, gauge)
#' # RIDW_cv  <- RIDW(sat, gauge, cross.val=TRUE)
#'
#' @export
ridw <- function(grd,obs,only_obs=FALSE){

        ## check grd
    if(!("SpatPointsDataFrame" %in% class(grd))){ stop("grd should be a SpatialPointsDataFrame object") }
    
    ## check obs
    if(!("SpatPointsDataFrame" %in% class(obs))){ stop("obs should be a SpatialPointsDataFrame") }
    
    ## check there are some common names
    nm <- intersect(names(grd),names(obs))
    if(length(nm)==0){
        if(only_obs){
            stop("No common series")
        }else{
            warning("No common series")
            return(grd)
        }
    }

       ## ensure the project information matches
    if( crs(obs) != crs(grd) ){
        warning("Reprojecting observation locations")
        obs <- project(obs,crs(grd))
    }

    ## initialise the output - we will rewrite grd
    if(only_obs){
        grd <- grd[,nm]
    }
    

Zs 	<- sat[[1]]
Zg 	<- gauge[[1]]
Tdata 	<- sat[[2]]
Gdata 	<- gauge[[2]]

# normalize vgm

#x-min(x)/min(x)-max(x)

# ensure Gaussian distribution (log-normal transformation, Box-Cox,normal score transform)

# get location of Zg in Zs
loc <- numeric()
for (i in 1:length(Gdata)) loc[i] <- which.min(spDists(Tdata,Gdata[i,] 
                                                       ,longlat))

Zg       <- as.data.frame(t(Zg))
Zs_field <- as.data.frame(t(Zs))
Zs_trend <- as.data.frame(t(Zs[,loc]))


names(Zg)       <- gsub("*-","_",paste("Gauge_",names(Zs_trend),sep=""))
names(Zs_trend) <- gsub("*-","_",paste("Trend_",names(Zs_trend),sep=""))
names(Zs_field) <- names(Zs_trend)

gaugename<- names(Zg)
trendname<- names(Zs_trend)

## merge maps
data_df  <- cbind(Gdata, Zg, Zs_trend) 
data <- data_df

coordinates(data)     <- coordinates(Tdata[loc,])  
proj4string(data)     <- proj4string(Tdata)

coordinates(Zs_field) <- coordinates(Tdata) 
proj4string(Zs_field) <- proj4string(Tdata)

if (cross.val==FALSE){
  
  #resids_g  <- matrix(ncol=ncol(Zs_trend),nrow=nrow(Zs_trend))
  resids    <- matrix(ncol=ncol(Zs),nrow=nrow(Zs))
  Zs        <- matrix(ncol=nrow(Zs_field),nrow=length(gaugename))

  for (i in 1:length(gaugename)){
  
  print(i)
  
  # Get data for time step and exclude gauges with missing data
  data_sub1   <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]
  
  #resids_g[,i]<- Zg[,i] - Zs_trend[,i]
  resids_g <- Zg[,i] - Zs_trend[,i]
  resids_g <- resids_g[which(is.finite(resids_g))] 
  
  data_sub1@data[,gaugename[i]]  <- resids_g#[,i]
  formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
  resids[i,] <- idw(formula,locations=data_sub1,newdata=Zs_field,idp=2)$var1.pred 
  Zs[i,] <- as.numeric(as.matrix(Zs_field[,i]@data)) + resids[i,]
  
  Zs[i,Zs[i,]<0] =0

  }


  # tag dates 
       Zs  <- as.zoo(Zs); time(Zs) <- time(gauge[[1]])

  return(Zs) 

} else if(cross.val==TRUE){
  
  resids_g    <- matrix(ncol=ncol(Zs_trend),nrow=(nrow(Zs_trend)-1))
  resids       <- matrix(ncol=ncol(Zs),nrow=nrow(Zs))
  crossval        <-  matrix(ncol=nrow(Zs_trend),nrow=ncol(Zs_trend)) 
  
  for(i in 1:length(gaugename)){
    
    # Get data for time step and exclude gauges with missing data
    data_sub2    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]
    data_sub2 <- data_sub2[,c(gaugename[i],trendname[i])]
    
    for(j in 1:nrow(Zg)){
      
      print(paste(i,".",j,sep=""))
      
      data_sub2x <- data_sub2[-j,]
      
      resids_g[,i]<- data_sub2x@data[,1] - data_sub2x@data[,2]
      data_sub2x@data[,1]  <- resids_g[,i]
      formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
      resids[i,] <- idw(formula,locations=data_sub2x,newdata=Zs_field,idp=2)$var1.pred 
      crossval[i,j] <- as.numeric(as.matrix(Zs_field[,i]@data))[loc[j]] + resids[i,loc[j]]
 
    }
    
    crossval[i,crossval[i,]<0] =0
  }
  
  # tag dates 
  crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])
  
  return(crossval)
  
}
}

