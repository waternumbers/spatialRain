#' Double Kernel Residual Smoothing
#'
#' @description This function merges gridded and rain gauge fields using 
#' the double kernel smoothing of the residuals following Li & Shao (2010)
#' References: Li, M. and Shao, Q. (2010). 
#' An improved statistical approach to merge satellite rainfall
#' estimates and raingauge data. Journal of Hydrology, 385(1-4):51-64.
#' Silverman, B. (1986). Density Estimation for Statistics and Data Analysis. 
#' Chapman and Hall/CRC Monographs on Statistics and Applied Probability Series. 
#' Chapman and Hall/CRC.
#'
#' @param rst is a SpatRaster of gridded rainfall data.
#' @param obs_location is a SpatVector of observation locations.
#' @param obs is a zoo object of observations
#' @param name_attr the column of attributes in obs_location used to name the columns in obs
#'
#' @return A SpatRaster of merged fields of the same dimension as rst
#'
#' @examples
#' # out <- dk(sat, gauge_loc, gauge)
#'
#' @export

dk <- function(rst,obs_location,obs,name_attr="name"){
    
    ## check rst
    if(!("SpatRaster" %in% class(rst))){ stop("rst should be a SpatRaster object") }
    if( length(time(rst)) != nlyr(rst) ){ stop("rst should have time attribute for each layer") }

    ## check obs_loc
    if(!("SpatVector" %in% class(obs_location))){ stop("obs_location should be a SpatialVector") }
    if(!is.points(obs_location)){ stop("observed locations should be points") }
    if(!(name_attr %in% names(obs_location))){ stop("name_attr is missing") }
    
    ## check obs
    if(!("zoo" %in% class(obs))){ stop("obs should be a zoo object") }

    ## ensure the project information matches
    if( crs(obs_location) != crs(rst) ){
        warning("Reprojecting observation locations")
        obs_location <- project(obs_location,crs(rst))
    }

    ## ######################################################
    ## TODO - currently this used spDists - it should use terra::distance
    ## There is a bug in terra 1.4.19 which means some distance
    ## results are not ordered correctly
    
    ## convert observation locations to a matrix for spDists
    obs_loc <- crds(obs_location)
    rownames(obs_loc) <- unlist( obs_location[[name_attr]] )
    
    ## ensure just use observation with data and locations
    obs_name <- intersect(names(obs),rownames(obs_loc))
    obs_loc <- obs_loc[obs_name,]
    
    ## convert raster and get information
    rst_points <- crds(rst,na.rm=FALSE) ## cells as points leave in NA so order matches values
    rst_times <- time(rst) ## times
    lonlat <- is.lonlat(rst) ## is the projection in lonitude and latitude

    ## work out which layers to process
    lyrs <- (1:nlyr(rst))[ rst_times %in% index(obs) ]
    ## if nothing to do return
    if(length(lyrs)==0){
        warning("No matching time stamps")
        return(rst)
    }
    
    ## compute distances between raster points and gauges
    dro <- sp::spDists(rst_points,obs_loc,longlat=lonlat)
    obs_cell <- setNames( apply(dro,2,which.min), obs_name) ## location of gauges as cell numbers

    ## distance between raster points - can be very big!
    drr <- sp::spDists(rst_points,longlat=lonlat)
    bgrid <- ((4/(3*length(drr)))^(1/5))*sd(drr)    ## kernal bandwidth for the grid
    drr <- exp(-0.5*(drr/bgrid)^2)     ## change distance to weight
    
    ## initialise the output
    out <- rst

    ## initialise some tempory variables
    eSS <- eDS <- w <- rep(0,ncell(rst))

    ## compute repeatedly used values
    rsdrr <- rowSums(drr) ## row sums use in

    ## loop layers
    for(lyr in lyrs){

        ## extract values from raster and observations
        x <- values(rst[[lyr]],mat=FALSE) ## vector of values from raster
        y <- as.numeric( obs[rst_times[lyr],obs_name] ) ## vector of observations
                            
        ## work out error
        eta <- setNames( x[obs_cell] - y, obs_name ) ## named residuals

        ## skip is all observations missing
        idx <- is.finite(eta)
        if(!any(idx)){ next}
            
        ## estimate station weights depending on combination
        dd <- dro[,idx,drop=FALSE] ## keep just info for gauges that have obs
        bobs <- ((4/(3*length(dd)))^(1/5))*sd(dd) ## kernal bandwidth
        dd <- exp(-0.5*(dd/bobs)^2) ## distance
            
        ## compute eSS
        eSS <- dd %*% eta[idx]
        w <- rowSums(dd)
        eDS <- eSS
        eSS <- eSS / w
            
        ## compute eDS
        eDS <- eDS + drr %*% eSS
        w <- w + rsdrr
        eDS <- eDS / w
            
        x <- x - eDS
        x[x<0] <- 0
        values(out)[,lyr] <- x
    }
    
    return(out)
}
