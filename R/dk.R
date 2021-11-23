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
#' @param grd is a sf of gridded rainfall data.
#' @param obs is a sf of observed rainfall data
#' @param nm names of series to utilise
#'
#' @return A sf of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs.
#'
#' @examples
#' # out <- dk(sat, obs)
#'
dk <- function(grd,obs,nm){
    
    ## compute distances between raster points and gauges
    dro <- st_distance(grd,obs)
    units(dro) <- NULL
    ## compute distances between grid points - can be very big
    drr <- st_distance(grd)
    units(drr) <- NULL

    ## location of gauges as nearest point
    obs_pnt <- apply(dro,2,which.min)

    ## kernal for grid point weighting
    bgrid <- ((4/(3*length(drr)))^(1/5))*sd(drr)    ## kernal bandwidth for the grid
    drr <- exp(-0.5*(drr/bgrid)^2)     ## change distance to weight
    
    ## initialise some tempory variables
    eSS <- eDS <- w <- rep(0,nrow(grd))

    ## compute repeatedly used values
    rsdrr <- rowSums(drr) ## row sums use in

    ## loop layers
    for(ii in nm){

        ## extract values from raster and observations
        x <- grd[[ii]] ## vector of values from raster
        y <- obs[[ii]] ## vector of observations
                            
        ## work out error
        eta <- x[obs_pnt] - y

        ## skip is all observations missing
        idx <- is.finite(eta)
        if(!any(idx)){ next }
            
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
        grd[[ii]] <- x
    }
    
    return(grd)
}
