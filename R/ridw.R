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
#' @param nm names of series to utilise
#'
#' @return A SpatialPointsDataFrame of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs.
#' @examples
#' # RIDW_out <- RIDW(sat, gauge)
#' # RIDW_cv  <- RIDW(sat, gauge, cross.val=TRUE)
#'
ridw <- function(grd,obs,nm){
    
    ## compute distances between raster points and gauges
    dro <- st_distance(grd,obs)
    obs_pnt <- apply(dro,2,which.min) ## location of gauges as nearest point
    
    for (ii in nm){
        res <- obs[,ii]
        res[[ii]] <- res[[ii]] - grd[[ii]][obs_pnt] ## get rid of units
        res <- res[is.finite(res[[ii]]),]
        if(length(res)==0){ next }

        ## perform idw
        r <- idw(formula(paste0(ii," ~ 1")),locations=res,newdata=grd,idp=2)[["var1.pred"]]
        x <-  grd[[ii]] + r
        x[x<0] <- 0
                
        grd[[ii]] <- x
    }
    
    return(grd)
}
