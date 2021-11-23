#' Mean Bias Correction
#'
#' @description  This function calculates a spatial mean multiplicative bias 
#' correction of the satellite field at every time step 
#'
#' @param grd is a SpatialPointsDataFrame of gridded rainfall data.
#' @param obs is a SpatialPointsDataFrame of observed rainfall data
#' @param nm names of columns to evaluate
#'
#' @return A SpatialPointsDataFrame of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs#'
#' @examples
#' # MBC_out <- MBC(sat, gauge)
#' # MBC_cv  <- MBC(sat, gauge, cross.val=TRUE)
#'
#' @export
mbc <- function(grd,obs,nm){
    
    ## compute distances between raster points and gauges
    dro <- st_distance(grd,obs)
    obs_pnt <- apply(dro,2,which.min) ## location of gauges as nearest point

    for (ii in nm){
        idx <- is.finite(obs[[ii]])
        k <- sum( obs[[ii]][idx] ) / sum( grd[[ii]][obs_pnt[idx]] )
        grd[[ii]] <- grd[[ii]] * k
    }
    return(grd)
}
