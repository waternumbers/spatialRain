#' Kriging with External Drift
#'
#'
#' @description This function performs a Kriging interpolation of rain gauge 
#' estimates using gridded estimates as an external drift
#'
#' @param grd is a SpatialPointsDataFrame of gridded rainfall data.
#' @param obs is a SpatialPointsDataFrame of observed rainfall data
#' @param nm names of series to utilise
#'
#' @return A SpatialPointsDataFrame of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs.
#' @examples
#' # KED_out <- KED(sat, gauge)
#' # KED_cv  <- KED(sat, gauge, cross.val=TRUE)
#'
ked <- function(grd,obs,nm){
    
    ## compute distances between raster points and gauges
    dro <- st_distance(grd,obs)
    obs_pnt <- apply(dro,2,which.min) ## location of gauges as nearest point

    for (ii in nm){
        
        o <- obs[,ii]
        o[["obs"]] <- o[[ii]]
        o[["trend"]] <- grd[[ii]][obs_pnt]
        o <- o[is.finite(o[["obs"]]),]

        fc <- grd[,ii]
        fc[["trend"]] <- fc[[ii]]
        
        ## Model semivariogram
        vm <- autofitVariogram(obs ~ trend, as(o,"Spatial"),
                               model = c("Sph", "Exp", "Gau"))

        ## Perform Kriging
        x <- krige(obs ~ trend, locations=obs, newdata=fc,
                   model = vm)[["var1.pred"]]
        x[x<0] <- 0
        
        grd[[ii]] <- x
    }
    
    return(grd)
}
