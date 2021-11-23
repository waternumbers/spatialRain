#' Ordinary Kriging
#'
#' @description Ordinary Kriging
#'
#' @param grd is a SpatialPoints object of predicted locations
#' @param obs is a SpatialPointsDataFrame of observed rainfall data
#'
#' @return A SpatialPointsDataFrame of Ordinary Kriged feilds. The number of points is setermined by grd with the number of series from obs. The Spatial Correlation is determined for each series seperately
#'
#' @examples
#' # out <- ok(sat, obs)
#'
#' @export
ok <- function(grd,obs){

    ## loop layers
    for(ii in names(obs)){
        
        ## Get data for time step and exclude gauges with missing data
	o <- obs[is.finite(obs[[ii]]),ii]
        o[["obs"]] <- o[[ii]]

        ## Model semivariogram
        ## blow up - try specifying the crs - see gstat bug reports
        vm <- autofitVariogram(obs ~ 1, o, model = c("Sph", "Exp", "Gau"))
        
        x <- krige(obs ~ 1, locations=obs, newdata=grd, 
                   model = vm$var_model)$var1.pred
        x[x<0] <- 0

        grd[[ii]] <- x
    }
    return(grd)
}
