#' Merge Series and observations
#'
#' @description This function merges gridded fields and gauge data. A number of techniques are available: Double Kernal Smoothing (dk); Kriging with drift (ked); Bayesian combination (bc); multiplicative bias correction (mbc) and inverse distance weighting of the residuals (ridw). Ordinary Kriging (ok) is also provided for the interpolation of gauge data. See details for more informations about the methods.
#'
#' @param grd is a sf object of gridded rainfall data.
#' @param obs is a sf object of observed rainfall data
#' @param method combination method (see details)
#' @param nm names of series to evaluate
#'
#' @return A sf of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection of the names of grd and obs.
#'
#' @details the double kernel smoothing of the residuals following Li & Shao (2010)
#' References: Li, M. and Shao, Q. (2010). 
#' An improved statistical approach to merge satellite rainfall
#' estimates and raingauge data. Journal of Hydrology, 385(1-4):51-64.
#' Silverman, B. (1986). Density Estimation for Statistics and Data Analysis. 
#' Chapman and Hall/CRC Monographs on Statistics and Applied Probability Series. 
#' Chapman and Hall/CRC.

#' @examples
#' # out <- combine(sat, obs,method="dk")
#'
#' @export
combine <- function(grd, obs, method=c("dk","ok","ked","mbc","ridw","bc"),
                    nm = names(grd)){

    ## check method
    method <- match.arg(method)
    
    ## check grd is of correct type
    trgt <- ifelse( method %in% c("ok"), c("sfc","sf"), "sf")
    if(!any(trgt %in% class(grd))){
        stop(paste0("grd should be a ",trgt," object"))
    }

    ## check obs
    if(!("sf" %in% class(obs))){
        stop("obs should be a sf object")
    }

    ## check all objects are points
    if(!all(st_geometry_type(obs)=="POINT")){
        stop("All obs geometries should be points")
    }
    if(!all(st_geometry_type(grd)=="POINT")){
        stop("All grd geometries should be points")
    }
    
    ## check there are some common names
    if( method %in% c("ok") ){
        ## convert grid to an sfc with correct names
        
        grd <- cbind(grd, matrix(numeric(NA),length(grid),ncol(obs)))
        names(grd) <- names(obs)
    }else{
        ## drop geometry from nm
        nm <- setdiff(nm,"geometry")
        ## names must be in grd
        nm <- intersect(nm,names(grd))
        grd <- grd[,c("geometry",nm)]
        ## only need to evaluate those in obs
        nm <- intersect(nm,names(obs))
        if(length(nm)==0){
            warning("grd and obs have no common series") 
            return(grd)
        }
    }
    
    ## ensure the projection information matches
    if( st_crs(obs) == st_crs(grd) ){
        ## they seem identical by may fail in gstat dur to different
        ## froms in the wkt so...Fixed if gstat 2.0.9
        st_crs(obs) <- st_crs(grd)
    }else{
        warning("Reprojecting observation locations")
        obs <- st_transform(obs, st_crs(grd))
    }

    ## call correct function
    grd <- switch(method,
                  bc = bc(grd,obs,nm),
                  dk = dk(grd,obs,nm),
                  ked = ked(grd,obs,nm),
                  mbc = mbc(grd,obs,nm),
                  ok = ok(grd,obs),
                  ridw = ridw(grd,obs,nm),
                  stop("Invalid method")
                  )

    return(grd)
}
