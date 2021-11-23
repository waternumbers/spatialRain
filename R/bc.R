#' Bayesian Combination 
#'
#'
#' @description This function merges satellite and Ordinary-Kriged rain 
#' gauge fields using Bayesian inference
#' References: Todini, E. (2001). A Bayesian technique for 
#' conditioning radar precipitation estimates to rain-gauge measurements. 
#' Hydrology and Earth System Sciences, 5(2):187-199. 
#' PROGEA Srl (2009). RAINMUSIC, User manual & references. PROGEA Srl, Bologna.
#'
#' @param grd is a SpatialPointsDataFrame of gridded rainfall data.
#' @param obs is a SpatialPointsDataFrame of observed rainfall data
#' @param nm names of series to utilise
#'
#' @return A SpatialPointsDataFrame of merged fields of the same dimension as grd, unless \code{only_obs == TRUE} in which case the columns are determined by the intersection
#'
#' 
bc <- function(grd,obs,nm){

    ## compute spatial distance between grid points
    D <- spDists(grd)
    
    for (ii in nm){
        browser()
        ## ##############################################################
        ## Krige observations (as ok)
        ## Get data for time step and exclude gauges with missing data
        o <- obs[is.finite(obs[[ii]]),ii]
        names(o) <- "obs"
        ## Model semivariogram
        vmo <- autofitVariogram(obs ~ 1, o,
                                model = c("Sph", "Exp", "Gau"))

        ## Perform Kriging 
        kgo <- krige0(obs ~ 1, data=o, newdata=grd,
                      model= vmo$var_model, computeVar=TRUE,
                      fullCovariance = TRUE) ## ouch might be big...

        ## ##############################################################
        ## compute residuals kriging
        ## get Kriging output
        G <- kgo[[1]] # TODO use name
        varG <- kgo[[2]] # TODO use name
        varG[varG < 0] <- 0 ## why? should it ever be <0?

        ## calcualte erros
        S <- grd[,ii]
	errS <- S - G; ## shoule be spatial point df
        names(errS) <- "errS"
        ## errS  <- data.frame(t(errS)); names(errS) <- "errS"
	## coordinates(errS) <- coordinates(Tdata)
	## proj4string(errS) <- proj4string(Tdata)

        ## Model error semivariogram
      	vme <- autofitVariogram(errS ~ 1, errS, model="Exp")

        ## Generate semivariance matrix between all pixels, then convert to covariance
      	Veps        <- variogramLine(object=vme$var_model, dist_vector=D)
      	Veps        <- vme$var_model[2,"psill"] - Veps ##TODO ??
      	Veps[Veps < 0] <- 0# TODO why??

        ## ###################################################
        ## Run Kalman Filter
        
        ## A priori - Pstar
	mu_S  <- mean( S[G>0] - G[G>0] )
	Pstar <- S - mu_S
        ## Kalman gain
	K  <- solve((varG+Veps),Veps) 
        ## Innovation
	Nu <- G - Pstar
        ## A posteriori
	x <- Pstar + K %*% Nu       # prior at grid + correlation*observation
	#x[x<0] <- 0

        grd[,ii] <- x
    }
    grd[grd<0] <- 0

    return(grd)
}


##         ##
##         res <- obs[,ii] - grd[obs_pnt,ii]
##         names(res) <- "res"
##         res <- res[is.finite(res[,"res"]),]

##         ## perform idw
##         r <- idw(res ~ 1,locations=res,newdata=grd,idp=2)$var1.pred
        
##         grd[,ii] <- grd[,ii] + r
##     }
## bc <- function(sat,gauge,longlat=TRUE,cross.val=FALSE){

## #library(rgdal)
## #library(zoo)
## #library(gstat)
## #library(automap)

## if(!cross.val){

## Zs 	<- sat[[1]]
## Zg 	<- gauge[[1]]
## Tdata 	<- sat[[2]]
## Gdata 	<- gauge[[2]]

## # first create gridded rain gauge using OK 

## Zg         <- data.frame(t(Zg))
## gaugename  <- names(Zg)


## ## log transform rain gauge data (comment out if not needed)
## #Zg <- log(Zg +0.01)

## ## merge maps
## data <- cbind(data.frame(Gdata), Zg)
## coordinates(data) <- coordinates(Gdata) #~coords.x1+coords.x2 
## proj4string(data) <- proj4string(Gdata)

## vm.fit    <-  crossval <- maps <- list()

## for (i in 1:length(gaugename)){

## # Get data for time step and exclude gauges with missing data
##   data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

## # Model semivariogram
##   formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,sep=""))
##   vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
##                                   model = c("Sph", "Exp", "Gau"))

## # Perform Kriging 
##   maps[[i]] <-  krige0(formula, data=data_sub, newdata=sat[[2]],
##                        model= vm.fit[[i]]$var_model, computeVar=TRUE,
##                        fullCovariance = TRUE)
	
## }

## # back-log transform output
## #	maps@data  			<- exp(maps@data) - 0.01

## ################################################################################

## vm.fit2 <- maps2 <- list()

## # calculate distances for generating covariance matrix
## distances <- spDists(Tdata,Tdata,longlat)

## for (i in 1:length(gaugename)){

## #get Kriging output
##   G 	     <- maps[[i]][[1]] 
##   varG     <- maps[[i]][[2]]
##   varG [varG < 0] <- 0

##   S     <- Zs[i,]

## # Calculate error
## 	errS <- S - t(G); errS  <- data.frame(t(errS)); names(errS) <- "errS"
## 	coordinates(errS) <- coordinates(Tdata)
## 	proj4string(errS) <- proj4string(Tdata)

## # Model error semivariogram
##       	formula     <- as.formula("errS ~ 1")
##       	vm.fit2[[i]]<- autofitVariogram(formula, errS, model="Exp")

## # Generate semivariance matrix between all pixels, then convert to covariance
##       	Veps        <- variogramLine(object=vm.fit2[[i]]$var_model, dist_vector=distances)
##       	Veps        <- vm.fit2[[i]]$var_model[2,"psill"] - Veps
##       	Veps [Veps < 0] <- 0


## # Run Kalman Filter

## # A priori - Pstar
## 	mu_S  <- mean(t(S)[G>0] - G[G>0])
## 	Pstar <- t(S) - mu_S

## # Kalman gain
## 	K  <- solve((varG+Veps),Veps) 

## # Innovation
## 	Nu <- G - Pstar

## # A posteriori
## 	BC <- Pstar + K %*% Nu       # prior at grid + correlation*observation
## 	BC[BC<0] <- 0
	
## 	Zs[i,] <- BC
## }

##         Zs  <- as.zoo(Zs); time(Zs) <- time(gauge[[1]])

## return(Zs)

## } else {

  
## ################################################################################

## # BC - cross validation

## ################################################################################

## crossval  <- matrix(ncol=nrow( gauge[[2]]), nrow=nrow(gauge[[1]]))

## # initiate progress bar
## pb <- txtProgressBar()
## print("Bayesian combination - cross validation")

## for (p in 1:length(gauge[[2]])){
## setTxtProgressBar(pb, p/nrow(gauge[[2]]))

## Zs 	<- sat[[1]]
## Zg 	<- gauge[[1]]
## Tdata 	<- sat[[2]]
## Gdata 	<- gauge[[2]]

## loc_p   <- which.min(spDists(Tdata,Gdata[p,] ,longlat))

## Zg    	<- Zg[,-p]
## Gdata 	<- Gdata[-p,]

## # first create gridded rain gauge using OK 

## Zg         <- data.frame(t(Zg))
## gaugename  <- names(Zg)

## ## log transform rain gauge data (comment out if not needed)
## #Zg <- log(Zg +0.01)

## ## merge maps
## data <- cbind(Gdata, Zg)
## coordinates(data) <- coordinates(Gdata) #~coords.x1+coords.x2 
## proj4string(data) <- proj4string(Gdata)

## vm.fit     <-  maps <- list()

## for (i in 1:length(gaugename)){

## # Get data for time step and exclude gauges with missing data
## 	data_sub    <- data[is.finite(unlist(as.data.frame(data)[gaugename[i]])),]

## # Model semivariogram
##       	formula     <- as.formula(paste(as.character(gaugename[i])," ~ ", 1 ,
##                                         sep=""))
##       	vm.fit[[i]] <- autofitVariogram(formula, data_sub, 
##                                         model = c("Sph", "Exp", "Gau"))

## # Perform Kriging 
##       	maps[[i]] <-  krige0(formula, data=data_sub, newdata=sat[[2]],
##       	                     model= vm.fit[[i]]$var_model, computeVar=TRUE,
##       	                     fullCovariance = TRUE)
## }


## # back-log transform output
## #	maps@data  			<- exp(maps@data) - 0.01


## ################################################################################

## vm.fit2 <- maps2 <- list()

## # calculate distances for generating covariance matrix
## distances <- spDists(Tdata,Tdata,longlat)

## for (i in 1:length(gaugename)){

## #get Kriging output
##   G 	     <- maps[[i]][[1]] 
##   varG     <- maps[[i]][[2]]
##   varG [varG < 0] <- 0

##   S     <- Zs[i,]

## # Calculate error
## 	errS <- S - G; errS  <- data.frame(t(errS)); names(errS) <- "errS"
## 	coordinates(errS) <- coordinates(Tdata)
## 	proj4string(errS) <- proj4string(Tdata)

## # Model error semivariogram
##        	formula     <- as.formula("errS ~ 1")
##        	vm.fit2[[i]]<- autofitVariogram(formula, errS, model="Exp")

## # Generate semivariance matrix between all pixels, then convert to covariance
##       	Veps        <- variogramLine(object=vm.fit2[[i]]$var_model, 
##                                      dist_vector=distances)
##       	Veps        <- vm.fit2[[i]]$var_model[2,"psill"] - Veps
##       	Veps [Veps < 0] <- 0

## # Run Kalman Filter

## # A priori - Pstar
## 	mu_S  <- mean(t(S)[G>0] - G[G>0])
## 	Pstar <- t(S) - mu_S

## # Kalman gain
##   K  <- solve((varG+Veps),Veps)  

## # Innovation
## 	Nu <- G - Pstar

## # A posteriori
## 	BC <- Pstar + K %*% Nu       # prior at grid + correlation*observation
## 	BC[BC<0] <- 0
	
## 	crossval[i,p] <- BC[loc_p]
## }

## }
## close(pb)

##       	crossval  <- as.zoo(crossval); time(crossval) <- time(gauge[[1]])


## return(crossval)
## }

## }



