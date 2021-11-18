rm(list=ls())
library(terra)


devtools::load_all("./spatialRain")

## read in test data from file - in format for new code
obs <- read.csv.zoo("./spatialRain/inst/extdata/obs.csv")
obs_location <- vect("./spatialRain/inst/extdata/obs_location.geojson")
rst <- rast("./spatialRain/inst/extdata/imerg.nc")[[5]]

## tidy up obs since old code can;t handle missing values.....
obs$Chitreghat <- NULL
obs <- obs[which(rowSums(is.finite(obs))==4), ]
obs <- obs[time(rst),]


## convert to format for old code
sat <- list(zoo(t(values(rst)),terra::time(rst)), ##zoo(matrix(values(rst[[iits]],mat=FALSE),nrow=1),ts),
            crds(rst,na.rm=FALSE))

tmp <- crds(obs_location)
rownames(tmp) <- obs_location$name
tmp <- tmp[names(obs),]
gauge <- list(obs,tmp)

## #####################################################
## double kernal
new <- dk(rst,obs_location,obs)
old <- DS(sat,gauge)
e <- t(values(new)) - old
range(e)

## #####################################################
## Ordinary Kriging
old <- OK(sat,gauge)
