## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './spatialRain'
devtools::document(pacPath)
devtools::check(pacPath,remote=TRUE)

## check documentation build
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)
pkgdown::clean_site(pacPath)

## build, populate drat
## linux
dratPath <- "~/Documents/Software/drat"
tmp <- devtools::build(pacPath)
install.packages(tmp)
drat::insertPackage(tmp,dratPath)

## mac and windows
rhub::validate_email() # for first time that session
pkgName <- sub('\\.tar.gz$', '', basename(tmp)) 
## rhub::platforms()[,1] # lists platforms
mch <- rhub::check(path = tmp,
                   platform = c("macos-highsierra-release-cran","windows-x86_64-release",
                                "macos-m1-bigsur-release"))

ext <- c(".tgz",".zip",".tgz")
for(ii in 1:2){ ## m1 not fixed yet in drat
    tmp <- paste0(pkgName,ext[ii])
    download.file(file.path(mch$urls()$artifacts[ii],tmp),tmp)
    drat::insertPackage(tmp,dratPath)
}

## tidy up drat
drat::pruneRepo(dratPath,pkg="FKF",remove="git")


## worth runnign prior to r submission is not done locally
mch <- rhub::check_with_valgrind(path = tmp)
