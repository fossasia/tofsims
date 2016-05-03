##' ToF-SIMS Toolbox
##'
##' \tabular{ll}{
##' Package: \tab tofsims\cr
##' Type: \tab Package\cr
##' Version: \tab 0.99.2\cr
##' Date: \tab 15-01-2016\cr
##' License: \tab GPL-3\cr
##' LazyLoad: \tab yes \cr
##' }
##'
##' Toolbox for Time-of-Flight Secondary Ion Mass-Spectrometry
##' (ToF-SIMS) data processing and analysis. The package facilitates 
##' importing of raw data files, loading preprocessed data and a range of 
##' multivariate analysis methods that are most commonly applied in the 
##' ToF-SIMS community.
##'
##' @name tofsims-package
##' @aliases tofsims-package
##' @docType package
##' @title ToF-SIMS Toolbox (tofsims)
##' @author Lorenz Gerber \email{lorenz.gerber@slu.se}
##' @author Viet Mai Hoang \email{hviet.0906@gmail.com}
##' @keywords package
##' @import methods
##' @importFrom graphics axis polygon locator
##' @importFrom stats cov prcomp runif sd var princomp
##' @importFrom utils read.csv unzip head
##' @importFrom ProtGenerics mz 'mz<-'
##' @importFrom grDevices blues9 colorRampPalette rainbow
##' @importFrom graphics abline box legend par
##' @importFrom stats embed lm predict smooth.spline
##' @importFrom utils tail
##' 
##' @importFrom Rcpp evalCpp
##' @useDynLib tofsims, .registration=TRUE
NULL