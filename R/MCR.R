##' opaMCR is a MCR-ALS function
##' 
##' opaMCR is a MCR-ALS function using the Orthogonal Projection Approach from 
##' 
##' opaMCR uses the function \code{ChemometricsWithR::opa()} (Orthogonal Projection Approach,
##' CRAN package 'ChemometricsWithR') for start estimates of pure spectras and
##' \code{ALS::als()} (CRAN package 'ALS') as MCR-ALS implementation. This method is doing
##' fine with images up to 256x256 pixels. For larger images, memory usage 
##' becomes unreasonably high. 
##' @export
##' @param dataObject object of class MassImage
##' @param opaComps numeric number of components for the opa method
##' @param maxiter numeric how many iterations
##' @return object of class MCR
##' @rdname MCR
##' @import ChemometricsWithR
##' @import ALS
##' @examples
##' testImage<-MassImage('dummy')
##' testImage<-opaMCR(testImage, 2, 2)
##' image(analysis(testImage,1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-MCR(testImage, 5, 5)
##' image(analysis(testImage,1), comp = 1)
##' }
opaMCR <- function(dataObject, opaComps, maxiter = 10) {
    if (xy(dataObject)[1] * xy(dataObject)[2] * zdim(dataObject) > 
        3.5e+07) 
        stop("Reduce the dataset to below 35 million\n 
             datapoints (x * y * M/z)")
    ### needs the require call as otherwise the dependencies of ALS::als
    ### are not loaded
    requireNamespace('ALS')
    message("Calculating ", opaComps, " initial guess by OPA. ")
    opaOut <- ChemometricsWithR::opa(nz(dataObject), opaComps)
    message("Calculating MCR-ALS ")
    opamcrals <- ALS::als(CList = list(matrix(0, xdim(dataObject) * 
        ydim(dataObject), opaComps)), PsiList = list(nz(dataObject)), 
        S = opaOut$pure.comp, nonnegS = TRUE, nonnegC = TRUE, 
        uniC = FALSE, maxiter = maxiter, forcemaxiter = TRUE, 
        optS1st = FALSE, baseline = TRUE)
    
    
    newAnalysisObjectPosition <- (length(dataObject@analysis) + 
        1)
    message(newAnalysisObjectPosition)
        dataObject@analysis[[newAnalysisObjectPosition]] <- new(
            "MCR", 
            pcaLoadings = opamcrals$S, 
            pcaScores = opamcrals$CList[[1]], 
            nComp = opaComps, 
            imageDim = dim(dataObject), 
            classOfData = as.character(class(dataObject)), 
            RSS = opamcrals$rss, 
            resids = opamcrals$resid[[1]], 
            iters = opamcrals$iter)
        return(dataObject)
    
}