##' PCAnalysis is a Principal Component Analysis Function 
##' 
##' PCAnalysis is a PCA constructor function
##' 
##' PCAnalysis constructor function uses call by reference. The new object is 
##' put into the \code{analysis} slot of the dataObject on which PCA was 
##' calculated. 
##' @export
##' @param dataObject object of type MassImage
##' @param nComp integer number of components
##' @param ... further args
##' @return PCAnalysis class object
##' @rdname PCAnalysis
##' @examples
##' testImage<-MassImage('dummy')
##' testImage<-PCAnalysis(testImage, 4)
##' image(analysis(testImage, 1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage, nComp = 4)
##' image(analysis(testImage, 1), comp = 1)}
PCAnalysis <- function(dataObject, nComp, ...) {
    if (nrow(as.matrix(dataObject@nz)) < ncol(as.matrix(dataObject@nz))) {
        tNz <- t(dataObject@nz)
        p <- eigen(dataObject@nz %*% tNz)$vectors[, 1:nComp]
        if (nComp > 1) {
            p <- apply(p, 2, function(p) tNz %*% p)
            p <- apply(p, 2, function(p) p/sqrt(sum(p^2)))
        } else {
            p <- tNz %*% p
            p <- p/sqrt(sum(p^2))
        }
        
    } else p <- eigen(t(dataObject@nz) %*% dataObject@nz)$vectors[, 1:nComp]
    
    vec <- dataObject@nz %*% p
    
    newAnalysisObjectPosition <- (length(dataObject@analysis) + 1)
    message(newAnalysisObjectPosition)
    
        dataObject@analysis[[newAnalysisObjectPosition]] <- new(
            "PCAnalysis", 
            pcaLoadings = p, 
            pcaScores = vec, 
            nComp = as.numeric(nComp), 
            imageDim = dim(dataObject), 
            classOfData = as.character(class(dataObject)))
        return(dataObject)
}
