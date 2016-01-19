##' PrComp is a Principal Component Analysis Function 
##' 
##' PrComp is a PCA constructor function
##' 
##' PrComp constructor function uses call by reference. The new object is put
##' into the \code{analysis} slot of the dataObject on which PCA was calculated. 
##' @export
##' @param dataObject object of class MassSpectra
##' @param ... additional args for prcomp
##' @return object of class PrComp
##' @rdname PrComp
##' @examples
##' testImage<-MassImage('dummy')
##' testImage<-prComp(testImage)
##' image(analysis(testImage, 1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-prComp(testImage)
##' image(analysis(testImage, 1), comp = 1)}
prComp<-function(dataObject, ...){             
  prcompOut<-prcomp(nz(dataObject), ...)
  
  newAnalysisObjectPosition<-(length(dataObject@analysis)+1)
  # print(newAnalysisObjectPosition)
  dataObject@analysis[[newAnalysisObjectPosition]]<-
                           new('PrComp',
                               pcaLoadings = matrix(prcompOut$rotation,
                                                    dim(prcompOut$rotation)[1],
                                                    dim(prcompOut$rotation)[2]),
                               pcaScores   = prcompOut$x,
                               nComp       = dim(prcompOut$x)[2],
                               imageDim    = dim(dataObject),
                               classOfData = as.character(class(dataObject)),
                               sdev        = prcompOut$sdev,
                               center      = prcompOut$center,
                               scale       = prcompOut$scale)
  return(dataObject)
}

##' constructor for PrComp
##' @param object object of class 
##' @return object of class PrComp
setMethod(baseObject,
          "PrComp",
          function(object){
            base <- list(sdev     = object@sdev, 
                         rotation = object@pcaLoadings, 
                         center   = object@center,
                         scale    = object@scale,
                         x        = object@pcaScores)
            attr(base, "class") <- "prcomp"
            return(base)
          }
)


