##' PrinComp is a Principal Component Analysis Function 
##' 
##' PrinComp is a PCA constructor function
##' 
##' PrinComp constructor function uses call by reference. The new object is put
##' into the \code{analysis} slot of the dataObject on which PCA was calculated. 
##' @export
##' @param dataObject object of class MassSpectra
##' @param ... additional args
##' @return object of class prinComp
##' @rdname PrinComp
##' @examples
##' testImage <- MassImage('dummy')
##' testImage<-prinComp(testImage)
##' image(analysis(testImage, 1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-prinComp(testImage)
##' image(analysis(testImage), 1), comp = 1)}
prinComp<-function(dataObject, ...){             
            princompOut<-princomp(nz(dataObject), ...)
            
            newAnalysisObjectPosition<-(length(dataObject@analysis)+1)
            # print(newAnalysisObjectPosition)
                dataObject@analysis[[newAnalysisObjectPosition]]<-new(
                    'PrinComp',
                    pcaLoadings = matrix(princompOut$loadings,
                                         dim(princompOut$loadings)[1],
                                         dim(princompOut$loadings)[2]),
                    pcaScores   = princompOut$scores,
                    nComp       = dim(princompOut$scores)[2],
                    imageDim    = dim(dataObject),
                    classOfData = as.character(class(dataObject)),
                    sdev        = princompOut$sdev,
                    center      = princompOut$center,
                    n.obs       = princompOut$n.obs,
                    scale       = princompOut$scale,
                    call        = princompOut$call)
            return(dataObject)
}

##' constructor for PrinComp
##' @param object object with class
##' @return object of class PrinComp
setMethod(baseObject,
          "PrinComp",
          function(object){
            base <- list(sdev     = object@sdev, 
                         loadings = object@pcaLoadings, 
                         center   = object@center,
                         scale    = object@scale,
                         n.obs    = object@n.obs,
                         scores   = object@pcaScores,
                         call     = object@call)
            attr(base, "class") <- "princomp"
            return(base)
          }
)



