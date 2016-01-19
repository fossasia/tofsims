#######################################################
######
### getters/accessors, ordered by Class, alphabetically
######
#######################################################










######
### MassImage
######

##' Getter, method definition "xdim" on "MassImage"
##' @param object objet of class MassImage
##' @return numeric x dimension of slot xy
setMethod(xdim, signature(object = "MassImage"),
          function(object) object@xy[1])

##' Getter, method definition "ydim" on "MassImage"
##' @param object object of class MassImage
##' @return numeric y dimension of slot xy
setMethod(ydim, signature(object = "MassImage"),
          function(object) object@xy[2])


##' @rdname xy
setMethod(xy, signature(object="MassImage"),
          function(object) object@xy)






######
### MassSpectra
######

##' @rdname analysisName
setMethod("analysisName","MassSpectra",
          function(object){
            object@analysisName
          })

##' @rdname instrument
setMethod("instrument", "MassSpectra",
          function(object){
            object@instrument
          })


##' @rdname nz
setMethod(nz, signature("MassSpectra", "missing"),
          function(object, mzRange){
            object@nz})

##' @rdname nz
setMethod(nz, signature("MassSpectra","numeric"),
          function(object, mzRange){
            if(!length(mzRange)==2      ||
                 !is.numeric(mzRange)   ||
                 !mzRange[1]<mzRange[2])
                stop('check input for mzRange')
            
            low.mz<-min(which(mz(object) >= mzRange[1]))
            high.mz<-max(which(mz(object) <= mzRange[2]))
              
            object@nz[,low.mz:high.mz]
            
          })
##' mz getter method
##' @rdname mz
##' @param object of type MassSpectra
##' @return MassSpectra object with updated mz slot
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' ## access the mz values fo each spectra point
##' mz(testSpectra)[1:100]
##' ## replace a mz value
##' mz(testSpectra)[1] <- 0.000025
##' mz(testSpectra)[1:100]
setMethod(mz, signature("MassSpectra"),
          function(object){
            object@mz
          })

##' @rdname analysis
setMethod(analysis, signature(object="MassSpectra", noAccess="missing"),
          function(object){
            noAnalysis<-length(object@analysis)
            if(noAnalysis==0)
                stop('analysis slot is empty')
            analysisVec<-unlist(lapply(1:noAnalysis, 
                                       function(i) class(object@analysis[[i]])))
            message(data.frame('analysis' = analysisVec))
            
          })


##' @rdname analysis
setMethod(analysis, signature(object="MassSpectra", noAccess="numeric"),
          function(object, noAccess){
            if(length(noAccess) == 1)
                object@analysis[[noAccess]]
            else
                object@analysis[noAccess]
          })

##' @rdname calibration
setMethod(calibration, signature(object="MassSpectra"),
          function(object){
                object@calibration
          })

##' @rdname calibPoints
setMethod(calibPoints, signature(object="MassSpectra"),
          function(object){
                object@calibPoints
          })



######
### PCA
######

##' PCA accessor \code{nComp}, number of component
##' @return nuemric number of components
##' @rdname nComp
setMethod(nComp, signature(object = "PCA"),
          function(object) object@nComp)

##' PCA accessor \code{pcaLoadings}, loading matrix
##' @return matrix numeric with loadings
##' @rdname pcaLoadings
setMethod(pcaLoadings, signature(object = "PCA", comps = "missing"),
          function(object) object@pcaLoadings)

##' PCA accessor \code{pcaLoadings}, loading matrix
##' @return vector or matrix numeric with loadings according comps
##' @rdname pcaLoadings
setMethod(pcaLoadings, signature(object = "PCA", comps = "numeric"),
          function(object, comps) object@pcaLoadings[,as.integer(comps)])

##' PCA accessor \code{pcaScores}, pcaScores matrix
##' @rdname pcaScores
setMethod(pcaScores, signature(object="PCA"),
          function(object) object@pcaScores)

##' PCA accessor \code{pcaScores}, pcaScores matrix
##' @return vector or matrix numeric with scores according comps
##' @rdname pcaScores
setMethod(pcaScores, signature(object = "PCA", comps = "numeric"),
          function(object, comps) object@pcaScores[,as.integer(comps)])

##' method xdim() for PCA class object
##' @param object object of class PCA
##' @return numeric x dimension of image
setMethod(xdim,signature(object = "PCA"),
          function(object){
            if(object@classOfData!='MassImage')
                stop('PCA analysis not based on image data')
            object@imageDim[1]           
          })

##' method ydim() for PCA class object
##' @param object object of class PCA
##' @return numeric y dimension of image
setMethod(ydim,signature(object = "PCA"),
          function(object){
            if(object@classOfData!='MassImage')
                stop('PCA analysis not based on image data')
            object@imageDim[2]           
          })






######
### PeakList
######

##' nPeaks accessor/getter \code{nPeaks} for PeakList Class
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' testSpectra<-calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra<-calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra<-recalibrate(testSpectra)
##' testSpectra<-unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' nPeaks(testSpectra)
##' @rdname nPeaks
setMethod(nPeaks, signature(object="PeakList"),
          function(object){
            if(is.null(peakIDs(object))){
              nPeaks<-0
            } else nPeaks<-dim(peakIDs(object))[2]
                return(nPeaks)
          })

##' @rdname peakIDs
setMethod(peakIDs, signature(object='PeakList'),
          function(object){
                object@peakIDs             
          })

##' @rdname peakMzs
setMethod(peakMzs, signature(object='PeakList'),
          function(object){
                object@peakMzs             
          })


##' MCR accessor iters,
##' @param object object of class MCR
##' @return iters from object
setMethod(iters, signature(object="MCR"),
          function(object) object@iters)

##' MCR accessor RSS,
##' @param object object of type MCR
##' @return RSS from object
setMethod(RSS, signature(object="MCR"),
          function(object) object@RSS)

##' MCR accessor resids,
##' @param object object of class MCR
##' @return resids from object
setMethod(resids, signature(object="MCR"),
          function(object) object@resids)



