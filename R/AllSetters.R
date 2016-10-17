#########################################################
######
### setters/replacement, ordered by Class, alphabetically
######
#########################################################






######
### MassImage
######

##' @rdname xy
setMethod("xy<-", signature(object="MassImage"),
          function(object, value){
            object@xy<-value
            if(validObject(object))
                return(object)
          })






######
### MassSpectra
######

##' @rdname analysisName
setMethod("analysisName<-", signature(object="MassSpectra"),
          function(object, value){
            object@analysisName<-value
            if(validObject(object))
                return(object)
          })

##' @rdname instrument
setMethod("instrument<-", "MassSpectra",
          function(object,value){
            object@instrument<-value
            if(validObject(object))
                return(object)
          })


##' @rdname nz
setMethod("nz<-", signature(object="MassSpectra"),
          function(object, value){
            object@nz<-value
            if(validObject(object))
                return(object)
          })
##' mz setter method
##' @param value double mass to charge ratio
##' @rdname mz
setMethod("mz<-", signature(object="MassSpectra"),
          function(object, value){
            object@mz<-value
            if(validObject(object))
                return(object)
          })


##' @rdname calibration
setMethod("calibration<-", signature(object="MassSpectra"),
          function(object, value){
            object@calibration<-value
            if(validObject(object))
                return(object)
          })

##' @rdname calibPoints
setMethod("calibPoints<-", signature(object="MassSpectra"),
          function(object, value){
            object@calibPoints<-value
            if(validObject(object))
                return(object)
          })

##' @rdname analysis
setMethod("analysis<-", signature(object="MassSpectra"),
          function(object, value){
            object@analysis<-value
            if(validObject(object))
                return(object)
          })


#######
### PeakList
#######

##' @rdname peakIDs
setMethod("peakIDs<-", signature(object='PeakList'),
          function(object, value){
            object@peakIDs<-value
            if(validObject(object))
                return(object)
          })

##' @rdname peakMzs
setMethod("peakMzs<-", signature(object='PeakList'),
          function(object, value){
            object@peakMzs<-value
            if(validObject(object))
                return(object)
          })





