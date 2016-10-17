##' @title PeakList Constructor
##' @description \code{PeakList} class constructor
##' @details The \code{PeakList} class constructor is used to construct a new
##' PeakList object. Input are currently all needed variables.
##' @export
##' @param analysisName character vector with the import filename
##' @param instrument character vector type of instrument used in the experiment
##' @param calibration data frame for numerics slope and intercept of the mass calibration
##' @param calibPoints data frame for time of flight to maass to charge calibration
##' @param nz matrix numeric containing ion counts, rows are image points, column toftimes/mass to charge ratios
##' @param mz vector same length as columns in \code{nz} for mass to charge values
##' @param peakIDs matrix integer ID for peaks
##' @param peakMzs matrix with mass to charge values for lower, middle and upper peak values
##' @param ... additional args
##' @return object of class PeakList
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
PeakList <- function(analysisName = NULL, 
                     instrument   = NULL,
                     nz           = NULL,
                     calibration  = NULL,
                     calibPoints  = NULL,
                     mz           = NULL,
                     peakIDs      = NULL,
                     peakMzs      = NULL,...){
  
  if(!exists('analysisName', mode = 'character')  ||
       !exists('instrument', mode = 'character')  ||
       !exists('nz', mode = 'numeric')            ||
       !exists('mz', mode = 'numeric')            ||
       !(exists('peakIDs', mode = 'numeric')      ||
         exists('peakIDs', mode = 'NULL'))        ||
       !(exists('peakMzs', mode = 'numeric')      ||
         exists('peakMzs', mode = 'NULL')))
    stop('missing parameters')
  
  
  new("PeakList",
      analysisName  = as.character(analysisName),
      instrument    = as.character(instrument),
      nz            = as.matrix(nz),
      calibration   = as.vector(calibration),
      calibPoints   = as.vector(calibPoints),
      mz            = as.vector(mz),
      peakIDs       = peakIDs,
      peakMzs       = peakMzs)
}

##' method defining \code{show()} for the \code{MassSpectra} class
##' show has a generic by default
##' @param object object of class PeakList
##' @return data.frame character
setMethod("show", signature(object = "PeakList"),
          function(object){
            print(data.frame(analysisName               = analysisName(object),
                             instrument                 = instrument(object),
                             'spectra.points'           = zdim(object),
                             nPeaks                     = nPeaks(object))
            )
          }
)

##' Method \code{plot()} for \code{MassSpectra}
##' 
##' Method defining \code{plot()} for the \code{MassSpectra} class
##' plot has no generic by default
##' 
##' The output of this method is adapted for plotting mass spectra. Uncalibrated
##' data is plotted as xy plot while uncalibrated data is plotted as barplot. The parameter
##' \code{mzRange} allows choosing the plot range directly according to the mz number
##' (when calibrated).
##' @param x object of type PeakList
##' @param y missing
##' @param ... further args
##' @param mzRange vector or length two, indicating the mz range to be plotted
##' @param plotDeriv boolean plot derivate if available
##' @param plotPeaks boolean plot peaks if available
##' @param plotWidths boolean plot peak widths if available
##' @return plot spectra with peaks and peak widths
setMethod("plot", signature(x = "PeakList",y="missing"),
          function(x,y,..., mzRange=c(0,200), plotDeriv=FALSE, plotPeaks=TRUE, plotWidths=TRUE){
            requireNamespace('RColorBrewer')
            low.mz<-1
            high.mz<-length(mz(x))
            
            if(!is.null(mzRange)){
              low.mz<-min(which(mz(x) >= mzRange[1]))
              high.mz<-max(which(mz(x) <= mzRange[2]))
            }
            
            
            
            toploty<-as.vector(nz(x))[low.mz:high.mz]
            toplotx<-mz(x)[low.mz:high.mz]
            
            plot(toplotx,toploty,ylab='ioncounts',xlab='M/z',...)
            
            
            if(!is.null(peakIDs(x))){
              
              mzvaluelow<-mz(x)[low.mz]
              peaks<-peakMzs(x)[2,]
              
              if(plotPeaks){
                abline(v=peaks, lwd=2)
              }
              
              if(plotWidths){
                if(max(peakMzs(x)[1,])>0){
                  lower<-peakMzs(x)[1,]
                  abline(v=lower,col=RColorBrewer::brewer.pal('Set1',n=9))
                }
                
                if(max(peakMzs(x)[3,])>0){
                  upper<-peakMzs(x)[3,]
                  abline(v=upper,col=RColorBrewer::brewer.pal('Set1',n=9))
                }
              }
              
              if(dim(peakMzs(x))[2]>1){
                if(plotDeriv){
                  toploty<-as.vector(peakIDs(x)[4,low.mz:high.mz])
                  toplotx<-peakMzs(x)[4,low.mz:high.mz]
                  par(new=TRUE)
                  plot(toplotx, toploty,,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
                  #points(toplotx,toploty,type='l',col='red')
                }
              }
            }
          }
)


##' @rdname peaks2Spectra
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-reduceSpectrumResolution(testSpectra,everyN = 4, mode = 'keep')
##' peakPickSpectra<-testSpectra
##' peakPickSpectra<-calibPointNew(peakPickSpectra, mz = 15, value = 15.01551)
##' peakPickSpectra<-calibPointNew(peakPickSpectra, mz = 181, value = 181.0228)
##' peakPickSpectra<-recalibrate(peakPickSpectra)
##' peakPickSpectra<-unitMassPeaks(peakPickSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' par(mfcol = c(1,2))
##' plot(testSpectra, mzRange = c(38.5, 40.5), type = 'l')
##' testSpectra<-peaks2Spectra(peakPickSpectra, testSpectra)
##' plot(testSpectra, mzRange = c(38.5, 40.5), type = 'l')
setMethod(peaks2Spectra, signature(objectPeaks="PeakList", objectSpectra="MassSpectra"),
          function(objectPeaks,objectSpectra){
            message('Matching lower peak Mz\'s.')
            lowerIDs<-unlist(lapply(peakMzs(objectPeaks)[1,], function(i) findClosestMatch(i, mz(objectSpectra), 'lower')))
            message('Matching center peak Mz\'s.')
            centerIDs<-unlist(lapply(peakMzs(objectPeaks)[2,], function(i) findClosestMatch(i, mz(objectSpectra), 'upper')))
            message('Matching upper peak Mz\'s.')
            upperIDs<-unlist(lapply(peakMzs(objectPeaks)[3,], function(i) findClosestMatch(i, mz(objectSpectra), 'upper')))
            
            
            
            PeakList(analysisName = analysisName(objectSpectra),
                                                           instrument   = instrument(objectSpectra),
                                                           nz           = nz(objectSpectra),
                                                           calibration  = calibration(objectSpectra),
                                                           calibPoints  = calibPoints(objectSpectra),
                                                           mz           = mz(objectSpectra),
                                                           peakIDs      = rbind(lowerIDs,centerIDs,upperIDs),
                                                           peakMzs      = peakMzs(objectPeaks))
            #return(tempPeakList)
                        
          })

##' method findPeakWidth
##' 
##' method findPeakWidth
##' 
##' This method uses signal processing to determine lower and upper peak width
##' limits based on local max/min detection of the first derivate next to
##' peak center values. The initial code for local min/max detection is adapted
##' from the CRAN package 'ChemometricsWithR'.
##' @rdname findPeakWidth
##' @import signal
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testPeakList<-PeakList(analysisName = analysisName(testSpectra),
##' instrument = instrument(testSpectra),
##' nz = nz(testSpectra),
##' calibration = calibration(testSpectra),
##' calibPoints = calibPoints(testSpectra),
##' mz = mz(testSpectra),
##' peakIDs = NULL,
##' peakMzs = NULL)
##' par(mfcol=c(1,2))
##' plot(testPeakList, mzRange=c(25,32), type = 'l')
##' testPeakList<-addPeaks(testPeakList, mzs=26:31, width=0.4)
##' testPeakList<-findPeakWidth(testPeakList, p = 3, n = 199, 
##' span = 100, widthExtLower = 2, widthExtUpper = 2)
##' plot(testPeakList, mzRange=c(25,32), type = 'l')
setMethod(findPeakWidth, signature(object="PeakList"),
          function(object, p=3, n=199, span=100, widthExtLower=1.7, widthExtUpper=2, ...){
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            message('')
            message('taking first derivate of spectra')
            firstDiff<-diff(nz(object)[1,])
            
            ###filtering the first diff with Savitzky-Golay
            message('applying Savitzky-Golay filter')
            requireNamespace('signal')
            sgolay.filter<-signal::sgolay(p=p,n=n,m=0,ts=1)
            filtered<-signal::filter(filt = sgolay.filter, x = as.vector(firstDiff))
            
            ###make vector with local minima
            message('determine local maxima of first derivate')
            spanWidth <- span * 2 + 1
            localMax <- spanWidth + 1 - apply(embed(filtered, spanWidth), 1, which.max)
            localMax[localMax == 1 | localMax == spanWidth] <- NA
            max.pks <- localMax + 0:(length(localMax) - 1)
            maxPeaks<-unique(max.pks[!is.na(max.pks)])
            
            ###make vector with local maxima
            message('determine local minima of first derivate')
            spanWidth <- span * 2 + 1
            localMin <- spanWidth + 1 - apply(embed(filtered, spanWidth), 1, which.min)
            localMin[localMin == 1 | localMin == spanWidth] <- NA
            min.pks <- localMin + 0:(length(localMin) - 1)
            minPeaks<-unique(min.pks[!is.na(min.pks)])
            
            message('determine minima and maxima next to peak center')
            tempUpperWidth<-NULL
            tempLowerWidth<-NULL
            peakCounter<-1
            
            for(i in peakIDs(object)[2,]){
              tempUpperWidth<-c(tempUpperWidth,minPeaks[head(which(minPeaks>i),1)])
              tempLowerWidth<-c(tempLowerWidth,maxPeaks[tail(which(maxPeaks<i),1)])
              
              if(length(tempUpperWidth)<peakCounter){tempUpperWidth<-c(tempUpperWidth,NaN)}
              if(length(tempLowerWidth)<peakCounter){tempLowerWidth<-c(tempLowerWidth,NaN)}
              peakCounter<-peakCounter+1
            }
            
            message('Removing peaks with missing width values.')
            lackingWidth<-unique(c(which(is.na(tempUpperWidth)),which(is.na(tempLowerWidth))))
            if(length(lackingWidth)>0){
              tempUpperWidth<-tempUpperWidth[-lackingWidth]
              tempLowerWidth<-tempLowerWidth[-lackingWidth]
              peakIDs(object)<-peakIDs(object)[,-lackingWidth]
              peakMzs(object)<-peakMzs(object)[,-lackingWidth]
            }
            
            message('applying the peak width extension factors')
            ### Width extension is a factor multiplied wiht the M/z values
            ### Then the closest available M/z values in the dataset have to be found (slow)
            ### From looking up the closest M/z values, ID values result
            extendedUpper<-mz(object)[tempUpperWidth]+(mz(object)[tempUpperWidth]-peakMzs(object)[2,])*widthExtUpper
            upperIDs<-unlist(lapply(extendedUpper, function(i) findClosestMatch(i,mz(object),'upper')))
            extendedLower<-mz(object)[tempLowerWidth]-(peakMzs(object)[2,]-mz(object)[tempLowerWidth])*widthExtLower
            lowerIDs<-unlist(lapply(extendedLower, function(i) findClosestMatch(i,mz(object),'lower')))
            
            
            peakIDs(object)[1,]<-lowerIDs
            peakIDs(object)[3,]<-upperIDs
            peakMzs(object)[1,]<-mz(object)[lowerIDs]
            peakMzs(object)[3,]<-mz(object)[upperIDs]
            return(object)
                        
            
          } )     

##' removePeaks for PeakList Class allows removing peaks below a certain treshhold of ioncounts.
##' the threshhold is not calculated as area, but just from the peak height (ion count at peak
##' center)
##' @rdname removePeaks
setMethod(removePeaks, signature(object="PeakList", mzs = 'missing', 
                                 operator = 'missing', limit='numeric', 
                                 nLocator="missing"),
          function(object, mzs, operator, limit, nLocator, ...){
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            toRemove<-which(nz(object)[match(peakMzs(object)[2,],mz(object))]<=limit)
            peakIDs(object)<-peakIDs(object)[,-toRemove]
            peakMzs(object)<-peakMzs(object)[,-toRemove]

            message('',length(toRemove),' peaks removed')
                        
            return(object)
          })

##' removePeaks for PeakList Class allows removing peaks manually
##' @rdname removePeaks
setMethod(removePeaks, signature(object="PeakList", mzs = 'missing',
                                 operator = 'missing', limit='missing', 
                                 nLocator="numeric"),
          function(object, mzs, operator, limit, nLocator, ...){
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            coords<-locator(n=nLocator)
            peaksMZ<-peakMzs(object)[2,]
            
            peakNr<-sapply(1:length(coords$x), function(i) which(abs(peaksMZ-coords$x[[i]])==min(abs(peaksMZ-coords$x[[i]]))))
            message(peakNr)
            peakIDs(object)<-peakIDs(object)[,-peakNr]
            peakMzs(object)<-peakMzs(object)[,-peakNr]
                       
            return(object)
            
          })

##' removePeaks for PeakList Class allows removing peaks manually
##' @rdname removePeaks
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testPeakList<-PeakList(analysisName = analysisName(testSpectra),
##' instrument = instrument(testSpectra),
##' nz = nz(testSpectra),
##' calibration = calibration(testSpectra),
##' calibPoints = calibPoints(testSpectra),
##' mz = mz(testSpectra),
##' peakIDs = NULL,
##' peakMzs = NULL)
##' par(mfcol=c(1,2))
##' testPeakList<-addPeaks(testPeakList, mzs = 26:31, width=0.4)
##' plot(testPeakList, mzRange = c(25,32), type = 'l')
##' testPeakList<-removePeaks(testPeakList, mzs = 27)
##' plot(testPeakList, mzRange = c(25,32), type = 'l')
setMethod(removePeaks, signature(object="PeakList", mzs = 'numeric', 
                                 operator = 'missing', limit='missing', 
                                 nLocator="missing"),
          function(object, mzs, operator, limit, nLocator, ...){
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            peaksMZ<-peakMzs(object)[2,]
            
            peakNr<-sapply(1:length(mzs), function(i) which(abs(peaksMZ-mzs[[i]])==min(abs(peaksMZ-mzs[[i]]))))
            message(peakNr)
            peakIDs(object)<-peakIDs(object)[,-peakNr]
            peakMzs(object)<-peakMzs(object)[,-peakNr]
            return(object)
          })

##' @rdname removePeaks
setMethod(removePeaks, signature(object="PeakList", mzs = 'missing', operator="character", 
                                 limit="numeric", nLocator = 'missing'),
          function(object, mzs, operator, limit, nLocator, ...){
            
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            if(operator %in% c(">", ">=", "==", "!=", "<", "<=")){
              peakWidths <- peakWidths(object)
              condition <- paste("which(peakWidths",operator,"limit)")
              removePeaks <- eval(parse(text=condition))
              upperPos <- peakIDs(object)[3,]
              lowerPos <- peakIDs(object)[1,]
              #               nz <- nz(object)
              #               for(peak in removePeaks){
              #                 nz[lowerPos[peak]:upperPos[peak]] <- 0
              #               }
              peakIDs(object)[,removePeaks] <- 0
              peakMzs(object)[,removePeaks] <- 0
              #               nz(object) <- nz
            } else {
              message("Invalid operator.")
            }
            return(object)
          }
)


##' Find single value 'toMatch' in vector 'MatchIn'
##' @param toMatch numeric
##' @param matchIn vector numeric
##' @param twoMatch character 'upper' or 'mean'
##' @return numeric ID of match
findClosestMatch<-function(toMatch,matchIn, twoMatch){
  matched<-which(abs(matchIn-toMatch)==min(abs(matchIn-toMatch)))
  if(length(matched)==2){
    if(twoMatch=='lower'){
      matched<-matched[1]
    } else if(twoMatch=='upper'){
      matched<-matched[2]
    } else if(twoMatch=='mean'){
      matched<-mean(matched)
    }
    
  }
  return(matched)
}

##' peakWidths
##' 
##' peakWidths
##' 
##' Method to return the peakWidth values of all peaks. On plot=TRUE the
##' width values are ploted against the M/z of the corresponding peak.
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testPeakList<-PeakList(analysisName = analysisName(testSpectra),
##' instrument = instrument(testSpectra),
##' nz = nz(testSpectra),
##' calibration = calibration(testSpectra),
##' calibPoints = calibPoints(testSpectra),
##' mz = mz(testSpectra),
##' peakIDs = NULL,
##' peakMzs = NULL)
##' testPeakList<-addPeaks(testPeakList, mzs=26:31, width=0.4)
##' testPeakList<-findPeakWidth(testPeakList, p = 3, n = 199, 
##' span = 100, widthExtLower = 2, widthExtUpper = 2)
##' testPeakList<-peakWidths(testPeakList, plot = FALSE)
##' @rdname peakWidths
setMethod(peakWidths, signature(object="PeakList"),
          function(object, plot=FALSE){
            
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            widths <- as.vector(peakMzs(object)[3,]) - as.vector(peakMzs(object)[1,])
            if(plot) plot(peakMzs(object)[2,],peakWidths(object), main='Peak Widths',
                          xlab='M/z',ylab='Peak width in Da')
            return(widths)
          }
)



##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-reduceSpectrumResolution(object = testSpectra, everyN = 4, mode = 'keep')
##' testSpectra<-smootherSpline(testSpectra, stepsize = 10, spar = 0.3)
##' testSpectra<-smootherGolay(testSpectra, p = 3, n = 5)
##' testSpectra<-peakPick(testSpectra, span = 100)
##' testSpectra<-addFixedWidth(testSpectra, 0.2, 0.2)
##' plot(testSpectra, , mzRange=c(38.5,40.5), type = 'l')
##' @rdname addFixedWidth
setMethod(addFixedWidth, signature(object="PeakList", lowerWidth="numeric", upperWidth="numeric"),
          function(object, lowerWidth, upperWidth){
            
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            mzs <- mz(object)
            peaksCenter <- peakMzs(object)[2,]
            ids <- sapply(peaksCenter, function(peak){
              lower <- peak - lowerWidth
              upper <- peak + upperWidth
              peakWidth <- c(findClosestMatch(lower, mzs, "mean"), findClosestMatch(upper, mzs, "mean"))
              return(peakWidth)
            })
            peakIDs(object)[1,] <- as.integer(ids[1,])
            peakIDs(object)[3,] <- as.integer(ids[2,])
            peakMzs(object)[1,] <- mz(object)[peakIDs(object)[1,]]
            peakMzs(object)[3,] <- mz(object)[peakIDs(object)[3,]]
            return(object)
          }
)

##' @rdname addPeaks
setMethod(addPeaks, signature(object = "PeakList", mzs = "missing", width = "numeric"),
          function(object, mzs, width, ...){
            peakCenters <- manualSelectPeaks(object, ...)
            halfWidth <- width / 2
            if(length(peakCenters) < 1)
              stop("Nothing to add.")
            peakMz <- matrix(0, 3, length(peakCenters))
            peakID <- matrix(0, 3, length(peakCenters))
            
            peakMz[2,] <- peakCenters
            peakMz[1,] <- peakCenters - halfWidth
            peakMz[3,] <- peakCenters + halfWidth
            peakID[2,] <- unlist(lapply(peakMz[2,], function(i) findClosestMatch(i,mz(object),'upper')))
            peakID[1,] <- unlist(lapply(peakMz[1,], function(i) findClosestMatch(i,mz(object),'lower')))
            peakID[3,] <- unlist(lapply(peakMz[3,], function(i) findClosestMatch(i,mz(object),'upper')))
            
            if(is.null(peakIDs(object))){
              peakMzs(object)<-peakMz
              peakIDs(object)<-peakID
            }else{
              peakMzs(object) <- cbind(peakMzs(object), peakMz)
              peakIDs(object) <- cbind(peakIDs(object), peakID)
            }
           plot(object, ...)
           return(object)
          }
)

##' @rdname addPeaks
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testPeakList<-PeakList(analysisName = analysisName(testSpectra),
##' instrument = instrument(testSpectra),
##' nz = nz(testSpectra),
##' calibration = calibration(testSpectra),
##' calibPoints = calibPoints(testSpectra),
##' mz = mz(testSpectra),
##' peakIDs = NULL,
##' peakMzs = NULL)
##' par(mfcol=c(1,2))
##' plot(testPeakList, mzRange=c(25,32), type = 'l')
##' testPeakList<-addPeaks(testPeakList, mzs=26:31, width=0.4)
##' plot(testPeakList, mzRange=c(25,32), type = 'l')
setMethod(addPeaks, signature(object = "PeakList", mzs = "numeric", width = "numeric"),
          function(object, mzs, width, ...){
            peakCenters <- mzs
            halfWidth <- width / 2
            if(length(peakCenters) < 1)
              stop("Nothing to add.")
            peakMz <- matrix(0, 3, length(peakCenters))
            peakID <- matrix(0, 3, length(peakCenters))
            peakMz[2,] <- peakCenters
            peakMz[1,] <- peakCenters - halfWidth
            peakMz[3,] <- peakCenters + halfWidth
            peakID[2,] <- unlist(lapply(peakMz[2,], function(i) findClosestMatch(i,mz(object),'mean')))
            peakID[1,] <- unlist(lapply(peakMz[1,], function(i) findClosestMatch(i,mz(object),'lower')))
            peakID[3,] <- unlist(lapply(peakMz[3,], function(i) findClosestMatch(i,mz(object),'upper')))
            
            if(is.null(peakIDs(object))){
              peakMzs(object)<-peakMz
              peakIDs(object)<-peakID
            }else{
              peakMzs(object) <- cbind(peakMzs(object), peakMz)
              peakIDs(object) <- cbind(peakIDs(object), peakID)
            }
          return(object)
          }
)


##' @rdname changePeakWidth
setMethod(changePeakWidth, signature(object="PeakList", selectMz = 'missing', 
                                     lowerWidth = 'missing', upperWidth = 'missing'),
          function(object, selectMz, lowerWidth, upperWidth, ...){
            
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            choosen <- manualSelectPeaks(object, 1, ...)
            if(length(choosen) < 1)
              stop("Choose 1 nearest peak.")
            mzs <- peakMzs(object)[2,]
            peakIndex <- findClosestMatch(choosen, mzs, "mean")
            peak <- peakMzs(object)[, peakIndex]
            
            mzRange <- list(...)$mzRange
            if(!is.null(mzRange)){
              low.mz<-min(which(mz(object) >= mzRange[1]))
              high.mz<-max(which(mz(object) <= mzRange[2]))
            } else {
              low.mz <- 1
              high.mz <- length(mz(object))
            }
            toploty<-as.vector(nz(object))[low.mz:high.mz]
            toplotx<-mz(object)[low.mz:high.mz]
            
            plot(toplotx,toploty, type='l',ylab='ioncounts',xlab='M/z')
            abline(v=peak[2], lwd=3)
            abline(v=peak[-2])
            newWidths <- locator(2)$x
            if(length(newWidths) < 2)
              stop("Choose 2 widths.")
            newWidths <- sort(newWidths, decreasing=FALSE)
            peakMzs(object)[1,peakIndex] <- newWidths[1]
            peakMzs(object)[3,peakIndex] <- newWidths[2]
            
            
            
            abline(v=newWidths, col=2)
            return(object)
          }
)


##' @rdname changePeakWidth
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testPeakList<-PeakList(analysisName = analysisName(testSpectra),
##' instrument = instrument(testSpectra),
##' nz = nz(testSpectra),
##' calibration = calibration(testSpectra),
##' calibPoints = calibPoints(testSpectra),
##' mz = mz(testSpectra),
##' peakIDs = NULL,
##' peakMzs = NULL)
##' par(mfcol=c(1,2))
##' testPeakList<-addPeaks(testPeakList, mzs=26:31, width=0.4)
##' peakWidths(testPeakList)
##' testPeakList<-changePeakWidth(testPeakList, selectMz = 27, lowerWidth = 0.2, upperWidth = 0.3)
##' peakWidths(testPeakList)
setMethod(changePeakWidth, signature(object="PeakList", selectMz = 'numeric', 
                                     lowerWidth = 'numeric', upperWidth = 'numeric'),
          function(object, selectMz, lowerWidth, upperWidth, ...){
            
            if(is.null(peakIDs(object)))
              stop('There are currently no peaks!')
            
            choosen <- selectMz
            
            mzs <- peakMzs(object)[2,]
            peakIndex <- findClosestMatch(choosen, mzs, "mean")
            peak <- peakMzs(object)[2, peakIndex]
            
            newWidths <- c(peak-lowerWidth, peak+upperWidth)
            newWidths <- sort(newWidths, decreasing=FALSE)
            peakMzs(object)[1,peakIndex] <- newWidths[1]
            peakMzs(object)[3,peakIndex] <- newWidths[2]
            
            
            peakID<-unlist(lapply(peakMzs(object)[1, peakIndex], function(i) findClosestMatch(i,mz(object), 'lower')))
            peakIDs(object)[1, peakIndex] <- peakID
                                   
            peakID<-unlist(lapply(peakMzs(object)[3, peakIndex], function(i) findClosestMatch(i,mz(object), 'upper')))                                                              
            peakIDs(object)[3, peakIndex] <- peakID
                                   
            return(object)
         
          }
)



##' This method is base method for plotting and manual select data
##' @param object object of type PeakList
##' @param n numeric
##' @param ... additional args
##' @return numeric x coordinates
manualSelectPeaks <- function(object, n=512, ...){
  plot(object, ...)
  selected <- locator(n)
  return(selected$x)
}


