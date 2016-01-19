##' \code{MassSpectra} class constructor
##'
##' \code{MassSpectra} is also the call to class constructor. It is used
##' for importing high-resolution mass spectra from raw data.
##'
##' \code{MassSpectra} is also the call to class constructor. It is used
##' for importing high-resolution mass spectra from raw data.
##' @export
##' @rdname MassSpectra
##' @param select \code{character}, 'ulvacraw', 'iontofgrd', 'dummy'
##' @param analysisName \code{character}, the (file)name of the dataset
##' @param ... additional args
##' @return object of class MassSpectra
##' @examples 
##' \dontrun{
##' ## access rawdata in tofsims package
##' library(tofsimsData)
##' importFile<-system.file("rawdata", "trift_test_001.RAW", package = "tofsimsData")
##' MassSpectra('ulvacraw', importFile)
##' }
##' ## create dummy MassSpectra object
##' MassSpectra('dummy')
MassSpectra <- function(select = c("ulvacraw", "iontofgrd", "dummy"), 
                        analysisName, ...) {
    
    #selection <- c("ulvacraw", "iontofgrd", "dummy")
    
    #choice <- which(selection == select)
    choice <- match.arg(select)
    
    switch(choice,
           ulvacraw = {
             import <- cReadRawPhi(analysisName, mode = "spectra", ...)
             
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             calibration <- import$calibration
             mz <- import$mz
           },
           
           iontofgrd = {
             import <- import.raw(analysisName, mode = "spectra", ...)
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             calibration <- import$calibration
             mz <- import$mz
           },
           
           dummy = {
             analysisName  <- "dummyAnalysisName"
             instrument    <- "dummyInstrument"
             nz            <- matrix(round(runif(256000,0,9)),256,1000)
             calibration   <- data.frame(intercept = 20000, slope = 20000)
             mz            <- seq(0.1,100,0.1)
           }
   
    )
    
#     if (length(choice) == 0) 
#         stop(select, " does not match with any available method.\n \n       
#              Please choose, \"ulvacraw\", or \"iontofgrd\".")
#     
#     if (choice == 1) {
#         import <- cReadRawPhi(analysisName, mode = "spectra", ...)
#         
#         analysisName <- import$analysisName
#         instrument <- import$instrument
#         nz <- import$nz
#         calibration <- import$calibration
#         mz <- import$mz
#     }
#     
#     
#     if (choice == 2) {
#         import <- import.raw(analysisName, mode = "spectra", ...)
#         analysisName <- import$analysisName
#         instrument <- import$instrument
#         nz <- import$nz
#         calibration <- import$calibration
#         mz <- import$mz
#     }
#     
#     if (choice == 3) {
#         analysisName  <- "dummyAnalysisName"
#         instrument    <- "dummyInstrument"
#         nz            <- matrix(round(runif(256000,0,9)),256,1000)
#         calibration   <- data.frame(intercept = 20000, slope = 20000)
#         mz            <- seq(0.1,100,0.1)
#     }
    
    new("MassSpectra", 
        analysisName = as.character(basename(analysisName)),
        instrument   = as.character(instrument),
        nz           = as.matrix(nz),
        calibration  = as.vector(calibration), 
        mz           = as.numeric(mz))
}


##' method definition 'dim' for 'MassSpectra'
##' dim is a primitive
##' @param x object object of type MassSpectra
##' @return numeric value
setMethod("dim", signature(x = "MassSpectra"), function(x) {
    dim(nz(x))
})

##' method definition 'ndim' on 'MassSpectra'
##' @param object object of type MassSpectra
##' @return numeric value
setMethod(ndim, signature(object = "MassSpectra"), 
    function(object) dim(nz(object))[1])

##' method definition 'zdim' on 'MassSpectra'
##' @param object object of class MassSpectra
##' @return numeric value
setMethod(zdim, signature(object = "MassSpectra"), 
    function(object) dim(nz(object))[2])

##' method defining \code{show()} for the \code{MassSpectra} class
##' show has a generic by default
##' @param object object of class MassSpectra
##' @return data.frame character
setMethod("show", signature(object = "MassSpectra"), 
    function(object) {
        if (length(object@analysis) > 0) 
            analysis <- paste(unlist(lapply(1:length(object@analysis), 
                function(i) 
                    paste(i, class(object@analysis[[i]]), sep = "-")
                )), 
                collapse = ",") 
        else 
            analysis <- "none"
        print(data.frame(
            analysisName = analysisName(object), 
            instrument = instrument(object), 
            no.spectra = format(ndim(object), big.mark = "'"), 
            spectra.points = format(zdim(object), big.mark = "'"), 
            `total ion count` = format(sum(nz(object)), big.mark = "'"), 
            analysis = analysis))
    })

##' Method \code{plot()} for \code{MassSpectra}
##' 
##' Method defining \code{plot()} for the \code{MassSpectra} class
##' plot has no generic by default
##' 
##' The output of this method is adapted for plotting mass spectra. Uncalibrated
##' data is plotted as xy plot while uncalibrated data is plotted as barplot. 
##' The parameter \code{mzRange} allows choosing the plot range directly 
##' according to the mz number (when calibrated). The argument \code{lineplot}, 
##' TRUE by default, allows to switch between line and barplot.
##' @param x object of type MassSpectra
##' @param y missing 
##' @param ... further args
##' @param mzRange vector or lenght two, indicating the mz range to be plotted
##' @param normalize boolean should the mass spectra be normalized
##' @return plot of mass spectra
##' @rdname plot
##' @examples 
##' ## plot method for MassSpectra objects
##' library(tofsimsData)
##' data(tofsimsData)
##' plot(testSpectra, mzRange=c(1,300),type='l')
setMethod("plot", signature(x = "MassSpectra", y = "missing"), 
    function(x, y, ..., mzRange = c(0, 200), normalize = FALSE) {
        low.mz <- 1
        high.mz <- length(mz(x))
        
        if (!is.null(mzRange)) {
            low.mz <- min(which(mz(x) >= mzRange[1]))
            high.mz <- max(which(mz(x) <= mzRange[2]))
        }
        
        if (ndim(x) == 1) {
            if (normalize) {
                nzNorm <- matrix(100/max(nz(x)) * nz(x), 1)
                toploty <- as.vector(nzNorm)[low.mz:high.mz]
            } else {
                toploty <- as.vector(nz(x))[low.mz:high.mz]
            }
            
            toplotx <- mz(x)[low.mz:high.mz]
            
            plot(toplotx, toploty, ylab = "ioncounts", xlab = "M/z", ...)
        } else {
            # toplotx<-mz(x)[low.mz:high.mz]
            
            if (normalize) {
                makeTIC <- apply(nz(x), 2, mean)
                nzNorm <- matrix(100/max(makeTIC) * makeTIC, 1)
                toploty <- nzNorm[low.mz:high.mz]
            } else {
                toploty <- apply(nz(x), 2, mean)[low.mz:high.mz]
            }
            
            toplotx <- mz(x)[low.mz:high.mz]
            plot(toplotx, toploty, ylab = "ioncounts", xlab = "M/z", ...)
        }
    })




##' scale
##' 
##' scale
##' autoscaling method for MassSpectra object. Scaling
##' is along the mass channels. Therefore more than one
##' spectra is needed for scaling. 
##' @param x object object of class MassSpectra
##' @param center boolean should data be centered
##' @param scale boolean should data be scaled
##' @return object of class MassSpectra
##' @rdname scale
##' @examples
##' ## autoscaling of dummy image data
##' testImage<-MassImage('dummy')
##' par(mfcol=c(2,2))
##' plot(testImage,type='l')
##' image(testImage)
##' testImage <- scale(testImage)
##' plot(testImage,type='l')
##' image(testImage)
##' \dontrun{
##' ## autoscaling of real spectral data
##' library(tofsimsData)
##' data(tofsimsData)
##' par(mfcol=c(2,2))
##' plot(testImage,type='l')
##' image(testImage)
##' testImage <- scale(testImage)
##' plot(testImage,type='l')
##' image(testImage)}
setMethod(scale, signature(x = "MassSpectra"), 
          function(x, center = TRUE, scale = TRUE) {
            if(ndim(x)==1) stop('can not autoscale a single spectra')
            scaled.dataset <- scale(nz(x), center, scale)
            colnames(scaled.dataset) <- mz(x)
            x@nz <- scaled.dataset
            return(x)
})



##' Possion scaling for data matrices.
##'
##' Possion scaling is proposed as the method of choice for ToF-SIMS data
##' see Keenan and Kotula 2004. This implementation was done according to
##' a description in \code{Multivariate Analysis of SIMS spectra in
##' ToF-SIMS: Materials Analysis by Mass Spectrometry, Vickerman and
##' Briggs 2013} and the \code{eigenvector wiki}. The offset is described in
##' the \code{eigenvector wiki}.
##' @export
##' @param object object of class MassSpectra
##' @param offset numeric value for poisson scaling
##' @param ... further args
##' @return object of class MassSpectra
##' @title Poisson Scaling
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
##' @rdname poissonScaling
##' @examples 
##' ## poisson scaling of MassSpectra objects
##' testImage <- MassImage('dummy')
##' testImage <- poissonScaling(testImage)
##' \dontrun{
##' # poission scaling on real data
##' library(tofsimsData)
##' data(tofsimsData)
##' par(mfcol=c(2,2))
##' plot(testImage,type='l')
##' image(testImage)
##' testImage <- poissonScaling(testImage)
##' plot(testImage,type='l')
##' image(testImage)
##' }
setMethod(poissonScaling, signature(object = "MassSpectra"), 
    function(object, offset = 1, ...) {
        dataset <- nz(object)
        message("calculating offset")
        off.100 <- max(apply(dataset, 2, mean))
        dataset <- dataset + (off.100 * (offset/100))
        message("scaling first dimension...")
        sqrt.mean.row <- sqrt(apply(dataset, 1, mean))
        message("scaling second dimension...")
        sqrt.mean.col <- sqrt(apply(dataset, 2, mean))
        message("combining...")
        row.matrix <- matrix(rep(
            sqrt.mean.row, ncol(dataset)), 
            nrow(dataset), 
            ncol(dataset))
        col.matrix <- matrix(rep(
            sqrt.mean.col, nrow(dataset)), 
            nrow(dataset), 
            ncol(dataset), 
            byrow = TRUE)
        scaled.dataset <- dataset/col.matrix/row.matrix
        colnames(scaled.dataset) <- mz(object)
        
      
        nz(object) <- scaled.dataset
        return(object)
       
    })

##' reduceSpectrumResolution
##' 
##' reduceSpectrumResolution
##' 
##' The method reduceSpectrumResolution for MassSpectra is used
##' sometimes for performance reasons.
##' @param object object of class MassSpectra
##' @param everyN numeric act on every nth spectra point
##' @param mode character 'remove' or 'keep'
##' @return object of class MassSpectra
##' @rdname reduceSpectrumResolution
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' par(mfcol=c(1,2))
##' plot(testSpectra,mzRange = c(40,50),type='l')
##' testSpectra <- reduceSpectrumResolution(object = testSpectra, everyN = 2, mode = 'remove')
##' plot(testSpectra, mzRange = c(40,50), type='l')
setMethod(reduceSpectrumResolution, signature(object = "MassSpectra"), 
    function(object, everyN, mode) {
        if(everyN < 0)
            stop("Invalid everyN parameter. This should be > 0.")
        if (mode == "remove") {
            reducedResolution <- nz(object)[, -(seq(1, zdim(object), everyN))]
            reducedResolutionMzs <- mz(object)[-(seq(, zdim(object), everyN))]
        } else if (mode == "keep") {
            reducedResolution <- nz(object)[, seq(1, zdim(object), everyN)]
            reducedResolutionMzs <- mz(object)[seq(, zdim(object), everyN)]
        } else {
            stop("Unsupported mode.")
        }
        # colNamesReducedResolution<-names(reducedResolution)
        matReducedResolution <- matrix(reducedResolution, 1)
        # colnames(matReducedResolution)<-colNamesReducedResolution
        object@nz <- matReducedResolution
        object@mz <- reducedResolutionMzs
        return(object)
    })


##' Method \code{points()} for \code{MassSpectra}
##' 
##' Method defining \code{points()} for the \code{MassSpectra} class
##' points has no generic by default
##' 
##' This function can be used to visualize several spectra in the same
##' plot.
##' @param x vector with mz for mass spectra plot
##' @param y vector with ion counts for mass spectra plot
##' @param ... additional args
##' @param mzRange vector of length 2, indicating the mz range to be plotted
##' @param normalize boolean should the mass spectra be normalized
##' @return graphic output
##' @rdname points
##' @examples
##' library(tofsimsData) 
##' data("tofsimsData")
##' plot(testImage, type='l', normalize = TRUE, col = 'blue')
##' points(testSpectra, type = 'l', normalize = TRUE, col = 'red')
setMethod("points", signature(x = "MassSpectra"), 
          function(x, y, ..., mzRange = c(0, 200), normalize = FALSE) {
    low.mz <- 1
    high.mz <- length(mz(x))
    
    if (!is.null(mzRange)) {
        low.mz <- min(which(mz(x) >= mzRange[1]))
        high.mz <- max(which(mz(x) <= mzRange[2]))
    }
    
    if (ndim(x) == 1) {
        if (normalize) {
            nzNorm <- matrix(100/max(nz(x)) * nz(x), 1)
            toploty <- as.vector(nzNorm)[low.mz:high.mz]
        } else {
            toploty <- as.vector(nz(x))[low.mz:high.mz]
        }
        
        toplotx <- mz(x)[low.mz:high.mz]
        points(toplotx, toploty, ...)
        
    } else {
        toplotx <- mz(x)[low.mz:high.mz]
        if (normalize) {
            makeTIC <- apply(nz(x), 2, mean)
            nzNorm <- matrix(100/max(makeTIC) * makeTIC, 1)
            toploty <- nzNorm[low.mz:high.mz]
        } else {
            toploty <- apply(nz(x), 2, mean)[low.mz:high.mz]
        }
        
        toplotx <- mz(x)[low.mz:high.mz]
        points(toplotx, toploty, ...)
        
    }
})

##' Overlay plot of Mass Spectra
##' 
##' This function takes as input a list with objects of type MassSpectra. The
##' easiest way to obtain the input data, is to use mclapply from the parallel 
##' package.
##' @param objectList list with object of type MassSpectra
##' @param ... additional args
##' @param type character type of plot, usually 'l'
##' @param mzRange vector numeric lower and upper range for plotting the spectra
##' @param PeakListObj object a PeakList object can be provided to plot peaks
##' @param cex.legend numeric text size
##' @return graphical output
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
##' @rdname overlayPlot
##' @examples
##' library(tofsimsData)
##' data('tofsimsData')
##' overlayPlot(list(testImage, testSpectra))
setMethod(overlayPlot, signature(object = "list"), 
    function(objectList, ..., 
             type = "l", 
             mzRange = c(1, 200), 
             PeakListObj = NULL, 
             cex.legend = 0.5) {
        
        ### check for wether the list objects are of correct
        ### type
        if (sum(unlist(lapply(1:length(objectList), 
            function(i) !validMassSpectraObject(objectList[[i]])))) > 0) 
            stop("invalid object")
        
        requireNamespace('RColorBrewer')
        numberOfSpectra <- length(objectList)
        ndims <- lapply(1:length(objectList), 
                        function(i) ndim(objectList[[i]]))
        
        if (length(which(ndims > 1)) > 0) {
            ndimsLarger <- which(ndims > 1)
            makeTICs <- lapply(ndimsLarger, 
                               function(i) apply(nz(objectList[[i]]), 2, sum))
            counter <- 1
            for (i in ndimsLarger) {
                tempMzs <- mz(objectList[[i]])
                nz(objectList[[i]]) <- matrix(makeTICs[[counter]], 1)
                mz(objectList[[i]]) <- tempMzs
                counter <- counter + 1
            }
            
        }
        
        
        normalized <- lapply(1:numberOfSpectra, 
                             function(i) 100/max(nz(objectList[[i]])) * nz(objectList[[i]]))
        mzList <- lapply(1:numberOfSpectra, 
                         function(i) mz(objectList[[i]]))
        low.mz <- lapply(1:numberOfSpectra, 
                         function(i) min(which(mzList[[i]] >= mzRange[1])))
        high.mz <- lapply(1:numberOfSpectra, 
                          function(i) max(which(mzList[[i]] <= mzRange[2])))
        
        
        
        find.ylim <- max(unlist(lapply(1:numberOfSpectra, 
            function(i) max(normalized[[i]][low.mz[[i]]:high.mz[[i]]]))))
        
        suppressWarnings(colorlist <- matrix(
            RColorBrewer::brewer.pal(name = "Set1", n = 9), 1, length(objectList)))
        # colorlist<-replicate(length(objectList),
        # rgb(round(runif(1,0,255),0),round(runif(1,0,255),0),
        # round(runif(1,0,255),0),maxColorValue=255))
        
        if (sum(unlist(lapply(1:length(objectList), 
            function(i) !validMassSpectraObject(objectList[[i]])))) > 0) 
                stop("invalid object")
        
        plot(objectList[[1]], col = colorlist[1], mzRange = mzRange, 
            ylim = c(0, (find.ylim)), type = type, normalize = TRUE, ...)
        if (numberOfSpectra > 1) {
            for (i in 2:numberOfSpectra) {
                points(objectList[[i]], 
                        col = colorlist[i], 
                        mzRange = mzRange, type = type, normalize = TRUE, ...)
            }
        }
        
        ### if provided print peaks
        if (class(PeakListObj) == "PeakList") {
            peaks <- peakMzs(PeakListObj)[2, ]
            abline(v = peaks)
            
            if (dim(peakMzs(PeakListObj))[2] > 1) {
                lower <- as.numeric(peakMzs(PeakListObj)[1,])
                upper <- as.numeric(peakMzs(PeakListObj)[3,])
                abline(v = lower, col = RColorBrewer::brewer.pal("Set1", n = 9))
                abline(v = upper, col = RColorBrewer::brewer.pal("Set1", n = 9))
            }
        }
        
        
        ### print legend
        legend.low.mz <- mzList[[1]][low.mz[[1]]]
        legend.names <- unlist(lapply(1:length(objectList), 
            function(i) analysisName(objectList[[i]])))
        legend(legend.low.mz, find.ylim, legend.names, 
            lty = c(1, 1), col = colorlist, cex = cex.legend)
        
        
    })

##' Method makeTIC for MassSpectra Class
##' 
##' Method makeTIC sums up all Mass Spectra in the called Mass Spectra object
##' @param object object of class MassSpectra
##' @return object of class MassSpectra with just one spectra, the TIC
setMethod(makeTIC, signature(object = c("MassSpectra")), 
    function(object) {
        if (zdim(object) > 1) {
            temp <- matrix(apply(nz(object), 2, sum), 1)
            colnames(temp) <- mz(object)
            nz(object) <- temp
            return(object)
        } else {
            message("Already just one spectra present.")
        }
    })


##' Method smootherGolay for MassSpectra class
##' @param object object of class MassSpectra
##' @param p numeric parameter for savitzky-golay filter
##' @param n numeric parameter for savitzky-golay filter
##' @param ... additional args
##' @return object of class MassSpectra with smoothed TIC
##' @rdname smootherGolay
##' @import signal
##' @examples 
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectraSmooth <- smootherGolay(testSpectra, p = 3, n = 9)
##' overlayPlot(list(testSpectra, testSpectraSmooth), mzRange = c(38.5, 40.5), type = 'l')
setMethod(smootherGolay, signature(object = "MassSpectra"), 
    function(object, p = 3, n = 5, ...) {
        if (zdim(object) > 1) {
            requireNamespace('signal')
            sgolay.filter <- signal::sgolay(p = p, 
                                            n = n, 
                                            m = 0, 
                                            ts = 1)
            requireNamespace('signal')
            filtered <- matrix(signal::filter(filt = sgolay.filter, 
                                      x = as.vector(nz(object))), 1)
            filtered[which(is.na(filtered))]<-0
            #colnames(filtered) <- mz(object)
            nz(object) <- filtered
            return(object)
        } else {
            message("More than one spectra. Please use for 'makeTIC' method.")
        }
    })


##' method smootherSpline for TIC
##' @return object of class MassSpectra
##' @rdname smootherSpline
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' testSpectraSmooth <- smootherSpline(testSpectra)
##' overlayPlot(list(testSpectra, testSpectraSmooth), mzRange = c(38.5, 40.5), type = 'l')
setMethod(smootherSpline, signature(object = "MassSpectra"), 
    function(object, stepsize = 5, spar = 0.3, ...) {
        if (zdim(object)) {
            
            lowest <- as.numeric(head(mz(object), 1))
            highest <- as.numeric(tail(mz(object), 1))
            stepsInitial <- seq(lowest, highest, stepsize)
            
            ### matching the initial steps to actual m/z's
            ### available
            steps <- NULL
            for (i in stepsInitial) {
                steps <- c(steps, 
                           which(abs(mz(object) - i) == 
                                     min(abs(mz(object) - i))))
            }
            
            ### adding the last point again to the steps
            steps <- c(steps, highest)
            
            ### remove steps with less than 4 datapoints
            if (!length(which(diff(steps) < 4) + 1) == 0) {
                steps <- steps[-(which(diff(steps) < 4) + 1)]
            }
            
            ### actual smoothing process, looping through all
            ### steps
            smoothTIC <- NULL
            
            for (i in 1:(length(steps) - 1)) {
                splinefit <- smooth.spline(
                    x = mz(object)[steps[i]:(steps[i + 1] - 1)], 
                    y = as.vector(nz(object))[steps[i]:(steps[i + 1] - 1)], 
                    spar = spar)
                buildSegment <- matrix(splinefit$y, 1)
                colnames(buildSegment) <- splinefit$x
                smoothTIC <- c(smoothTIC, buildSegment)
            }
            
            smoothTIC <- matrix(smoothTIC, 1)
            colnames(smoothTIC) <- mz(object)[1:length(smoothTIC)]
            
            nz(object) <- smoothTIC
            mz(object) <- mz(object)[1:length(smoothTIC)]
            return(object)
        } else {
            message("More than one spectra. Please use for 'makeTIC' method.")
        }
    })

##' method peakPick
##' 
##' method peakPick
##' 
##' Method peakPick for MassSpectra class, works as a constructor for PeakList class.
##' The local min/max detection implementation is adapted from the CRAN package
##' 'ChemometricsWithR'. 
##' @param object object of class MassSpectra
##' @param span numeric parameter for local max/min detection
##' @param ... additional args
##' @return object of class PeakList
##' @rdname peakPick
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra <- reduceSpectrumResolution(object = testSpectra, everyN = 4, mode = 'keep')
##' testSpectra <- smootherSpline(testSpectra, stepsize = 10, spar = 0.3)
##' testSpectra <- smootherGolay(testSpectra, p = 3, n = 5)
##' testSpectra <- peakPick(testSpectra, span = 100)
##' plot(testSpectra, , mzRange=c(38.5,40.5), type = 'l')
setMethod(peakPick, signature(object = "MassSpectra"), 
    function(object, span = 100, ...) {
        if (zdim(object) > 1) {
            low.mz <- 1
            high.mz <- length(mz(object))
            
            x <- as.vector(nz(object))
            
            spanWidth <- span * 2 + 1
            localMax <- spanWidth + 1 - apply(embed(x, spanWidth), 
                                              1, 
                                              which.max)
            localMax[localMax == 1 | localMax == spanWidth] <- NA
            pks <- localMax + 0:(length(localMax) - 1)
            peaks <- matrix(unique(pks[!is.na(pks)]), 1)
            colnames(peaks) <- mz(object)[peaks]
            peakIDs <- peakMzs <- matrix(rep(0, length(peaks) * 3), 3)
            peakIDs[2, ] <- peaks
            peakMzs[2, ] <- mz(object)[peaks]
            # newMzs<-mz(object)[peaks] print(peaks[1:10])
            
            object <- PeakList(
                analysisName = analysisName(object), 
                instrument = instrument(object), nz = nz(object), 
                calibration = calibration(object), 
                calibPoints = calibPoints(object), 
                mz = mz(object), peakIDs = peakIDs, 
                peakMzs = peakMzs)
            return(object)
        } else {
            message("More than one spectra. Please use for 'makeTIC' method.")
        }
    })


##' method getTOFs
##' @param object object of class MassSpectra
##' @return vector numeric with TOF values
##' @rdname getTOFs
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' timeOfFlight <- getTOFs(testSpectra)
##' head(timeOfFlight)
setMethod(getTOFs, signature(object = "MassSpectra"), 
    function(object) {
        calibration(object)$slope * sqrt(mz(object)) + 
            calibration(object)$intercept
    })

##' \code{calibPointNew} is a method to set a new mass calibration point
##' 
##' \code{calibPointNew} is a method to set a new mass calibration point
##' 
##' \code{calibPointNew} ia a method to set a new mass calibration point. When 
##' \code{value} is not provided as arguemnt, the TOF for the chosen \code{mz} 
##' value has to be chosen interactively by mouse. 
##' @param object MassSpectra object
##' @param mz the m/z value to be specified with a TOF value
##' @param reset shall the list of calibration points be reset
##' @param value TOF value to be assigned to mz
##' @return object MassSpectra with added/updated calibration points
##' @rdname calibPointNew
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra <- calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra <- calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' calibPoints(testSpectra)
##' par(mfcol=c(1,2))
##' plot(testSpectra,mzRange=c(38.5,40.5),type='l')
##' testSpectra <- recalibrate(testSpectra)
##' plot(testSpectra, mzRange=c(38.5,40.5), type='l')
setMethod(calibPointNew, signature(object = "MassSpectra", 
    mz = "numeric"), function(object, mz, reset = FALSE, value = NULL) {
    
    if (is.null(value)) {
        plot(object, type = "l", mzRange = c(mz - 2, mz + 2))
        coords <- locator(n = 1)
        abline(v = coords$x, col = "red")
        message(coords$x)
    } else {
        coords <- list(x = value)
    }
    
    mzid <- sapply(1:length(coords$x), 
                   function(i) 
                       which(abs(mz(object) - coords$x[[i]]) == 
                                 min(abs(mz(object) - coords$x[[i]]))))
    if (reset) {
        calibPoints(object) <- data.frame(
            mz = as.numeric(mz), 
            TOF = as.numeric(getTOFs(object)[mzid]))
        return(object)
    } else {
        calibPoints(object) <- rbind(
            calibPoints(object), 
            data.frame(mz = as.numeric(mz), 
                       TOF = as.numeric(getTOFs(object)[mzid])))
        return(object)
    }
})

##' method recalibrate
##' @param object object of class MassSpctra
##' @return object of class MassSpectra, recalibrated mass values
##' @rdname recalibrate
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra <- calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra <- calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' calibPoints(testSpectra)
##' par(mfcol=c(1,2))
##' plot(testSpectra,mzRange=c(38.5,40.5),type='l')
##' testSpectra <- recalibrate(testSpectra)
##' plot(testSpectra, mzRange=c(38.5,40.5), type='l')
setMethod(recalibrate, signature(object = "MassSpectra"), 
    function(object) {
        message("Calculating recalibration.")
        dataPoints <- data.frame(mzs = sqrt(calibPoints(object)$mz), 
            TOFs = calibPoints(object)$TOF)
        dataFull <- data.frame(TOFs = getTOFs(object))
        new.calib <- lm(mzs ~ TOFs, dataPoints)
        predicted <- predict(new.calib, dataFull)^2
        message("Writing recalibration into variables.")
        
        mz(object) <- predicted
        reverseCalib <- lm(TOFs ~ mzs, dataPoints)
        calibration(object) <- data.frame(
                intercept = as.numeric(reverseCalib$coefficients[1]), 
                slope = as.numeric(reverseCalib$coefficients[2]))
        return(object)
        
    })


##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra <- calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra <- calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra <- recalibrate(testSpectra)
##' testSpectra <- unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' plot(testSpectra, mzRange = c(1,200), type = 'l')
##' @rdname unitMassPeaks
setMethod(unitMassPeaks, signature(object = "MassSpectra", 
                                   mzRange = "numeric", 
                                   widthAt = "numeric"), 
    function(object, mzRange, widthAt, factor, upper = NULL, 
             lower = NULL, ...) {
    if (!length(mzRange) == 2) 
        stop("mzRange needs two values, lower and upper peak range")
    
    ### choose upper and lower width of the two chosen
    ### peaks
    widthCoords <- NULL
    for (i in c(1, 2)) {
        
        if (is.null(upper) || is.null(lower)) {
            plot(object, type = "l", 
                 mzRange = c(widthAt[i] - 1, widthAt[i] + 1))
            abline(v = widthAt[i])
            message("Please choose upper and lower peak limit.")
            widthCoords <- c(widthCoords, locator(n = 2)$x)
        } else {
            widthCoords <- c(widthCoords, lower, upper)
        }
    }
    
    widthIDs <- unlist(lapply(widthCoords, 
                              function(i) findClosestMatch(i, mz(object))))
    
    ### make peak variables
    peakIDs <- peakMzs <- matrix(rep(0, length(
        seq(mzRange[1], mzRange[2], 1)) * 3), 3)
    
    message("Looking up unit mass id's.")
    peakIDs[2, ] <- unlist(lapply(seq(mzRange[1], mzRange[2], 1), 
                                  function(i) findClosestMatch(i, mz(object))))
    peakMzs[2, ] <- mz(object)[peakIDs[2, ]]
    
    message("Calculating peak widths.")
    lowerWidth <- abs(widthCoords[1] - widthCoords[2])
    upperWidth <- abs(widthCoords[3] - widthCoords[4])
    modelSet <- data.frame(width = c(lowerWidth, upperWidth), mz = mzRange)
    modelWidth <- lm(width ~ mz, modelSet)
    calculatedWidth <- predict(modelWidth, 
                               data.frame(mz = seq(mzRange[1], mzRange[2], 1)))
    lowerWidths <- seq(mzRange[1], mzRange[2], 1) - calculatedWidth * factor[1]
    upperWidths <- seq(mzRange[1], mzRange[2], 1) + calculatedWidth * factor[2]
    
    message("Looking up lower peak limit id's.")
    peakIDs[1, ] <- unlist(lapply(lowerWidths, 
                                  function(i) findClosestMatch(i, mz(object))))
    message("Looking up upper peak limit id's.")
    peakIDs[3, ] <- unlist(lapply(upperWidths, 
                                  function(i) findClosestMatch(i, mz(object)))
                           )
    peakMzs[c(1, 3), ] <- mz(object)[peakIDs[c(1, 3), ]]
    
    object <- PeakList(
        analysisName = analysisName(object), 
        instrument = instrument(object), nz = nz(object), 
        calibration = calibration(object), calibPoints = calibPoints(object), 
        mz = mz(object), peakIDs = peakIDs, peakMzs = peakMzs)
    return(object)
    
    if (is.null(upper) || is.null(lower)) 
        return(widthCoords)
})



##' @rdname bwApply
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage <- PCAnalysis(testImage, nComp = 2)
##' library(EBImage)
##' mask<-thresh(imageMatrix(analysis(testImage,noAccess = 1),comp = 1), w = 15, h = 15)
##' #inverse of mask
##' mask <- (mask-1)^2
##' par(mfcol=c(1,2), oma=c(0,0,0,0), mar=c(0,0,0,0))
##' image(testImage)
##' image(bwApply(testImage, mask))
setMethod(bwApply, signature("MassSpectra", "matrix"), 
    function(object, bwMatrix) {
        bwDim <- dim(bwMatrix)
        if (bwDim[1] != xdim(object) || bwDim[2] != 
            ydim(object)) 
            stop("B/W matrix doesn't have the same dimensions with input 
                 object.")
        bwVector <- as.vector(bwMatrix)
        nz(object)[which(bwVector == 1), ] <- 0
        return(object)
    })
