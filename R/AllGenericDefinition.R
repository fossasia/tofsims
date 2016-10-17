######
### A
######

##' Generic method to add/update peak width
##' 
##' This method will update current upper/lower width for all
##' peaks
##' @export
##' @param object PeakList object
##' @param lowerWidth numeric
##' @param upperWidth numeric
##' @return object PeakList with updated/new peak widths
##' @rdname addFixedWidth
setGeneric("addFixedWidth",
           function(object, lowerWidth, upperWidth) standardGeneric("addFixedWidth"))

##' generic method to add peaks
##' 
##' This method will allow user to plot and add peaks manually.
##' This method will take all parameters of PeakList plot method.
##' @export
##' @param object PeakList object
##' @param mzs numeric vector M/z's where peaks shall be added 
##' @param width fixed value to add (m/z)
##' @param ... further args
##' @return object updated PeakList object
##' @rdname addPeaks
setGeneric("addPeaks",
           function(object, mzs, width, ...) standardGeneric("addPeaks"))

##' @title \code{analysis}, slot of \code{MassSpectra} class objects
##' @description \code{analysis}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @param noAccess numeric access number to analysis slot
##' @param ... additional args
##' @return summary or content of analysis slot
##' @rdname analysis
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage, nComp = 3)
##' ## obtain summary of analysis slot content
##' analysis(testImage)
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{nz}} \code{\link{analysisName}} \code{\link{instrument}}
##' \code{\link{calibPoints}} \code{\link{calibration}}
setGeneric("analysis",
           function(object, noAccess, ...) standardGeneric("analysis"))

##' @export
##' @param value object to be put in analysis slot
##' @rdname analysis
setGeneric("analysis<-",
           function(object, value) standardGeneric("analysis<-"))

##' @title \code{analysisName}, slot of \code{MassSpectra} class objects
##' @description \code{analysisName}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @param ... further args
##' @return content of analysisName slot
##' @rdname analysisName
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' ## access name of analysis
##' analysisName(testSpectra)
##' ## replace name of analysis
##' analysisName(testSpectra) <- 'sample001_pos001_settings_default'
##' analysisName(testSpectra)
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{analysis}} \code{\link{nz}} \code{\link{instrument}}
##' \code{\link{calibPoints}} \code{\link{calibration}}
setGeneric("analysisName", 
           function(object, ...) standardGeneric("analysisName"))

##' @export
##' @param value character replacement value for slot analysisName
##' @rdname analysisName
setGeneric("analysisName<-",
           function(object, value) standardGeneric("analysisName<-"))






######
### B
######

##' generic accessor method baseObject
##' @param object helper for prcomp and princomp wrappers
##' @return baseObject
setGeneric("baseObject",
           function(object) standardGeneric("baseObject"))


##' binning
##' 
##' binning
##' 
##' bining is used to reduce the resolution/size of MassImage objects. 
##' Optionally \code{mclapply} from the parallel package is used to speed up 
##' processing time.
##' @export
##' @param object object of class MassImage
##' @param binningFactor numeric factor for binning (2, 4, etc)
##' @param ... additional args
##' @return binned object of class MassImage
##' @rdname binning
setGeneric("binning",
           function(object, binningFactor, ...) standardGeneric("binning"))

##' bwApply
##' 
##' bwApply allow to get new object from a black / white matrix
##' All NZs at black positions will be taken
##' @export
##' @param object object of class MassImage
##' @param bwMatrix matrix with boolean or numeric 1 and 0
##' @return object of class MassImage multiplied with B/W matrix
##' @rdname bwApply
setGeneric("bwApply",
          function(object, bwMatrix) standardGeneric("bwApply"))



######
### C
######

##' @title \code{calibration}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @return content of calibration slot
##' @rdname calibration
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' ## access calibration slot
##' calibration(testSpectra)
##' ## replacing the values of the 'calibration' slot is
##' ## possible but it makes at the moment no sense as it
##' ## doesn't change the actual mass calibration. The 
##' ## 'calibration' slot is just used to store the values
##' ## while 'recalibration' uses the values from 
##' ## 'calibPoints' slot.
##' calibration(testSpectra) <- data.frame(intercept = 21420, slope = 20480)
##' calibration(testSpectra)
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{analysis}} \code{\link{analysisName}} \code{\link{instrument}}
##' \code{\link{calibPoints}} \code{\link{nz}}
setGeneric("calibration",
           function(object) standardGeneric("calibration"))

##' Generic setter for slot calibration<- 
##' @export
##' @param value data.frame with replacement values for calibration slot
##' @rdname calibration
setGeneric("calibration<-",
           function(object,value) standardGeneric("calibration<-"))

##' Generic method calibPointNew that modifies slot calibPoints
##' @export
##' @return call by reference, hence MassSpectra object with new calib point
setGeneric("calibPointNew",
           function(object, mz, reset=FALSE, value=NULL) standardGeneric("calibPointNew"))

##' @title \code{calibPoints}, slot of \code{MassSpectra} class objects
##' @description  \code{calibPoints}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @return contents of slot calibPoints
##' @rdname calibPoints
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-calibPointNew(testSpectra, mz = 15, value = 15.0232)
##' testSpectra<-calibPointNew(testSpectra, mz = 181, value = 181.0569)
##' ## access 'calibPoint' slot of 'MassSpectra' object
##' calibPoints(testSpectra)
##' ## replacing values in the 'calibPoint' slot
##' calibPoints(testSpectra)[2,2]<-297000
##' calibPoints(testSpectra)
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{analysis}} \code{\link{analysisName}} \code{\link{instrument}}
##' \code{\link{nz}} \code{\link{calibration}}
setGeneric("calibPoints",
           function(object) standardGeneric("calibPoints"))

##' @export
##' @param value data.frame replacement values for calibPoints slot
##' @rdname calibPoints
setGeneric("calibPoints<-",
           function(object, value) standardGeneric("calibPoints<-"))





##' method changePeakWidth
##' 
##' method changePeakWidth
##' 
##' method changePeakWidth is used to modify the peak width of an individual peak
##' it should be called with the argument mzRange to zoom into the region of
##' interest for choosing the peak. Then two further clicks for choosing the
##' (new) lower and upper peak widths.
##' @param object PeakList object
##' @param selectMz numeric change width of peak closest to selectMz
##' @param lowerWidth numeric lower width value in mass units
##' @param upperWidth numeric upper width value in mass units
##' @param ... additional args
##' @return PeakList object with upated peak widths
##' @export
##' @rdname changePeakWidth
setGeneric("changePeakWidth",
           function(object, selectMz, lowerWidth, upperWidth, ...) standardGeneric("changePeakWidth"))

##' coordToPixel
##' coordToPixel translates xy coordinates from the locator() function
##' to cell coordinates from the image function. Origo is according to
##' ToF-SIMS images the upper left corner. 
##' @param object of class MassImage
##' @param xy numeric vector with x/y locator coordinate
##' @return xy coordinate of MassImage pixels
setGeneric("coordToPixel",
           function(object, xy) standardGeneric("coordToPixel"))

######
### D
######






######
### E
######




######
### F
######






######
### G
######

##' generic method to calculate and get TOFs
##' @export
##' @return vector with ToFs
##' @rdname getTOFs
setGeneric("getTOFs",
           function(object) standardGeneric("getTOFs"))


######
### H
######




######
### I
######

##' @title \code{instrument}, slot of \code{MassSpectra} class objects
##' @description \code{instrument}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @param ... additional args
##' @return content of instrument slot
##' @rdname instrument
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' ## access instrument slot in MassSpectra objects
##' instrument(testSpectra)
##' ## values for the 'instrument' slot can currently be
##' ## 'iontof' or 'ulvacphi'. It is not advisable to 
##' ## change those values manually
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{analysis}} \code{\link{analysisName}} \code{\link{nz}}
##' \code{\link{calibPoints}} \code{\link{calibration}}
setGeneric("instrument", 
           function(object, ...) standardGeneric("instrument"))

##' @export
##' @param value character name of instrument used in the experiment
##' @rdname instrument
setGeneric("instrument<-", 
           function(object,value) standardGeneric("instrument<-"))




##' set a generic method for image
##' @export
##' @param x object object with image data
##' @param ... additional args
##' @return graphical output
##' @rdname image
setGeneric("image",
           function(x, ...) standardGeneric("image"))

##' generic method to obtain imageMatrix
##' @export
##' @param object object of class MassImage
##' @param ... additional args
##' @return numeric matrix
##' @rdname imageMatrix
setGeneric("imageMatrix",
           function(object, ...) standardGeneric("imageMatrix"))

##' generic accessor for iters slot
##' @param object object of class MCR
##' @return content of iters slot
setGeneric("iters",
           function(object) standardGeneric("iters"))





######
### J
######






######
### K
######






######
### L
######






######
### M
######

##' generic for makeTIC
##' @param object object of type MassSpectra
##' @return object of class MassSpectra with TIC
setGeneric("makeTIC",
           function(object) standardGeneric("makeTIC"))



######
### N
######

##' defining generic accessor method for "itzipName"
##' @param object internal
##' @return content of itzipName
setGeneric("itzipName",
           function(object) standardGeneric("itzipName"))

##' generic for setter itzipName
##' @param object internal
##' @param value internal
##' @return object with updated itzipName slot
setGeneric("itzipName<-",
           function(object, value) standardGeneric("itzipName<-"))

##' generic accessor method for slot nComp 
##' @export
##' @param object object of class PCA
##' @return contents of nComp slot
##' @rdname nComp
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage,4)
##' nComp(analysis(testImage,1))
setGeneric("nComp",
           function(object) standardGeneric("nComp"))

##' generic accessor method for slot ndim
##' @param object object of class MassSpectra
##' @return contents of slot ndim
setGeneric("ndim",
           function(object) standardGeneric("ndim"))

##' generic method for 'noPlottingData' aka 'is.null'
##' @param object object of class PCA
##' @return boolean validity check of PCA object
setGeneric("noPlottingData", function(object) standardGeneric("noPlottingData"))

##' generic method for nPeaks
##' @export
##' @param object object of class PeakList
##' @return integer value for number of peaks
##' @rdname nPeaks
setGeneric("nPeaks",
           function(object) standardGeneric("nPeaks"))

##' @title \code{nz}, slot of \code{MassSpectra} class objects
##' @description \code{nz}, slot of \code{MassSpectra} class objects
##' @export
##' @param object object of class MassSpectra
##' @param mzRange vector numeric mass values for nz matrix
##' @return numeric matrix, content of nz
##' @rdname nz
##' @examples
##' library(tofsimsData) 
##' data(tofsimsData)
##' ## access main data slot
##' nz(testSpectra)[,1:1000]
##' @seealso object \code{\link{MassSpectra}} other slots \code{\link{mz}}
##' \code{\link{analysis}} \code{\link{analysisName}} \code{\link{instrument}}
##' \code{\link{calibPoints}} \code{\link{calibration}}
setGeneric("nz",
           function(object, mzRange=NULL) standardGeneric("nz"))

##' @export
##' @param value matrix replacement values for nz
##' @rdname nz
setGeneric("nz<-",
           function(object, value) standardGeneric("nz<-"))






######
### O
######

##' generic overlayPlot
##' 
##' @export
##' @return graphical output
setGeneric("overlayPlot",
           function(objectList, ...) standardGeneric("overlayPlot"))






######
### P
######

##' generic accessor for slot pcaLoadings
##' @export
##' @param object object of class PCA
##' @param comps numeric number of components
##' @return contents of slot pcaLoadings
##' @rdname pcaLoadings
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage,4)
##' plot(pcaLoadings(analysis(testImage,1), comps = c(1,2)))
setGeneric("pcaLoadings",
           function(object, comps=c(1,2)) standardGeneric("pcaLoadings"))

##' generic accessor for slot pcaScores
##' @export
##' @param object object of class PCA
##' @param comps numeric number of components
##' @return contents of slot pcaScores
##' @rdname pcaScores
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage,4)
##' plot(pcaScores(analysis(testImage,1), comps = c(1,2)))
setGeneric("pcaScores",
           function(object, comps=c(1,2)) standardGeneric("pcaScores"))

##' @title \code{peakIDs}, slot of \code{PeakList} class objects
##' @description \code{peakIDs}, slot of \code{PeakList} class objects
##' @export
##' @param object object of class PeakList
##' @return content of slot peakIDs
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra<-calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra<-recalibrate(testSpectra)
##' testSpectra<-unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' peakIDs(testSpectra)[,1:10]
##' @rdname peakIDs
setGeneric("peakIDs",
           function(object) standardGeneric("peakIDs"))

##' @export
##' @param value data.frame
##' @rdname peakIDs
setGeneric("peakIDs<-",
           function(object, value) standardGeneric("peakIDs<-"))

##' @title \code{peakMzs}, slot of \code{PeakList} class objects
##' @description \code{peakMzs}, slot of \code{PeakList} class objects
##' @export
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra<-calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra<-recalibrate(testSpectra)
##' testSpectra<-unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' peakMzs(testSpectra)[,1:10]
##' @param object object of class PeakList
##' @return contents of slot peakMzs
##' @rdname peakMzs
setGeneric("peakMzs",
           function(object) standardGeneric("peakMzs"))


##' @export
##' @param value data.frame
##' @rdname peakMzs
setGeneric("peakMzs<-",
           function(object, value) standardGeneric("peakMzs<-"))

##' generic method peak.pick
##' @export
##' @return object of class PeakList with updated slots PeakIDs and peakMzs
setGeneric("peakPick",
           function(object, span=100, ...) standardGeneric("peakPick"))

##' generic method peaks2Spectra
##' 
##' peaks2Spectra allows to transfer the peaks from a PeakList object onto
##' a MassSpectra object. By this, the MassSpectra object is promoted into
##' a PeakList object
##' @export
##' @param objectPeaks object object of class PeakList
##' @param objectSpectra object object of class MassSpectra
##' @return object of class PeakList
##' @rdname peaks2Spectra 
setGeneric("peaks2Spectra",
           function(objectPeaks,objectSpectra) standardGeneric("peaks2Spectra"))

##' generic method findPeakWidth
##' @export
##' @param object object of class PeakList
##' @param p numeric value for savitzky-golay filter on first derivate
##' @param n numeric value for savitzky-golay filter on first derivate
##' @param span numeric smoothing for determining local minima/maxima values
##' @param widthExtLower numeric factor to extend lower peak width
##' @param widthExtUpper numeric factor to extend upper peak width
##' @param ... additional args
##' @return object of class PeakList with updated peaks
setGeneric("findPeakWidth",
           function(object, p=3, n=5, span=100, widthExtLower=1.5,widthExtUpper=1.75, ...) standardGeneric("findPeakWidth"))

##' Generic method peakWidths
##' 
##' Generic method peakWidths
##' 
##' This method will calculate peak widths (m/z) based on
##' lower and upper widths.
##' @export
##' @param object PeakList object
##' @param plot boolean should there be graphical output
##' @return vector of peak widths
setGeneric("peakWidths",
           function(object, plot=FALSE) standardGeneric("peakWidths"))

##' Generic method for plot
##' @export
##' @return graphical output
setGeneric("plot",
           function(x,y,...) standardGeneric("plot"))

##' generic method points
##' generic method points
##' @export
##' @return graphical output
setGeneric("points",
           function(x, ...) standardGeneric("points"))

##' generic method for "poissonScaling"
##' @export
##' @return object of class MassSpectra with poission 
##' scaled mass spectra in slot nz
setGeneric("poissonScaling",
           function(object, offset=1, ...) standardGeneric("poissonScaling"))






######
### Q
######






######
### R
######

##' Generic method recalibrate
##' @export
##' @return object of class MassSpectra, recalibrated using the data from
##' slots calibPoints
setGeneric("recalibrate",
           function(object) standardGeneric("recalibrate"))

##' generic method reduceSpectrumResolution
##' @export
##' @return object of class MassSpectra with 
##' reduced spectral resolution
setGeneric("reduceSpectrumResolution",
           function(object,everyN=2,mode='remove') standardGeneric("reduceSpectrumResolution"))

##' generic method removePeaks
##' @export
##' @param object object of class PeakList
##' @param mzs M/z's of peaks to be removed
##' @param limit numeric limit for peaks to be removed
##' @param nLocator numeric how many peaks to remove with visual selection
##' @param operator Accept ">", "<", "==", "<=", ">=", "!="
##' @param ... additional args
##' @return object of class PeakList with removed/updated peaks
setGeneric("removePeaks",
           function(object, mzs, operator, limit, nLocator, ...) standardGeneric("removePeaks"))


##' generic accessor method for resids
##' @param object object of class MCR
##' @return content of slot resids
setGeneric("resids",
           function(object) standardGeneric("resids"))


##' generic accessor for RSS
##' @param object object of class MCR
##' @return content of slot RSS 
setGeneric("RSS",
           function(object) standardGeneric("RSS"))





######
### S
######

##' generic method smootherGolay
##' @export
##' @return object of class MassSpectra with updated mass spectra
setGeneric("smootherGolay",
           function(object, p = 3, n = 5, ...) standardGeneric("smootherGolay"))

##' generic smootherSpline
##' @export
##' @param object MassSpectra
##' @param stepsize numeric arg for spline smoother
##' @param spar numeric arg for spline smoother
##' @param ... additional args
##' @return object of class MassSpectra with updated mass spectra
setGeneric("smootherSpline",
           function(object,stepsize=5,spar=0.3,...) standardGeneric("smootherSpline"))

##' generic for smoothScatter
##' @export
##' @param x object of class PCA
##' @param y numeric usually NULL
##' @param nbin numeric
##' @param bandwidth numeric vector length 1 or 2
##' @param colramp numeric
##' @param nrpoints numeric
##' @param ret.selection logical
##' @param pch character
##' @param cex numeric
##' @param col character
##' @param transformation function
##' @param postPlotHook box
##' @param xlab NULL
##' @param ylab NULL
##' @param xlim numeric
##' @param ylim numeric
##' @param xaxs par
##' @param yaxs par
##' @param ... additional args
##' @return graphical output
##' @rdname smoothScatter
setGeneric("smoothScatter")

##' generic for scale
##' @export
##' @return object of class MassSpectra with scaled mass spectra
setGeneric("scale",
           function(x, center = TRUE, scale = TRUE) standardGeneric("scale"))

##' Generic method for subset
##' @export
##' @param x object of class MassImage
##' @param ... additional args
##' @return object of class MassImage a subest of the in-object
##' @rdname subset
setGeneric("subset",
           function(x, ...) standardGeneric("subset"))






######
### T
######






######
### U
######

##' Generic method for unitMassPeaks
##' @export
##' @param object object of class MassSpectra
##' @param mzRange vector numeric with lower and upper mass 
##' range limit for which to set unit mass peaks
##' @param widthAt vector numeric two mass values at which to sample for peak width
##' @param factor vector numeric two values summing up to 1 for 
##' setting assymetric peak width limits
##' @param upper vector numeric upper peak width limits
##' @param lower vectpr numeric lower peak width limits
##' @param ... additional args
##' @return object of class PeakList with unit mass peaks
setGeneric("unitMassPeaks",
           function(object, mzRange, widthAt, factor, upper=NULL, lower=NULL, ...) standardGeneric("unitMassPeaks"))




######
### V
######






######
### W
######






######
### X
######

##' generic accessor method for "xdim"
##' @param object object of class MassImage
##' @return numeric value x dimension of mass image
setGeneric("xdim",
           function(object) standardGeneric("xdim"))

##' generic setter method for "xdim"
##' @param object object of class MassImage
##' @param value numeric x dimension of image
##' @return object of class MassImage with updated x dimension
setGeneric("xdim<-",
           function(object, value) standardGeneric("xdim<-"))

##' @title \code{xy}, slot of \code{MassImage} class objects
##' @description \code{xy}, slot of \code{MassImage} class objects
##' @export
##' @param object object of class MassImage
##' @return vector numeric with xy dimensions of image
##' @rdname xy
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' xy(testImage)
setGeneric("xy",
           function(object) standardGeneric("xy"))

##' @export
##' @param value vector numeric two values for 
##' x and y dimension of image
##' @rdname xy
setGeneric("xy<-",
           function(object, value) standardGeneric("xy<-"))

##' Generic method \code{xySpec}
##'
##' Selection of Spectra
##'
##' @export
##' @param object object of class MassImage
##' @param x numeric x coordinate from where to sample a mass spectra
##' @param y numeric y coordinate from where to sample a mass spectra
##' @return object of class MassSpectra with selected mass spectra
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @rdname xySpec
setGeneric("xySpec",
           function(object, x = NULL, y = NULL) standardGeneric("xySpec"))



######
### Y
######

##' generic accessor method for "ydim"
##' @param object object of class MassImage
##' @return numeric integer, y dimension of image
setGeneric("ydim",
           function(object)
             standardGeneric("ydim"))

##' generic setter method for "ydim"
##' @param object object of class MassImage
##' @param value numeric y dimension of image
##' @return updated object of type MassImage
setGeneric("ydim<-",
           function(object, value) standardGeneric("ydim<-"))



######
### Z
######

##' generic accessor method for "zdim"
##' @param object object of class MassImage
##' @return numeric, number of mass channels / peaks
setGeneric("zdim",
           function(object) standardGeneric("zdim"))





