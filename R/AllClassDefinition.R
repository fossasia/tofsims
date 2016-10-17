################ Class Unions

setClassUnion("Vector", c("vector", "NULL"))
setClassUnion("Matrix", c("matrix", "Vector"))




######################################################
## Basic Class System, sorted according
## class hierarchy
######################################################





###### MassSpectra

##' Validation method function for class MassImage objects
##' @param object object of class MassSpectra
##' @return boolean class validity test
validMassSpectraObject <- function(object) {
    
    if (!is.character(analysisName(object)) || 
        !is.character(instrument(object)) || 
        !is.numeric(nz(object)) || !is.matrix(nz(object)) || 
        !is.vector(mz(object)) || !is.data.frame(calibration(object))) 
        stop("Invalid parameters")
    TRUE
}

##' Class \code{MassSpectra}
##'
##' Class \code{MassSpectra} is the main data container in the tofsims
##' package as it contains the individual mass spectra.
##'
##' Class \code{MassSpectra} is the main data container of the \code{tofsims} 
##' package, containing the individual mass spectra in the slot \code{nz}. 
##' Additional metadata about the analysis can be found in the slots 
##' \code{analysisName} and \code{instrument}. Values for slope and intercept 
##' of the linear mass calibration equation are stored in the slot 
##' \code{calibration}. The M/z values can be found in \code{mz}. 
##' \code{calibration} allows calculating from M/z values back to 
##' times-of-flight.
##' The slot \code{calibPoints} is used to recalibrate the dataset. It contains 
##' a data.frame with the columns \code{mz} and \code{TOF}. The slot 
##' \code{analysis} of type \code{list}, is used as a container for data 
##' analysis objects. Typically, object of the class \code{MassSpectra} are 
##' constructed during data import using the user constructor function with the 
##' same name as the class, \code{MassSpectra}.
##' @export
##' @slot analysisName character vector with the import filename
##' @slot instrument character vector type of instrument used in the experiment
##' @slot calibration data frame for numerics slope and intercept of the mass 
##' calibration
##' @slot calibPoints data frame for time of flight to maass to charge 
##' calibration
##' @slot nz matrix with rows of ion counts and columns as toftimes or mass to 
##' charge ratios
##' @slot mz vector same length as columns in \code{nz} for mass to charge 
##' values
##' @rdname MassSpectra
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
MassSpectra <- setClass("MassSpectra", 
                        slots = c(analysisName = "character", 
                                  instrument = "character", 
                                  calibration = "data.frame", 
                                  calibPoints = "data.frame", 
                                  nz = "matrix", 
                                  mz = "vector", 
                                  analysis = "list"), 
                        validity = validMassSpectraObject)












###### MassImage


##' Validation method functionf for class MassImage objects
##' @param object object of class MassImage
##' @return boolean class validity test
validMassImageObject <- function(object) {
    if (!is.numeric(xy(object))) 
        stop("Invalid paramteres")
    TRUE
}


##' Class MassImage
##'
##' Class MassImage contains the information to shape a number of mass spectra
##' into an image.
##'
##' Class \code{MassImage} inherits from the classes \code{MassAnalysis} and
##' \code{MassSpectra}. It contains the information to shape a number of mass 
##' spectra into an image.
##' @rdname MassImage
##' @export
##' @slot xy vector giving the pixel dimension of the image
##' @inheritParams MassSpectra
MassImage <- setClass("MassImage", 
                      slots = c(xy = "vector"), 
                      contains = c("MassSpectra"), 
                      validity = validMassImageObject)






###### PeakList

##' Validation method function for class PeakList objects
##' @param object object of class PeakList
##' @return boolean class validity test
validPeakListObject <- function(object) {
    #if (!ndim(object) == 1 || !is.numeric(peakIDs(object)) || 
    #    !is.numeric(peakMzs(object))) 
    #    stop("Invalid paramteres")
    TRUE
}

##' Class PeakList
##'
##' Class PeakList is an extension of TIC class that can hold information about 
##' peaks.
##'  
##' Class \code{PeakList} inherits from the classes \code{MassAnalysis}, 
##' \code{MassSpectra}
##' and \code{TIC}. 
##' @export
##' @slot peakIDs matrix integer ID for peaks
##' @slot peakMzs matrix with mass to charge values for lower, middle and upper 
##' peak values
##' @inheritParams MassSpectra
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @examples 
##' # The typical way to obtain a PeakList object is by
##' # applying some peak picking method to a MassSpectra
##' # below an example using the 'unitMassPeaks' method
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra<-calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra<-calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra<-recalibrate(testSpectra)
##' testSpectra<-unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' show(testSpectra)
PeakList <- setClass("PeakList", 
                     slots = c(peakIDs = "Matrix", 
                               peakMzs = "Matrix"), 
                     contains = c("MassSpectra"), 
                     validity = validPeakListObject)








######################################################### 





#################### Analysis Classes





###### PCA (Virtual Class)

##' Validation method function for class PCA objects
##' @param object object of class PCA
##' @return boolean class validity test
validPCAObject <- function(object) {
    if (
        !(
            (nComp(object) == 0 && 
             is.null(pcaLoadings(object)) && 
             is.null(pcaScores(object))) | 
            (nComp(object) != 0 && 
             !is.null(pcaLoadings(object)) && 
             !is.null(pcaScores(object)) && 
             nComp(object) == dim(pcaScores(object))[2])
        )
        ) {
        stop("Invalid object.")
    } else TRUE
}

##' Class PCA
##'
##' Class \code{PCA} is a virtual class for PCA that will be inherited.
##' 
##' Class \code{PCA} is a virtual class for PCA that will be inherited.
##' @export
##' @slot pcaLoadings matrix that holds the loadings of a principal component 
##' like analysis
##' @slot pcaScores matrix that holds the scores of a principal component like 
##' analysis
##' @slot nComp numeric number of components in the principal component like 
##' analysis
##' @slot imageDim vector x and y values of the image dimension
##' @slot classOfData character a more detailed description of the analysis type
PCA <- setClass("PCA", 
                representation(pcaLoadings = "Matrix", 
                          pcaScores = "Matrix", 
                          nComp = "numeric", 
                          imageDim = "Vector", 
                          classOfData = "character"), 
                prototype(nComp = 0,
                          imageDim = c()),
                contains = "VIRTUAL", 
                validity = validPCAObject)






###### PrinComp

##' Class PrinComp
##'
##' Class \code{PrinComp} is a wrapper for the S3 function princomp
##' 
##' Class \code{PrinComp} is a wrapper for the S3 function princomp
##' @export
##' @slot scale vector see description of \code{stats::princomp}
##' @slot n.obs numeric see description of \code{stats::princomp}
##' @slot call language see description of \code{stats::princomp}
##' @slot center center see description of \code{stats::princomp}
##' @slot sdev vector see description of \code{stats::princomp}
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @rdname PrinComp
PrinComp <- setClass("PrinComp", 
                     slots = c(scale = "vector", 
                               n.obs = "numeric", 
                               call = "language", 
                               center = "vector", 
                               sdev = "vector"), 
                     contains = "PCA", 
                     validity = validPCAObject)






###### PrComp

##' Class PrComp
##'
##' Class \code{PrComp} is a wrapper for the S3 function prcomp
##' 
##' Class \code{PrComp} is a wrapper for the S3 function prcomp
##' @export
##' @slot scale logical see description of \code{stats::prcomp}
##' @slot center vector see description of \code{stats::prcomp}
##' @slot sdev vector see description of \code{stats::prcomp}
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @rdname PrComp
PrComp <- setClass("PrComp", 
                   slots = c(scale = "logical", 
                             center = "vector", 
                             sdev = "vector"), 
                   contains = "PCA", 
                   validity = validPCAObject)






###### PCAnalysis

##' Class PCAnalysis
##'
##' Class \code{PCAnalysis} contains methods for simple PCA analyis
##'
##' Class \code{PCAnalysis} contains methods for simple PCA analysis
##' @export
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @rdname PCAnalysis
PCAnalysis <- setClass("PCAnalysis", 
                       contains = "PCA", 
                       validity = validPCAObject)






###### MAF

##' Class MAF
##' 
##' Class \code{MAF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' 
##' Class \code{MAF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' @export
##' @rdname MAF
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' \dontrun{data(tofsimsData)
##' MAF(testImage,5,TRUE)
##' image(analysis(testImage,1),comp = 1)}
MAF <- setClass("MAF", contains = "PCA", validity = validPCAObject)






###### MNF

##' Class MNF
##' 
##' Class \code{MNF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' 
##' Class \code{MNF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' @export
##' @rdname MNF
MNF <- setClass("MNF", contains = "PCA", validity = validPCAObject)






###### nnMNF

##' Class nnMNF
##' 
##' Class \code{nnMNF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' 
##' Class \code{nnMNF} contains methods for Maximum Autocorrelation Factors 
##' analysis
##' @export
##' @rdname nnMNF
nnMNF <- setClass("nnMNF", contains = "PCA", validity = validPCAObject)


###### MCR

##' Class MCR
##'
##' Class \code{MCR} contains methods for 'Multivariate Curve 
##' Resolution by Alternate Least Squares'
##' 
##' Class \code{MCR} contains methods for 'Multivariate Curve
##' Resolution by Alternate Least Squares'
##' @export
##' @slot RSS numeric residual sum of squares
##' @slot resids matrix with residuals
##' @slot iters numeric number of iterations
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @rdname MCR
MCR <- setClass("MCR", 
                slots = c(RSS = "numeric", 
                          resids = "matrix", 
                          iters = "numeric"), 
                contains = "PCA", 
                validity = validPCAObject)
