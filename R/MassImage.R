##' \code{MassImage} is also the call to the class constructor. It is used for 
##' importing both BIF/BIF6 and raw image data. 
##' 
##' \code{MassImage}  is also the call to the class constructor. It is used for 
##' importing both BIF/BIF6 and raw image data. 
##' 
##' \code{MassImage} is the user class constructor to obtain a MassImage object.
##' Data can be imported from BIF or raw data files (Iontof or Ulvacphi). To 
##' import raw data, a MassSpectra object with a valid PeakList object has to 
##' be provided as argument.
##' @export
##' @rdname MassImage
##' @param select \code{character}, 'ulvacbif', 'iontofbif',
##'   'iontofgrdpeaks', 'ulvacrawpeaks', 'dummy'
##' @param analysisName \code{character}, name of analysis
##' @param PeakListobj \code{PeakList} class object, used as peaklist for 
##' rawdata import
##' @param untilScan integer or NULL to determine number of ToF-SIMS scans to 
##' import
##' @param ... additional args
##' @return object of class MassImage
##' @author Lorenz Gerber <lorenz.gerber@@slu.se>
##' @examples 
##' # creating dummy data
##' testImage<-MassImage('dummy')
##' image(testImage)
##' \dontrun{
##' # import of rawdata
##' # first a PeakList object has to be created
##' library(tofsimsData)
##' data(tofsimsData)
##' testSpectra <- calibPointNew(testSpectra, mz = 15, value = 15.01551)
##' testSpectra <- calibPointNew(testSpectra, mz = 181, value = 181.0228)
##' testSpectra <- recalibrate(testSpectra)
##' testSpectra <- unitMassPeaks(testSpectra, mzRange = c(1,200), widthAt = c(15, 181), 
##' factor = c(0.4, 0.6), lower = c(14.97, 15.05), upper = c(180.84, 181.43))
##' # obtaining the path to the raw data file in 'tofsims' package
##' importFile<-system.file("rawdata", "trift_test_001.RAW", package = "tofsimsData")
##' rawImportedImage <- MassImage('ulvacrawpeaks', importFile, 
##' PeakListobj = testSpectra)
##' image(rawImportedImage)
##' }
MassImage <- function(select = c("ulvacbif", "iontofbif", "iontofgrdpeaks", 
                               "ulvacrawpeaks", "dummy"), analysisName,
                      PeakListobj = c(), untilScan = NULL, ...) {
    
    choice <- match.arg(select)
    
    switch(choice,
           ulvacbif = {
             import <- readBIF(analysisName, "ulvacphi", mode = "image", ...)
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             mz <- import$mz
             xy <- import$xy
           },
           
           iontofbif = {
             import <- readBIF(analysisName, "iontof", mode = "image")
             
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             mz <- import$mz
             xy <- import$xy
           },
           
           iontofgrdpeaks = {
             import <- import.raw(analysisName, 
                                  mode = "imagepeaks", PeakListobj, untilScan, ...)
             
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             mz <- import$mz
             xy <- import$xy
           },
           
           ulvacrawpeaks = {
             import <- cReadRawPhi(analysisName, 
                                   mode = "imagepeaks", PeakListobj, ...)
             
             analysisName <- import$analysisName
             instrument <- import$instrument
             nz <- import$nz
             mz <- import$mz
             xy <- import$xy
           },
           
           dummy = {
             analysisName  <- "dummyAnalysisName"
             instrument    <- "dummyInstrument"
             nz            <- matrix(round(runif(25600,0,9)),256,100)
             calibration   <- data.frame(intercept = 20000, slope = 20000)
             mz            <- seq(1,100,1)
             xy            <- c(32,8)
           }
           
    )
    
    new("MassImage", analysisName = as.character(basename(analysisName)), 
        instrument = as.character(instrument), 
        nz = as.matrix(nz), mz = as.vector(mz), 
        xy = as.vector(xy))
}

##' method dim for MassImage
##' @param x object of class MassImage
##' @return vector numeric
setMethod("dim", signature(x = "MassImage"), 
    function(x) c(x@xy, dim(x@nz)[2]))


##' method definition 'show' on 'MassImage'
##' show has a generic by default
##' @param object object of class MassImage
##' @return data.frame character
setMethod(show, signature(object = "MassImage"), 
    function(object) {
        if (length(object@analysis) > 0) 
            analysis <- paste(unlist(lapply(1:length(object@analysis), 
                function(i) paste(i, class(object@analysis[[i]]), 
                  sep = "-"))), collapse = ",") else analysis <- "none"
        print(data.frame(
            analysisName = analysisName(object), 
            instrument = instrument(object), 
            no.spectra = format(ndim(object), big.mark = "'"), 
            spectra.points = format(zdim(object), big.mark = "'"), 
            `total ion count` = format(sum(nz(object)), big.mark = "'"), 
            x.image.dimension = format(xy(object)[1], big.mark = "'"), 
            y.image.dimention = format(xy(object)[2], big.mark = "'"), 
            analysis = analysis))
    })

##' Method \code{plot()} for \code{MassImage}
##' 
##' Method defining \code{plot()} for the \code{MassImage} class
##' plot has no generic by default
##' 
##' This method will call \code{plot} method of \code{MassSpectra} class.
##' @param x object of type MassImage
##' @param y missing
##' @param ... additional args
##' @param mzRange vector or lenght two, indicating the mz range to be plotted
##' @param normalize should the mass spectra be normalized
##' @return scatter plot with loading or scores
setMethod("plot", signature(x = "MassImage", y = "missing"), 
          function(x, y, ..., mzRange = c(0, 200), normalize = FALSE) {
    callNextMethod(x, y, ..., mzRange = mzRange, normalize = normalize)
})



##' Method for image on MassImage
##'
##' Method to visualize an IMS Mass Image of class \code{MassImage}
##' @param mzSelect vector, which m/z to combine for visualization. if none are 
##' chosen, the TIC is shownhel
##' @return image plot of the ToF SIMS image data
##' @examples
##' testImage<-MassImage('dummy')
##' image(testImage)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' image(testImage)}
##' @rdname image
setMethod(image, signature(x = "MassImage"), 
    function(x, ..., mzSelect = NULL) {
        if (length(mzSelect) > 0) {
            if (length(mzSelect) == 1) {
                image(matrix(nz(x)[, mzSelect], xdim(x), ydim(x)), 
                      xlab = 'X', ylab = 'Y', axes = FALSE,
                      main = paste('Mass Image m/z ',mzSelect),...)
                axis(side = 1, at = c(0,1), labels = c(1,xdim(x)))
                axis(side = 2, at = c(0,1), labels = c(ydim(x),1))
            } else {
                mainCompound<-'Mass Image m/z\'s: \n'
                for(iii in mzSelect){
                  mainCompound<-paste(mainCompound, iii,', ',sep='')
                }
                mainCompound<-substr(mainCompound,1,nchar(mainCompound)-2)
                graphics::image(matrix(apply(nz(x)[, mzSelect], 1, sum),
                                xdim(x), ydim(x)), xlab = 'X',
                                ylab = 'Y', axes = FALSE,
                                main = mainCompound, ...)
                axis(side = 1, at = c(0,1), labels = c(1,xdim(x)))
                axis(side = 2, at = c(0,1), labels = c(ydim(x),1))
            }
        } else {
            graphics::image(matrix(apply(nz(x), 1, sum), xdim(x), 
                                   ydim(x)), xlab = 'X',
                                   ylab = 'Y', axes = FALSE,
                                   main = 'Mass Image TIC', ...)
            axis(side = 1, at = c(0,1), labels = c(1,xdim(x)))
            axis(side = 2, at = c(0,1), labels = c(ydim(x),1)) 
        }
    })



##' method xySpec on MassImage class objects
##'
##' method xySpec extracts the mass spectra of positon x/y and puts
##' them in a MassSpectar class object
##'
##' Selection of mass spectra by vectors of equal length for x and y.
##' @rdname xySpec
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' spectra100100<-xySpec(testImage, 100,100)
##' plot(spectra100100, type = 'l')
setMethod(xySpec, signature(object = "MassImage"), 
    function(object, x=NULL, y=NULL) {
        
      if(is.null(x) || is.null(y)){
            image(object)
            message('Press \'Esc\' when finished')
            coords <- locator(type = "p")
            
            n <- length(coords$x)
            xyMat<-matrix(0, 2, n)
            
            ### determining pixel coordinates
            for(i in 1:n){
                xyMat[,i] <- coordToPixel(object, c(coords$x[i],coords$y[i]))
            }
            
            x <- xyMat[1,]
            y <- xyMat[2,]
        }
      
      
        if (length(x) != length(y)) 
            stop("not matching number of x and y coordinates")
        if (max(x) > xdim(object) || max(y) > 
            ydim(object)) 
            stop("one or several x/y value(s) larger than image dimensions")
        
        # counter clockwise rotation of x/y
        xRotated <- ydim(object) + 1 - y
        yRotated <- (x)
        n <- length(x)
        
        massspec <- nz(object)[(xRotated[1] - 1) * xdim(object) + yRotated[1], ]
        
        if (n >= 2) {
            for (i in 2:n) {
                massspec <- rbind(massspec, 
                                  nz(object)[(xRotated[i] - 1) * xdim(object) + yRotated[i], ]
                                  )
            }
        } else {
            massspec <- matrix(massspec, 1)
        }
        
        rownames(massspec) <- NULL
        
        new("MassSpectra", 
            analysisName = analysisName(object), 
            instrument = instrument(object), 
            nz = massspec,
            mz = mz(object))
    })



##' \code{Subset} method for objects of class \code{MassImage}
##'
##' @param xyUpperLeft vector of length two with x and y for the upper left 
##' subset corner
##' @param xyLowerRight vector of length two with x and y for the lower right 
##' subset corner
##' @return object of class MassImage
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' subsetTestImage<-subset(testImage, xyUpperLeft = c(1,1), xyLowerRight = c(50,50))
##' image(subsetTestImage)
##' @rdname subset
setMethod(subset, signature(x = "MassImage"), 
    function(x, ..., xyUpperLeft = NULL, xyLowerRight = NULL) {
      
      ## choose subset interactive  
      if (length(xyUpperLeft) != 2 || length(xyLowerRight) != 2) {
            image(x)
            coords <- locator(n = 2, type = "p")
            
            ### translate locator coords to pixel coordinates
            xyUpperLeft <- coordToPixel(x,c(coords$x[1], coords$y[1]))
            xyLowerRight <- coordToPixel(x, c(coords$x[2], coords$y[2]))
            
        }
      
        ### rotate pixel coordinates into correct orientation
        pixelRotation90CW<-function(xy, yLength){
            xyRotated<-c(0,0)
            xyRotated[1] <- yLength+1-xy[2]
            xyRotated[2] <- xy[1]
            return(xyRotated)
        }
        
        xyUpperLeft<-pixelRotation90CW(xyUpperLeft,ydim(x))
        xyLowerRight<-pixelRotation90CW(xyLowerRight, ydim(x))
      
        ## determine number of Rows needed in new nz slot
        newMatNRows <- length(
            seq(xyUpperLeft[1]:xyLowerRight[1])) * 
            length(seq(xyUpperLeft[2]:xyLowerRight[2]))
        
        imageDataNew <- matrix(0, newMatNRows, zdim(x))
      
        loopCols <- xyLowerRight[1]:xyUpperLeft[1]
        loopRows <- xyUpperLeft[2]:xyLowerRight[2]
        
        rowCounter <- 1
        
        for (i in loopCols) {
            for (j in loopRows) {
                imageDataNew[rowCounter, ] <- nz(x)[(i - 1) * xdim(x) + j, ]
                rowCounter <- (rowCounter + 1)
            }
        }
        
        imageDataNew <- matrix(as.integer(imageDataNew), , zdim(x))
        colnames(imageDataNew) <- colnames(nz(x))
        newxy <- c(length(loopRows), length(loopCols))
        
        nz(x) <- imageDataNew
        xy(x) <- newxy
        return(x)
        
    })





##' Method imageMatrix for class MassImage
##' @return matrix numeric
##' @rdname imageMatrix
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' ## the TIC matrix can be extracted 
##' dataMatrix <- imageMatrix(testImage)
##' dim(dataMatrix)
##' ## the matrix can be visualized with the 
##' ## normal image() function
##' image(dataMatrix)
setMethod(imageMatrix, signature(object = "MassImage"), 
    function(object) {
        out <- matrix(apply(nz(object), 1, sum), xdim(object), ydim(object))
        return(out)
    })


##' @rdname binning
##' @examples
##' library(BiocParallel)
##' testImage<-MassImage('dummy')
##' par(mfcol=c(1,2), oma=c(0,0,0,0), mar=c(0,0,0,0))
##' image(testImage)
##' ## the following param will cause to run non parallel
##' register(SerialParam(), default=TRUE)
##' testImage <- binning(testImage,binningFactor = 4)
##' image(testImage)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' par(mfcol=c(1,2), oma=c(0,0,0,0), mar=c(0,0,0,0))
##' image(testImage)
##' testImage <- binning(testImage,binningFactor = 4)
##' image(testImage)}
setMethod(binning, signature(object = "MassImage"), 
    function(object, binningFactor = 2) {
            requireNamespace('BiocParallel')
        
        
        applyBin <- function(object, binningFactor, i) {
            apply(array(as.vector(nz(object)[, i]), 
                        c(binningFactor, 
                          xdim(object)/binningFactor, 
                          binningFactor, 
                          ydim(object)/binningFactor)), 
                  c(2, 4), 
                  sum)
        }
        
        if(BiocParallel::bpworkers(BiocParallel::bpparam())==0)
            BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
            
        applyOut <- BiocParallel::bplapply(1:zdim(object),
            function(i) applyBin(object, binningFactor, i))

        
        newnz <- matrix(unlist(applyOut), 
                        ndim(object)/(binningFactor^2), 
                        zdim(object))
        
        xy(object) <- c(xdim(object)/binningFactor, ydim(object)/binningFactor)
        colnames(newnz) <- colnames(nz(object))
        nz(object) <- newnz
        return(object)
    })

##' coordToPixel
##' 
##' coordToPixel
##' 
##' coordToPixel translates xy coordinates from the locator() function
##' to cell coordinates from the image function. Origo is according to
##' ToF-SIMS images the upper left corner. 
##' @param object of class MassImage
##' @param xy numeric vector with x/y locator coordinate
##' @return xy coordinate of MassImage pixels
setMethod(coordToPixel, signature(object='MassImage',
                                  xy='numeric'),
          function(object, xy ){
            
            singleCoord<-function(coord, nPixel){
                # half a unit cell
                halfUnitCell <- 1/(nPixel-1)*0.5
                
                # full length
                lengthAllCells <- 1/(nPixel-1)*nPixel
                
                # shifted coord
                shiftedCoord <- coord + halfUnitCell
                
                # standardize
                standardCoord <- 1 / lengthAllCells * shiftedCoord
                
                # calculate cell
                pixel <- ceiling(standardCoord * nPixel)
                
                return(pixel)
            }
            
            xyPixel<-c(0,0)
            
            xyPixel[1] <- singleCoord(xy[1],xdim(object))
            xyPixel[2] <- ydim(object) + 1 - singleCoord(xy[2],ydim(object))
            
            return(xyPixel)
          })



