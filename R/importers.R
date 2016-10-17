##' Function to read ToF-SIMS data in the form of preprocessed BIF files
##'
##' This function imports BIF files from IONTOF Surface Lab or ULVAC-PHI's
##' WinCadence. This function reads the data sequential directly from the
##' binary stream. Therefore it's rather slow, but uses less memory
##' than the \code{readBIFParallel} function.
##' @title ToF-SIMS BIF/BIF6 file import
##' @param analysisName : filename of BIF/BIF6 file to read
##' @param instrument : character, 'iontof' or 'ulvacphi'
##' @param mode, 'spectra' or 'image'
##' @return object of type MassImage or MassSpectra
##' @author Lorenz Gerber
readBIF <- function(analysisName, instrument = c('iontof', 'ulvacphi'), mode = c('spectra', 'image')) {
    
    imported <- read_BIF(analysisName, instrument)
    
    message("adapt orientation... ")
    imported$imageData <- imported$imageData[imported$dim1Flip, 
        imported$dim2Flip, ]
    
    message("converting array back to matrix...")
    imported$imageData <- matrix(imported$imageData, 
        imported$nPixels, imported$nIntervals)
    
    #colnames(imported$imageData) <- imported$massValues
    massValues <- imported$massValues
    
    if (mode == "spectra") {
        
        readBIFOut <- list(analysisName = analysisName, 
            nz = imported$imageData, instrument = instrument, 
            mz = imported$massValues)
    }
    
    
    if (mode == "image") {
        
        
        readBIFOut <- list(analysisName = analysisName, 
            instrument = instrument, nz = imported$imageData, 
            mz = imported$massValues, xy = c(imported$nXPixels, 
                imported$nYPixels), preProcess = "bifimport", 
            id = imported$ids, middle = imported$middles, 
            lower = imported$lowers, upper = imported$uppers)
    }
    
    
    
    return(readBIFOut)
}




##' Extracting the data from a Ulvac-phi Tof-SIMS
##' raw header character string.
##'
##' This function takes a raw header character string read by
##' get.raw.header() as input and extracts variable names and values.
##' values are currently forwarded just as character string. This is
##' a helper function for read.raw.phi.
##' @title extract variable names and values from Ulvac-phi
##' ToF-SIMS datafile headers
##' @param header header as a raw character string
##' @return list with two vectors containing variable names
##' and values as characters
##' @author Lorenz Gerber
extract.header.data <- function(header) {
    fields <- strsplit(header, "\r\n")[[1]]
    start.ind <- which(fields == "SOFH") + 
        1
    end.ind <- which(fields == "EOFH") - 
        1
    data <- fields[start.ind:end.ind]
    term <- gsub(": [A-z0-9 .(),+:-]*", "", 
        data)
    value <- gsub("[A-z//()]+: ", "", data)
    header.parsed <- list(term = term, value = value)
    header.parsed
}

##' Function to read raw data from the ulvac-phi trift TOF-SIMS
##'
##' This import function works on data recorded on the ulvac-phi
##' trift ToF-SIMS with WinCadence software version V4.2. Other
##' versions mostl likley will not work. In the current version,
##' data has to be imported with 16bit word length, then converted
##' to 64bit binary and finally converted and read with the
##' word lenghts of the respective variables. Currently, the data
##' is unit mass binned with bins of size one from -0.5 to + 0.5.
##' @title Ulvac phi ToF-SIMS raw data import
##' @param analysisName character
##' @param mode character
##' @param PeakListobj object of class PeakList
##' @param ... additional args
##' @return parsed rawdata for further processing
##' @author Lorenz Gerber, Viet Mai Hoang
cReadRawPhi <- function(analysisName, mode = c('spectra', 'imagepeaks', 'image'), 
    PeakListobj = c(), ...){
    
    ####### Checks
    
    analysisName <- gsub("\\~", Sys.getenv("HOME"), 
        analysisName)
    con <- file(analysisName, "rb")
    message("Importing file ", analysisName, 
        " in ", mode, " mode.")
    
    
    
    ####### Read meta data
    
    message("reading header")
    headerRaw <- readBin(con, character(), size = 1, n = 1)
    headerData <- extract.header.data(headerRaw)
    MassOffset <- as.numeric(headerData$value[grep("MassOffset", 
        headerData$term)])
    MassTime <- as.numeric(headerData$value[grep("Mass/Time", 
        headerData$term)])
    SpecBinSize <- as.numeric(headerData$value[grep("SpecBinSize", 
        headerData$term)])
    SoftwareVersion <- as.character(headerData$value[grep("SoftwareVersion", 
        headerData$term)])
    AcqFilename <- as.character(headerData$value[grep("AcqFilename", 
        headerData$term)])
    AcqPrimGun <- as.character(headerData$value[grep("AcqPrimGun", 
        headerData$term)])
    LmigEmitter <- as.character(headerData$value[grep("LmigEmitter", 
        headerData$term)])
    ImagePixels <- as.numeric(headerData$value[grep("ImagePixels", 
        headerData$term)])
    close(con)
    
    message("Aquisition Filename:", AcqFilename, "")
    message("Software Version:", SoftwareVersion, "")
    message("Primary Ion Gun:", AcqPrimGun, LmigEmitter, "")
    message("Image size:", ImagePixels, "x", ImagePixels, "")
    
    
    
    ####### Mass calibration
    
    message("calculating mass calibration")
    slope <- 1/(MassTime * SpecBinSize * 0.001)
    intercept <- (-slope) * MassOffset
    message("slope:", slope, "")
    message("intercept:", intercept, "")
    
    
    
    ####### Import rawdata
    
    message("calling C++ subroutine to read the binary data...")
    rawdata <- readRawPhiC(analysisName, slope, intercept, ImagePixels)
    
    
    
    
    
    ####### Compose high-resolution spectra
    
    if (mode == "spectra") {
        message("calculate high-resolution mass spectra")
        # highresMzCounts<-table(rawdata[,2])
        # highresMzChannels<-round(as.numeric(names(highresMzCounts)),6)
        # highresMzCounts<-matrix(highresMzCounts,1)
        highresMzCounts <- cTable(rawdata[, 2])
        highresMzChannels <- round(as.numeric(highresMzCounts$vars), 6)
        highresMzCounts <- matrix(highresMzCounts$freqs, 1)
    }
    
    
    ####### Compose peak picked image dataset
    
    if (mode == "imagepeaks") {
        
        ### recalibrate
        message("Recalibrate the mass values with curve from PeakListobj.")
        rawdata[, 2] <- ((sqrt(rawdata[, 
            2]) * slope + intercept - unlist(calibration(PeakListobj)[1]))/
                unlist(calibration(PeakListobj)[2]))^2
        
        peaklist <- 1:nPeaks(PeakListobj)
        
        ### Making Mz lists for extracting peak
        ### start and end
        lowerMzs <- peakMzs(PeakListobj)[1, peaklist]
        upperMzs <- rev(peakMzs(PeakListobj)[3, peaklist])
        
        ### Peak start and end Mz lists need to be
        ### sorted, order needs to be stored for
        ### later re-ordering
        lowerMzsOrder <- c(1:length(peaklist))[order(lowerMzs)]
        lowerMzs <- lowerMzs[order(lowerMzs)]
        upperMzsOrder <- rev(c(1:length(peaklist))[order(upperMzs)])
        upperMzs <- rev(upperMzs[order(upperMzs)])
        
        ### Matrix to be filled with data
        imageData <- matrix(0, ImagePixels * 
            ImagePixels, length(peaklist))
        
        ### sorting the data
        message("Sorting the m/z vector.")
        rawdata <- rawdata[sort.list(rawdata[, 2], 
                                     partial = NULL, 
                                     method = "quick", 
                                     na.last = NA), ]
        
        requireNamespace('BiocParallel')
        
        #### Here need to correct the M/z from the
        #### PeakList object
        
        
        message("Extraction of peak limit indices from rawdata")
        startList <- list(rawVector = rawdata[, 2], 
                          mzs = lowerMzs, 
                          mzsOrder = lowerMzsOrder, 
                          startOrEnd = "start")
        endList <- list(rawVector = rawdata[, 2], 
                        mzs = upperMzs, 
                        mzsOrder = upperMzsOrder, 
                        startOrEnd = "end")
        inputList <- list(startList, endList)
        
        # check if it's possible to do parallel processing otherwise
        # load SerialParam() settings
        if(BiocParallel::bpworkers(BiocParallel::bpparam())==0)
            BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
        
        peakIndices <- BiocParallel::bplapply(inputList, function(i) 
                    do.call(cParIndicesSearch, i))
        
        peakIndices[[2]] <- rev(peakIndices[[2]])
        
        message("Calculation of count data per M/z")
        
        # check if it's possible to do parallel processing otherwise
        # load SerialParam() settings
        if(BiocParallel::bpworkers(BiocParallel::bpparam())==0)
            BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
        
        tableFromData <- BiocParallel::bplapply(seq(along.with = peaklist), 
            function(i) cTable(rawdata[peakIndices[[1]][i]:peakIndices[[2]][i], 1]))
        
        
        message("Compounding contingency tables")
        
        for (i in seq(along.with = peaklist)) {
            message(i, " ", appendLF = FALSE)
            imageData[as.numeric(tableFromData[[i]]$vars), 
                i] <- as.vector(tableFromData[[i]]$freqs)
        }
        
        
        message("reshape image data to standard orientation")
        nmz <- dim(imageData)[2]
        
        imageData <- array(as.vector(imageData), 
            c(ImagePixels, ImagePixels, nmz))
        
        imageData <- imageData[ImagePixels:1, 
            ImagePixels:1, ]
        
        message("transform image data from array to matrix")
        imageData <- apply(imageData, 3, 
            FUN = function(x) as.vector(x))
        
        massValues <- peakMzs(PeakListobj)[2, 
            peaklist]
        
        
        
    }
    
    
    
    ####### return data
    
    message("cleaning up...")
    if (mode == "image") {
        
        #colnames(imageData) <- 1:dim(imageData)[2]
        readRawPhiOut <- list(analysisName = analysisName, 
            instrument = "ulvacphi", nz = imageData, 
            calibration = data.frame(intercept = NULL, 
                slope = NULL), mz = massValues, 
            xy = c(ImagePixels, ImagePixels))
    }
    
    if (mode == "spectra") {
        
        #colnames(highresMzCounts) <- highresMzChannels
        readRawPhiOut <- list(analysisName = analysisName, 
            instrument = "ulvacphi", nz = highresMzCounts, 
            calibration = data.frame(intercept = intercept, 
                slope = slope), mz = highresMzChannels)
    }
    
    if (mode == "imagepeaks") {
        readRawPhiOut <- list(analysisName = analysisName, 
            instrument = "ulvacphi", nz = imageData, 
            calibration = data.frame(intercept = intercept, 
                slope = slope), mz = massValues, 
            xy = c(ImagePixels, ImagePixels))
    }
    
    
    
    
    readRawPhiOut
    
}


##' Function to read raw data.
##'
##' This import function works on GRD and ITZIP format
##' @title Raw data import
##' @param analysisName character
##' @param mode charcter
##' @param PeakListobj object of class PeakList
##' @param untilScan numeric read data up to which scan number
##' @param ... addtional args
##' @return parsed rawdata for further processing
##' @author Lorenz Gerber, Viet Mai Hoang
import.raw <- function(analysisName, 
                       mode = c('spectra', 'imagepeaks'),
                       PeakListobj = c(), untilScan = NULL, ...) {
    
    ####### Checks
    image_dim_prop <- "Registration.Raster.Resolution"
    upper_mass_prop <- "Measurement.UpperMass"
    rawdata <- NULL
    analysisName <- path.expand(analysisName)
    
    if (check.extension(analysisName, "grd")) {
        
        ## Looking for properties file in current
        ## directory
        currentDir <- dirname(file.path(analysisName))
        # library(tools)
        propFile <- paste(currentDir, "/", 
            paste(head(strsplit(basename(file.path(analysisName)), split="\\.")[[1]], 
                       -2),
                  collapse = "."), ".properties.txt", 
            sep = "")
        
        ## Read & get required info from
        ## properties file
        properties_data <- read.csv(propFile, 
            sep = "\t", header = FALSE, skipNul = TRUE, 
            stringsAsFactors = FALSE)
        image_size <- look.for.itzip.property(image_dim_prop[1], 
            properties_data)
        upper_mass <- look.for.itzip.property(upper_mass_prop[1], 
            properties_data)
        
        rawdata <- import(analysisName, "grd", 
            image_size, upper_mass)
        # return(rawdata)
    } else if (check.extension(analysisName, "zip")) {
        
        
        
        ####### Read meta data
        
        message("Unzip the file.")
        temporaryDir <- file.path("temporarydata/")
        unzip(analysisName, exdir = temporaryDir)
        properties_data <- read.csv(paste("temporarydata/", 
            list.files(temporaryDir, pattern = ".properties.txt$"), 
            sep = ""), sep = "\t", header = FALSE, 
            skipNul = TRUE, stringsAsFactors = FALSE)
        image_size <- look.for.itzip.property(image_dim_prop[1], 
            properties_data)
        upper_mass <- look.for.itzip.property(upper_mass_prop[1], 
            properties_data)
        
        ####### Read rawdata
        
        message("Reading the rawdata.")
        rawdata <- import(paste("temporarydata/", 
            list.files("temporarydata/"), 
            sep = ""), "itzip", image_size, 
            upper_mass)
        unlink("temporarydata/", recursive = TRUE)
    }
    
    highestTofs <- as.numeric(rawdata["highestTofs"])
    rawdata <- rawdata["calibratedMatrix"][[1]]
    message("Highest TOF :", highestTofs, 
        "")
    
    #### just have one file to test, and there
    #### is something wrong with the first scan
    #### values, 5 lines, Should be removed
    #### later
    rawdata <- rawdata[-c(1:5), ]
    
    #### select a range of scans
    if (length(untilScan) != 0) {
        message("Removing data above scan ", 
            untilScan, ".")
        rawdata <- rawdata[which(rawdata[, 
            3] <= untilScan), ]
    }
    
    ### remove scan and shot data from rawdata
    rawdata <- rawdata[, c(1, 2)]
    
    ####### Mass calibration
    
    intercept <- 0
    slope <- highestTofs/sqrt(upper_mass)
    message("slope: ", slope, "")
    
    
    
    
    ######## Compose high-resolution spectra
    mode <- match.arg(mode)
    
    switch(mode,
           spectra = {
              #if (mode == "spectra") {
              message("Calculating high-resolution mass spectra.")
              # highresMzCounts<-table(rawdata[,2])
              # highresMzChannels<-as.vector(round(as.numeric(names(highresMzCounts)),6))
              # highresMzCounts<-matrix(highresMzCounts,1)
             
             highresMzCounts <- cTable(rawdata[, 2])
             highresMzChannels <- as.vector(
                round(as.numeric(highresMzCounts$vars), 6))
             highresMzCounts <- matrix(highresMzCounts$freqs, 1)
             
             readRawIontofOut <- list(
                analysisName = analysisName, instrument = "iontof",
                nz = highresMzCounts, calibration = 
                    data.frame(intercept = intercept, slope = slope),
                mz = highresMzChannels)
          },
          
          
          ####### compose peak picked image data set
          imagepeaks = {
              #if (mode == "imagepeaks") {
              
              ### recalibrate
              message("Recalibrate the mass values with curve from PeakListobj.")
              rawdata[, 2] <- ((sqrt(rawdata[, 2]) * 
                                    slope + intercept - 
                                    unlist(calibration(PeakListobj)[1])) /
                                  unlist(calibration(PeakListobj)[2]))^2
              
              peaklist <- 1:nPeaks(PeakListobj)
              
              
              
        
        
              ### Making Mz lists for extracting peak
              ### start and end
              lowerMzs <- peakMzs(PeakListobj)[1, peaklist]
              upperMzs <- rev(peakMzs(PeakListobj)[3, peaklist])
              
              ### Peak start and end Mz lists need to be
              ### sorted, the but the order needs to be
              ### stored for later re-ordering
              lowerMzsOrder <- c(1:length(peaklist))[order(lowerMzs)]
              lowerMzs <- lowerMzs[order(lowerMzs)]
              upperMzsOrder <- rev(c(1:length(peaklist))[order(upperMzs)])
              upperMzs <- rev(upperMzs[order(upperMzs)])
              
              ### Matrix to be filled with data
              imageData <- matrix(0, image_size * image_size, length(peaklist))
              
              ### sorting the data
              message("Sorting the m/z vector.")
              rawdata <- rawdata[sort.list(rawdata[, 2], 
                                           partial = NULL, 
                                           method = "quick", 
                                           na.last = NA), ]
              
              
              
              
              message("Extraction of peak limit indices from rawdata")
              startList <- list(rawVector = rawdata[, 2], 
                                mzs = lowerMzs, 
                                mzsOrder = lowerMzsOrder, 
                                startOrEnd = "start")
              endList <- list(rawVector = rawdata[, 2], 
                              mzs = upperMzs, 
                              mzsOrder = upperMzsOrder, 
                              startOrEnd = "end")
              inputList <- list(startList, endList)
              
              # check if it's possible to do parallel processing otherwise
              # load SerialParam() settings
              if(BiocParallel::bpworkers(BiocParallel::bpparam())==0)
                BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
              
              peakIndices <- BiocParallel::bplapply(inputList, 
                                                    function(i) do.call(cParIndicesSearch, i))
              
              
              
              peakIndices[[2]] <- rev(peakIndices[[2]])
              
              message("Calculation of count data per M/z")
              
              # check if it's possible to do parallel processing otherwise
              # load SerialParam() settings
              if(BiocParallel::bpworkers(BiocParallel::bpparam())==0)
                BiocParallel::register(BiocParallel::SerialParam(), default=TRUE)
              
              tableFromData <- BiocParallel::bplapply(seq(along.with = peaklist), 
                                                      function(i) cTable(rawdata[peakIndices[[1]][i]:peakIndices[[2]][i], 1]))
              
              
              message("Compounding contingency tables")
              
              for (i in seq(along.with = peaklist)) {
                message(i, " ", appendLF = FALSE)
                imageData[as.numeric(tableFromData[[i]]$vars), i] <- 
                  as.vector(tableFromData[[i]]$freqs)
              }
              
              
              
              message("reshape image data to standard orientation")
              nmz <- dim(imageData)[2]
              imageData <- array(as.vector(imageData), c(image_size, image_size, nmz))
              imageData <- aperm(imageData, c(2, 1, 3))
              imageData <- imageData[1:image_size, image_size:1, ]
              
              
              message("transform image data from array to matrix")
              imageData <- apply(imageData, 3, FUN = function(x) as.vector(x))
              
              massValues <- peakMzs(PeakListobj)[2, peaklist]
              
              readRawIontofOut <- list(analysisName = analysisName, 
                                       instrument = "iontof", nz = imageData, 
                                       calibration = data.frame(intercept = intercept, slope = slope), 
                                       mz = massValues, 
                                       xy = c(image_size, image_size))
          }
    )
    message("cleaning up...")
    return(readRawIontofOut)
}

##' Function to check file extension
##'
##' This function is used for check the file extension
##' @title Check file extension
##' @param filepath character
##' @param extension character
##' @return boolean
##' @author Lorenz Gerber, Viet Mai Hoang
check.extension <- function(filepath, extension) {
    if (length(grep(pattern = paste(".", extension, "$", sep = ""), 
                    tolower(filepath))) != 0) 
        return(TRUE)
    return(FALSE)
}

##' Function to extract value by passing property name
##'
##' This function is used to get ITZIP property value by passing its name
##' @title Get ITZIP property value
##' @param itzipName character
##' @param itzipProperties character
##' @return character value from itzipProperties corresponding itzipName
##' @author Lorenz Gerber, Viet Mai Hoang
look.for.itzip.property <- function(itzipName, 
    itzipProperties) {
    factor <- itzipProperties$V3[which(itzipProperties$V1 == itzipName)]
    # return(as.numeric(levels(factor))[factor])
    return(as.numeric(factor))
}

##' helper function for parallel processing in rawdata import routines
##' @param rawVector unknown
##' @param mzs unknown
##' @param mzsOrder unknown
##' @param startOrEnd character 'start' or 'end'
##' @return numeric indicies of time of flight
parIndicesSearch <- function(rawVector, mzs, 
    mzsOrder, startOrEnd = "start") {
    if (startOrEnd == "start") {
        counter <- 1
        largeCounter <- 1
        peakIndices <- rep(0, length(mzs))
        
        
        for (i in rawVector) {
            if (i >= mzs[counter]) {
                peakIndices[mzsOrder[counter]] <- largeCounter
                counter <- counter + 1
                if ((counter - 1) == length(mzs)) {
                    break
                }
            }
            largeCounter <- largeCounter + 1
        }
    } else if (startOrEnd == "end") {
        counter <- 1
        largeCounter <- length(rawVector)
        peakIndices <- rep(0, length(mzs))
        
        for (i in rev(rawVector)) {
            if (i <= mzs[counter]) {
                peakIndices[mzsOrder[counter]] <- largeCounter
                counter <- counter + 1
                if ((counter - 1) == length(mzs)) {
                    break
                }
            }
            largeCounter <- largeCounter - 1
        }
    }
    return(peakIndices)
}
