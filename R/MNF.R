##' Minimum Noise Fraction Transform  (MNF) using 
##' a \code{nearest neighbour} estimate
##' 
##' This method calculates MNF transform using an 
##' \code{nearest neighbour} estimate as implemented 
##' in \code{mzImage} from Stone et al. (2012). 
##' 
##' Minimum Noise Fraction according Green et al. (1988) but
##' using a nearest neighbour estimate for the noise determination
##' as seen in the package \code{mzImage} from Stone et al. (2012).
##' As the mentioned package is no longer maintained, we used an 
##' archieved version as code base for a new version. The C code was
##' implemented through Rcpp (Eddelbuettel and Francois, 2011). The
##' present function is a user constructur that will create a new
##' analysis slot in the chosen MassSpectra/MassImage object.
##' @export
##' @param dataObject object of type MassImage
##' @param limitSNR numeric
##' @return object of class MNF
##' @rdname nnMNF
##' @examples
##' testImage<-MassImage('dummy')
##' testImage<-MNF(testImage)
##' image(analysis(testImage,1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-nnMNF(testImage)
##' image(analysis(testImage,1), comp = 1)}
nnMNF <- function(dataObject, limitSNR = 1.5) {
    
    vecX <- rep(1:ydim(dataObject), each = xdim(dataObject))
    vecY <- rep(1:ydim(dataObject), xdim(dataObject))
    
    nzData <- nz(dataObject)
    nzData[is.na(nzData)] <- 1
    
    message("Computing Summaries...")
    TIC <- apply(nzData, 2, sum, na.rm = TRUE)
    SD <- apply(nzData, 2, sd, na.rm = TRUE)
    SNR <- apply(nzData, MARGIN = 2, FUN = SNR, x = vecX, y = vecY)
    
    message("MNF transform...")
    noise <- apply(nzData, MARGIN = 2, FUN = computeNoise, x = vecX, y = vecY)
    MNF <- computeMNF(nzData = nzData, 
                        noise = noise, 
                        SNR = SNR, 
                        limitSNR = limitSNR)
    
    ### Compute MNF projections
    V <- MNF$vectors
    V <- V[, rev(1:ncol(V))]
    Z <- nzData[, MNF$ind] %*% V
    Z <- scale(Z)
    
    outMNF <- list(loadings = V, scores = Z)
    
    newAnalysisObjectPosition <- (length(dataObject@analysis) + 1)
    message(newAnalysisObjectPosition)
    dataObject@analysis[[newAnalysisObjectPosition]] <- new(
            "MAF", 
            pcaLoadings = V, 
            pcaScores = Z, 
            nComp = dim(Z)[2], 
            imageDim = dim(dataObject), 
            classOfData = as.character(class(dataObject))
            
    )
    return(dataObject)
    
}

##' Minimum Noise Fraction Transform  (MNF)
##' 
##' This method calculates MNF transform using the 
##' diagonal shift method from Switzer and Green (1984)
##' to estimate the noise. 
##' 
##' Minimum Noise Fraction according Green et al. (1988) using
##' diagonal shift method from Switzer and Green (1984) to estimate 
##' the noise. As the original package \code{mzImage} from 
##' Stone et al. 2012 is no longer maintained, we use it as code 
##' base for the present version. The C code was implemented 
##' through Rcpp (Eddelbuettel and Francois, 2011). Practically,
##' this method uses \code{covDiffCalc} from the MAF method. The
##' present function is a user constructur that will create a new
##' analysis slot in the chosen MassSpectra/MassImage object.
##' @export
##' @param dataObject object of type massImage
##' @return object of class MNF
##' @rdname MNF
##' @export
##' @rdname MNF
##' @examples
##' testImage<-MassImage('dummy')
##' testImage<-MNF(testImage)
##' image(analysis(testImage,1), comp = 1)
##' \dontrun{
##' library(tofsimsData)
##' data(tofsimsData)
##' MNF(testImage)
##' image(analysis(testImage,1), comp = 1)}
MNF <- function(dataObject) {
    
    nzData <- nz(dataObject)
    nzData[is.na(nzData)] <- 1
    
    message("MNF transform...")
    covNoise <- covDiffCalc(nzData = nzData, dataObject = dataObject)
    MNF <- computeMNF(nzData = nzData, covNoise = covNoise)
    
    ### Compute MNF projections
    message("compute MNF projections.")
    V <- MNF$vectors
    V <- V[, rev(1:ncol(V))]
    Z <- nzData[, MNF$ind] %*% V
    Z <- scale(Z)
    
    
    newAnalysisObjectPosition <- (length(dataObject@analysis) + 1)
    message(newAnalysisObjectPosition)
        dataObject@analysis[[newAnalysisObjectPosition]] <- new(
            "MNF", 
            pcaLoadings = V, 
            pcaScores = Z, 
            nComp = dim(Z)[2], 
            imageDim = dim(dataObject), 
            classOfData = as.character(class(dataObject))
            )
        return(dataObject)
    
    
}




##' Signal-to-Noise Ratio (SNR)
##' 
##' SNR function for MNF to calculate Signal to Noise Ratio
##' 
##' function from mzimage to calculate signal-to-noise ratio function
##' @param stat unknown
##' @param x unknown
##' @param y unknown
##' @return matrix numeric with signal-to-noise ratios
SNR <- function(stat, x, y) {
    
    xx <- sort(unique(x))
    yy <- sort(unique(y))
    z <- matrix(NA, nrow = max(xx) - min(xx) + 1, ncol = max(yy) - min(yy) + 1)
    z[cbind(x - min(xx) + 1, y - min(yy) + 1)] <- stat
    
    v1 <- var(c(z), na.rm = TRUE)
    
    signal <- nearestNeighbourMean(z)
    noise <- z - signal
    
    ## caluclate SNR
    v2 <- var(c(noise), na.rm = TRUE)
    ret <- v1/v2
    
    return(ret)
}

##' computeNoise
##' 
##' computeNoise determinates the noise by nearest neighbour estimate.
##' This is a helper function for the nnMNF method.
##'
##' computeNoise determinates the noise by nearest neighbour estimate.
##' This is a helper function for the nnMNF method and originates 
##' from the \code{mzImage} package.
##' @param stat unknown
##' @param x unknown
##' @param y unknown
##' @return matrix numeric noise
computeNoise <- function(stat, x, y) {
    xx <- sort(unique(x))
    yy <- sort(unique(y))
    z <- matrix(NA, nrow = max(xx) - min(xx) + 1, ncol = max(yy) - min(yy) + 1)
    z[cbind(x - min(xx) + 1, y - min(yy) + 1)] <- stat
    mtx <- nearestNeighbourMean(z)
    as.vector((z - mtx)[cbind(x - min(xx) + 1, y - min(yy) + 1)])
}

##' nearestNeighbourMean
##' 
##' nearestNeighbourMean helper for nnMNF 
##'
##' function from mzimage to calculate nearest neighbour means
##' @param x unknown see mzimage
##' @return matrix numeric nearest neighbours
nearestNeighbourMean <- function(x) {
    
    ## Add rows and columns of NAs - in C this makes
    ## dealing with edge effects easier
    x <- cbind(NA, x, NA)
    x <- rbind(NA, x, NA)
    
    nrows <- nrow(x)
    ncols <- ncol(x)
    
    pp <- nnMean(as.double(as.vector(x)), as.integer(nrows), as.integer(ncols))
    res <- matrix(pp, nrows, ncols)
    res[which(is.na(x))] <- NA
    res[which(is.nan(res))] <- NA
    res <- res[-c(1, nrows), ]
    res <- res[, -c(1, ncols)]
    res
}


##' compute MNF
##' 
##' compute MNF, helper for MNF/nnMNF
##' 
##' This is a helper function for the MNF/nnMNF function and originates
##' from the \code{mzImage} package.
##' @param nzData matrix
##' @param noise matrix
##' @param SNR numeric
##' @param ind numeric
##' @param iter boolean
##' @param limitSNR numeric
##' @param covNoise matrix
##' @return MNF transform
computeMNF <- function(nzData = NULL, noise = NULL, 
    SNR = NULL, ind = NULL, iter = TRUE, limitSNR = NULL, 
    covNoise = NULL) {
    
    if (is.null(covNoise)) {
        if (is.null(!limitSNR)) {
            if (is.null(ind)) 
                ind <- which(SNR > limitSNR)
        }
    }
    
    
    if (is.null(ind)) 
        ind <- 1:dim(nzData)[2]
    
    ## data covariance matrix
    A <- cov(nzData[, ind])
    ## noise covariance matrix
    if (is.null(covNoise)) {
        B <- covNoise <- cov(noise[, ind])
    } else {
        B <- covNoise
    }
    
    MNF <- LapackGenEigen(A, B, 1, nrow(A))
    
    ## on error in lapack
    if (MNF$info > nrow(A) & iter) {
        message("Trying once again with", MNF$info - 
            nrow(A) - 1, "instead of", nrow(A), "peaks")
        ## if INFO = N + i, for 1 <= i <= N, then the
        ## leading minor of order i of B is not positive
        ## definite. The factorization of B could not be
        ## completed and no eigenvalues or eigenvectors were
        ## computed.
        MNF$info - nrow(A) - 1  ## Choose this many peaks instead at most
        ind <- order(SNR, decreasing = TRUE)[1:(MNF$info - nrow(A) - 1)]
        MNF <- computeMNF(nzData = nzData, ind = ind, 
            noise = noise, SNR = SNR, limitSNR = limitSNR, 
            iter = FALSE, covNoise = covNoise)
    }
    
    MNF$ind <- ind
    return(MNF)
}


##' LapackGenEigen
##' 
##' LapackGenEigen is helper function for MNF and nnMNF
##' 
##' LapackGenEigen is adapted from the \code{mzImage} package. While
##' it initially used dsygvx from the \code{LAPACK} library, it is now
##' ported to \code{RcppArmadillo}, using the \code{eig_pair} function.
##' @param A matrix
##' @param B matrix 
##' @param IL int start index 
##' @param IU int end index
##' @return list with values, vectors and info 
LapackGenEigen <- function(A, B, IL = 1, IU = 3) {
    W <- double(nrow(A))
    Z <- double((IU - IL + 1) * nrow(A))
    info <- 0
    IL <- IL - 1
    IU <- IU - 1
    pp <- EigenDecompose(A, B, IL, IU)
    values <- rev(as.numeric(pp$eigval[1:(IU - IL + 1)]))
    vectors <- matrix(rev(as.numeric(pp$eigvec)), nrow(A), IU - IL + 1)
    return(list(values = values, vectors = vectors, info = info))
}
