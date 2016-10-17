##' MAF is a Maximum Autocorrelation Factor Analysis
##' 
##' MAF is a Maximum Autocorrelation Factor Analysis
##' 
##' MAF is a Maximum Autocorrelation Factor Analysis. The code is implemented from the
##' publication of 
##' @export
##' @param dataObject object of type MassImage 
##' @param nComp integer number of components
##' @param usePCA boolean use PCA
##' @return object of type MAF
##' @rdname MAF
MAF <- function(dataObject, nComp = 10, usePCA = TRUE) {
    ### (i) do a linear transformation of the
    ### initial dataset (PCA)
    
    if (usePCA) {
        message("calculating initial pca...")
        nzData <- pcaMAF(nz(dataObject))$vec
    } else {
        message("using image data matrix....")
        nzData <- nz(dataObject)
    }
    
    ### (ii) calculate x and y diffs for a 1
    ### unit shift.
    message("calculating difference matrices...")
    covDiff <- covDiffCalc(nzData, dataObject)
    
    ### (iii) obtain the principal components
    message("calculate eigen vectors and values...")
    eigenOut <- eigen(covDiff)
    loadings <- eigenOut$vectors
    eigenValues <- eigenOut$values
    scores <- nzData %*% loadings
    
    # results<-list(loadings=loadings,scores=scores,values=eigenValues,
    # covDiff=covDiff)
    # return(results)
    
    newAnalysisObjectPosition <- (length(dataObject@analysis) + 1)
    message(newAnalysisObjectPosition)
    dataObject@analysis[[newAnalysisObjectPosition]] <- new(
                "MAF", 
                pcaLoadings = loadings, 
                pcaScores = scores[, 1:nComp], 
                nComp = nComp, 
                imageDim = dim(dataObject), 
                classOfData = as.character(class(dataObject))
            )
    return(dataObject)
    
}

##' covDiffCalc calculates a x/y shift covariance matrix of a multispectral
##' image according to Switzer and Green 1984.
##' @param nzData unknown
##' @param dataObject unknown
##' @return shifted cov matrix
covDiffCalc <- function(nzData, dataObject) {
    xDiff <- array(as.vector(nzData), 
                   dim = c(ydim(dataObject), xdim(dataObject), zdim(dataObject))
                   )
    ### permutating that afterwards row-by-row
    ### vector is eaual to x direction in
    ### original data
    xDiff <- aperm(xDiff, c(2, 1, 3))
    ### calculate the Diff on the whole vector
    xDiff <- diff(as.vector(xDiff))
    xDiff[seq(1, length(xDiff), xdim(dataObject))] <- 0
    xDiff <- c(0, xDiff)
    xDiff <- array(xDiff, 
                   dim = c(xdim(dataObject), ydim(dataObject), zdim(dataObject))
                   )
    xDiff <- aperm(xDiff, c(2, 1, 3))
    xDiff <- matrix(xDiff, 
                    xdim(dataObject) * ydim(dataObject), 
                    zdim(dataObject))
    
    yDiff <- diff(as.vector(nzData))
    yDiff[seq(1, length(yDiff), ydim(dataObject))] <- 0
    yDiff <- c(0, yDiff)
    yDiff <- matrix(yDiff, 
                    xdim(dataObject) * ydim(dataObject))
    
    ### calculate covariance matrices for xDiff
    ### and yDiff and pool them
    message("calculating covariance matrices...")
    # covDiffx <- t(xDiff) %*% xDiff covDiffy
    # <- t(yDiff) %*% yDiff
    covDiffx <- cov(xDiff)
    covDiffy <- cov(yDiff)
    covDiff <- covDiffx + covDiffy
    rm(list = c("covDiffx", "covDiffy"))
    
    return(covDiff)
}


##' helper function for MAF calculation
##' @param X matrix numeric, matrix to calculate PCA from
##' @param nComp number of components
##' @return principal component analysis
pcaMAF <- function(X, nComp) {
    if (nrow(as.matrix(X)) < ncol(as.matrix(X))) {
        tX <- t(X)
        eigenResults <- eigen(X %*% tX)
        
        p <- tX %*% p
        p <- p/sqrt(sum(p^2))
        
        
    } else eigenResults <- eigen(t(X) %*% X)
    
    vec <- X %*% eigenResults$vectors
    pcaOut <- list(vec = vec, eigenResults = eigenResults)
    return(pcaOut)
}