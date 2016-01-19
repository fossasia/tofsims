##' Check NULL PCA object
##' @param object object of class PCA
##' @return boolean validity check of class PCA object
setMethod(noPlottingData, signature(object = "PCA"), function(object) {
    if (nComp(object) == 0 && 
        is.null(pcaLoadings(object)) && 
        is.null(pcaScores(object))) 
        return(TRUE)
    return(FALSE)
})

##' image for PCA class type loading plots
##' @param comp numeric which component to visualize
##' @rdname image
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage,3)
##' image(analysis(testImage, 1), comp = 1)
setMethod(image, signature(x = "PCA"), function(x, comp, ...) {
    if (is.null(comp)) 
        comp <- 1 else comp <- as.integer(comp)
    if (!validObject(x) || noPlottingData(x)) 
        stop("No data for imaging.")
    if (comp <= 0 || comp > nComp(x)) 
        stop(paste("Invalid number of components. It should be > 0 && <=", 
            nComp(x)))
    if (validObject(x)) {
        requireNamespace('graphics')
        graphics::image(matrix(pcaScores(x, comp), xdim(x), ydim(x)), ...)
    }
})


##' @param comp numeric which component
##' @rdname imageMatrix
setMethod(imageMatrix, signature(object = "PCA"), function(object, comp, ...) {
    out <- matrix(pcaScores(object, comp), xdim(object), ydim(object))
    return(out)
})

##' smoothScatter method for PCA class
##' 
##' smoothScatter method for PCA class
##' @param comps numeric
##' @param pcType character
##' @param label boolean
##' @param labelThreshold numeric
##' @import KernSmooth
##' @rdname smoothScatter
##' @examples
##' library(tofsimsData)
##' data(tofsimsData)
##' testImage<-PCAnalysis(testImage, nComp = 4)
##' smoothScatter(analysis(testImage, 1), comps = c(1,2), 
##' pcType = 'pcaScores', xlab = 'comp 1', ylab = 'comp 2')
setMethod("smoothScatter", 
          signature(x = "PCA"), 
          function(x, 
                   y = NULL, 
                   nbin = 128, 
                   bandwidth, 
                   colramp = colorRampPalette(c("white", blues9)), 
                   nrpoints = 100,
                   ret.selection = FALSE,
                   pch = ".", 
                   cex = 1, 
                   col = "black", 
                   transformation = function(x) x^0.25, 
                   postPlotHook = box, 
                   xlab = NULL, 
                   ylab = NULL, 
                   xlim, 
                   ylim, 
                   xaxs = par("xaxs"), 
                   yaxs = par("yaxs"), 
                   ..., 
                   comps = c(1, 2), 
                   pcType = "pcaScores", 
                   label = FALSE, 
                   labelThreshold = 1) {
    
    # print(dots <- eval(substitute(alist(...))))
    
    ### input validation
    nComp <- length(comps)
    if (length(pcType) != 1 || !pcType %in% c("pcaScores", "pcaLoadings")) 
        stop("Invalid plot type.")
    if (!validObject(x) || noPlottingData(x)) 
        stop("No data for plotting.")
    if (nComp > 2 || nComp <= 0) 
        stop("Invalid number of components. It should be > 0 <= 2.")
    if (max(comps) > nComp(x) || min(comps) <= 0) 
        stop(paste("Invalid component index. It should be > 0 && <=", nComp(x)))
    
    if (pcType == "pcaScores") {
        plotData <- pcaScores(x)[, comps]
    } else if (pcType == "pcaLoadings") {
        plotData <- pcaLoadings(x)[, comps]
    }
    
    requireNamespace('graphics')
    smoothScatter(plotData, 
                  y, 
                  nbin, 
                  bandwidth, 
                  colramp, 
                  nrpoints,
                  ret.selection,
                  pch, 
                  cex, 
                  col, 
                  transformation, 
                  postPlotHook, 
                  xlab, 
                  ylab, 
                  xlim, 
                  ylim, 
                  xaxs, 
                  yaxs, 
                  ...)
    
})



##' @param comps numeric vector of length two denominating the components to be 
##' plotted
##' @param pcType character 'pcaLoadings' or pcaScores'
##' @param label boolean plot label
##' @param labelThreshold numeric threshhold on which values to plot a label
##' @return scatter loading/score plot
##' @rdname plot
setMethod("plot", 
          signature(x = "PCA"), 
          function(x, 
                   ..., 
                   comps = c(1, 2), 
                   pcType = "pcaLoadings", 
                   label = FALSE, 
                   labelThreshold = 1) {
    
    ### input validation
    nComp <- length(comps)
    if (length(pcType) != 1 || !pcType %in% c("pcaScores", "pcaLoadings")) 
        stop("Invalid plot type.")
    if (!validObject(x) || noPlottingData(x)) 
        stop("No data for plotting.")
    if (nComp > 3 || nComp <= 0) 
        stop("Invalid number of components. It should be > 0 <= 3.")
    if (max(comps) > nComp(x) || min(comps) <= 0) 
        stop(paste("Invalid component index. It should be > 0 && <=", nComp(x)))
    
    if (pcType == "pcaScores") {
        main <- "Scores"
        method <- pcaScores
    } else if (pcType == "pcaLoadings") {
        main <- "Loadings"
        method <- pcaLoadings
    }
    
    xlab <- paste("Component", comps[1])
    if (nComp > 2) {
        plotData <- method(x)[, comps[1:2]]
        ylab <- paste("Component", comps[2])
        col <- method(x)[, comps[3]]
        color <- rainbow(length(col))[rank(col)]
        plot(plotData, col = color, xlab = xlab, ylab = ylab, main = main, ...)
        legend.col(rainbow(length(col)), col)
    } else {
        if (nComp == 1) 
            ylab <- "Magnitude" else ylab <- paste("Component", comps[2])
        plotData <- method(x)[, comps]
        plot(plotData, xlab = xlab, ylab = ylab, main = main, ...)
    }
    
    # if(label){ xData <- dfData[,comps[1]] shownData <-
    # which(abs(sqrt(xData^2 + lblYData^2)) >= labelThreshold *
    # sd(sqrt(xData^2 + lblYData^2))) dfData$ID <- '' dfData$ID[shownData]
    # <- shownData plot <- plot + geom_text(data=dfData, aes(label=ID),
    # vjust=-0.5, color='red') }
    
})




# chooseData <- function(){ require(sp) require(grid) coord <-
# data.frame(x=double(), y=double()) x <- 1 for(i in 1:1000){ locator
# <- grid.locator(unit='npc') grid.points(x=unit(locator$x, 'npc'),
# y=unit(locator$y, 'npc'), pch=16) coord[x, ] <- c(locator$x,
# locator$y) x <- x + 1 } collect data points coordinates<-locator(type
# = 'l') print(coordinates) make polygon object
# coordinates$x[length(coordinates$x)+1]<-coordinates$x[1]
# coordinates$y[length(coordinates$y)+1]<-coordinates$y[1]
# new.poly<-Polygon(coordinates, hole = TRUE) new.poly }

# addlinetoplot <- function(dataset, varx, vary) { list(
# geom_line(data=dataset, aes_string(x=varx, y=vary)),
# geom_point(data=dataset, aes_string(x=varx, y=vary)) ) }

# makeHighLightMask <- function(object, highLight){ x <- xdim(object) y
# <- xdim(object) mask <- matrix(0, x, y) mask[highLight] <- 1
# return(mask) }



# polygonSelection<-function (plotData){ require(sp) coordinates <-
# locator(type = 'l') coordinates$x[length(coordinates$x + 1)] <-
# coordinates$x[1] coordinates$y[length(coordinates$y + 1)] <-
# coordinates$y[1] new.poly <- Polygon(coordinates)

# in.out <- point.in.polygon(plotData[, 1], plotData[,2],
# new.poly@coords[, 1], new.poly@coords[, 2]) return(in.out) }


##' legend.col
##' 
##' legend.col is a helper for the plot function of Scoreplots. It 
##' allows to visualize a third component by a color range. legend.col
##' plots the color range as legend on the side of the plot
##' @param col character color
##' @param lev character levels
##' @return graphical output
legend.col <- function(col, lev) {
    opar <- par
    n <- length(col)
    bx <- par("usr")
    box.cx <- c(bx[2] + (bx[2] - bx[1])/1000, 
                bx[2] + (bx[2] - bx[1])/1000 + (bx[2] - bx[1])/50)
    box.cy <- c(bx[3], bx[3])
    box.sy <- (bx[4] - bx[3])/n
    
    xx <- rep(box.cx, each = 2)
    
    par(xpd = TRUE)
    for (i in 1:n) {
        
        yy <- c(box.cy[1] + (box.sy * (i - 1)), 
                box.cy[1] + (box.sy * (i)), 
                box.cy[1] + (box.sy * (i)), 
                box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
        
    }
    par(new = TRUE)
    plot(0, 
         0, 
         type = "n", 
         ylim = c(min(lev), max(lev)), 
         yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
    axis(side = 4, las = 2, tick = FALSE, line = 0.01)
    par <- opar
}
