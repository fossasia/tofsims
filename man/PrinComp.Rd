% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClassDefinition.R, R/princomp.R
\docType{class}
\name{PrinComp-class}
\alias{PrinComp-class}
\alias{PrinComp}
\alias{prinComp}
\title{Class PrinComp}
\usage{
prinComp(dataObject, ...)
}
\arguments{
\item{dataObject}{object of class MassSpectra}

\item{...}{additional args}
}
\value{
object of class prinComp
}
\description{
Class \code{PrinComp} is a wrapper for the S3 function princomp

PrinComp is a PCA constructor function
}
\details{
Class \code{PrinComp} is a wrapper for the S3 function princomp

PrinComp constructor function uses call by reference. The new object is put
into the \code{analysis} slot of the dataObject on which PCA was calculated.
}
\section{Slots}{

\describe{
\item{\code{scale}}{vector see description of \code{stats::princomp}}

\item{\code{n.obs}}{numeric see description of \code{stats::princomp}}

\item{\code{call}}{language see description of \code{stats::princomp}}

\item{\code{center}}{center see description of \code{stats::princomp}}

\item{\code{sdev}}{vector see description of \code{stats::princomp}}
}}

\examples{
testImage <- MassImage('dummy')
testImage<-prinComp(testImage)
image(analysis(testImage, 1), comp = 1)
\dontrun{
library(tofsimsData)
data(tofsimsData)
testImage<-prinComp(testImage)
image(analysis(testImage), 1), comp = 1)}
}
\author{
Lorenz Gerber <lorenz.gerber@slu.se>
}
