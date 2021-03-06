% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenericDefinition.R, R/AllGetters.R,
%   R/AllSetters.R
\docType{methods}
\name{nz}
\alias{nz}
\alias{nz<-}
\alias{nz,MassSpectra,missing-method}
\alias{nz,MassSpectra,numeric-method}
\alias{nz<-,MassSpectra-method}
\title{\code{nz}, slot of \code{MassSpectra} class objects}
\usage{
nz(object, mzRange = NULL)

nz(object) <- value

\S4method{nz}{MassSpectra,missing}(object, mzRange = NULL)

\S4method{nz}{MassSpectra,numeric}(object, mzRange = NULL)

\S4method{nz}{MassSpectra}(object) <- value
}
\arguments{
\item{object}{object of class MassSpectra}

\item{mzRange}{vector numeric mass values for nz matrix}

\item{value}{matrix replacement values for nz}
}
\value{
numeric matrix, content of nz
}
\description{
\code{nz}, slot of \code{MassSpectra} class objects
}
\examples{
library(tofsimsData) 
data(tofsimsData)
## access main data slot
nz(testSpectra)[,1:1000]
}
\seealso{
object \code{\link{MassSpectra}} other slots \code{\link{mz}}
\code{\link{analysis}} \code{\link{analysisName}} \code{\link{instrument}}
\code{\link{calibPoints}} \code{\link{calibration}}
}
