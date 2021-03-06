% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixGenoError.R
\name{fixGenoError}
\alias{fixGenoError}
\title{Fix potential genotype error}
\usage{
fixGenoError(x, fix.size = 10)
}
\arguments{
\item{x}{binmapr object}

\item{fix.size}{define the neighbor size to fix the error}
}
\value{
binmapr object
}
\description{
Fix potential genotype error
}
\examples{
x <- fixGenoError(x, fix.size = 15)

}
\author{
Zhougeng Xu, Guangwei Li
}
