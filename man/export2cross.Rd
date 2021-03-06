% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{export2cross}
\alias{export2cross}
\title{Generate cross object for R/qtl}
\usage{
export2cross(x, BC.gen = 0, F.gen = 0, alleles = c("A", "B"), parents = NULL)
}
\arguments{
\item{x}{}

\item{BC.gen}{Used only for cross type "bcsft"}

\item{F.gen}{Used only for cross type "bcsft"}

\item{alleles}{A vector of two one-letter character strings 
(or four, for the four-way cross), to be used as labels for the two alleles.}

\item{parents}{a vector with the position or name of two parents}
}
\value{
a binmapr object
}
\description{
Generate cross object for R/qtl
}
\author{
Zhougeng Xu
}
