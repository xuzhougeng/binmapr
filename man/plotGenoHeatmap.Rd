% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotGenoHeatmap}
\alias{plotGenoHeatmap}
\title{genotype heatmap}
\usage{
plotGenoHeatmap(
  obj,
  samples = NULL,
  markers = NULL,
  chrom = NULL,
  start = NULL,
  end = NULL,
  order = NULL,
  show.name = FALSE,
  cols = c("red", "green", "blue")
)
}
\arguments{
\item{obj}{binmapr object}

\item{samples}{samples}

\item{markers}{markers}

\item{chrom}{chromosome/contig/scaffold name, vector}

\item{start}{start position, vector}

\item{end}{end position, vector}

\item{order}{chromosome order, vector equal length to chrom}

\item{show.name}{show sample name in heatmap}
}
\description{
genotype heatmap
}
