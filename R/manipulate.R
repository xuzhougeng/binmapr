#' subset the binmaprt data
#' 
#' @param x a binmapr object
#' @param by subset data by "ind"(sample) or "mar"(marker)
#' @method subset binmapr
#' @export
subset.binmapr <- function(x, by, ...){
  
  
  #calculate the missing ratio of individual and marker
  if ( is.na(match("mar.miss", names(x))) ||
       is.na(match("ind.miss", names(x))) ){
    x <- calcMissRatio(x)
  }
  
  if (by == "ind"){
    df <- data.frame(
      index    = seq(1,x$n.ind),
      ind.name = x$ind.name,
      ind.miss = x$ind.miss
    )
    df <- subset(df, ...)
    x$refmat <- x$refmat[, df$index, drop=FALSE]
    x$altmat <- x$altmat[, df$index, drop=FALSE ]
    x$n.ind <- length(df$index)
    x$ind.name <- df$ind.name
    if ( !is.null(x$geno)){
      x$geno <- x$geno[, df$index, drop=FALSE]
    }
    
  } else if ( by == "mar"){
    df <- data.frame(
      index    = seq(1,x$n.mar),
      CHROM    = x$CHROM,
      POS      = x$POS,
      mar.name = x$mar.name,
      mar.miss = x$mar.miss
    )
    df <- subset(df, ...)
    x$refmat <- x$refmat[df$index, , drop=FALSE]
    x$altmat <- x$altmat[df$index, , drop=FALSE]
    x$n.mar <- length(df$index)
    x$CHROM <- df$CHROM
    x$POS <- df$POS
    x$mar.name <- df$mar.name
    if ( !is.null(x$geno)){
      x$geno <- x$geno[df$index, , drop=FALSE]
    }
  } else{
    stop("by paramter is illegal")
  }
  x <- calcMissRatio(x)
  x
}



#' Subset method for binmapr objects
#'
#' This function subsets binmapr objects either by individual or by marker.
#'
#' @param x A binmapr object.
#' @param i A numeric vector or logical vector to subset rows.
#' @param j A numeric vector or logical vector to subset columns.
#' @param by A character string to specify the subset method. Can be either "ind" (individual) or "mar" (marker).
#' @param ... Other arguments to be passed to the function.
#' @return A subsetted binmapr object.
#' @export
#' @examples
#' # Create a binmapr object
#' bm <- binmapr(...)
#' # Subset the binmapr object by individual
#' bm1 <- bm[1:5, , by="ind"]
#' # Subset the binmapr object by marker
#' bm2 <- bm[ , 1:5, by="mar"]
`[.binmapr` <- function(x, i, j, by, ...){
  
  #calculate the missing ratio of individual and marker
  if ( is.na(match("mar.miss", names(x))) ||
       is.na(match("ind.miss", names(x))) ){
    x <- calcMissRatio(x)
  }
  
  if (by == "ind"){
    df <- data.frame(
      index    = seq(1,x$n.ind),
      ind.name = x$ind.name,
      ind.miss = x$ind.miss
    )
    df <- df[i,]
    x$refmat <- x$refmat[, df$index, drop=FALSE]
    x$altmat <- x$altmat[, df$index, drop=FALSE ]
    x$n.ind <- length(df$index)
    x$ind.name <- df$ind.name
    if ( !is.null(x$geno)){
      x$geno <- x$geno[, df$index, drop=FALSE]
    }
    
  } else if ( by == "mar"){
    df <- data.frame(
      index    = seq(1,x$n.mar),
      CHROM    = x$CHROM,
      POS      = x$POS,
      mar.name = x$mar.name,
      mar.miss = x$mar.miss
    )
    df <- df[i,]
    x$refmat <- x$refmat[df$index, , drop=FALSE]
    x$altmat <- x$altmat[df$index, , drop=FALSE]
    x$n.mar <- length(df$index)
    x$CHROM <- df$CHROM
    x$POS <- df$POS
    x$mar.name <- df$mar.name
    if ( !is.null(x$geno)){
      x$geno <- x$geno[df$index, , drop=FALSE]
    }
  } else{
    stop("by paramter is illegal")
  }
  x <- calcMissRatio(x)
  x
}



#' Calculate the missing ratio of each sample on all site
#'
#' @param x a binmapr object
#'
#' @return x a binmapr object with ind.miss and mar.miss
#'
#' @export
#' @author Zhou-geng Xu
calcMissRatio <- function(x){
  
  x$ind.miss <- apply(x$geno, 2, function(x){
    sum(is.na(x)) / length(x)
  })
  
  x$mar.miss <- apply(x$geno, 1, function(x){
    sum(is.na(x)) / length(x)
  })
  return(x)
}