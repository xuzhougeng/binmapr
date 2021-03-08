#' subset the binmaprt data
#' 
#' @param x a binmapr object
#' @param by subset data by "sample" or "marker"
#' @method subset binmapr
#' @export
subset.binmapr <- function(x, by, ...){
  
  #calculate the missing ratio of individual and marker
  if ( is.na(match("mar.miss", names(x))) ||
       is.na(match("ind.miss", names(x))) ){
    x <- calcMissRatio(x)
  }
  
  if (by == "sample"){
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
    
  } else if ( by == "marker"){
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


#' Calculate the missing ratio of each sample on all site
#'
#' @param x a binmapr object
#'
#' @return x a binmapr object with ind.miss and mar.miss
#'
#' @examples
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