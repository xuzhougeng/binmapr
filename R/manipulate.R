#' subset the binmaprt data
#' 
#' @param x a binmapr object
#' @param by subset data by "sample" or "marker"
#' @method subset binmapr
#' @export
subset.binmapr <- function(x, by, ...){
  
  if (by == "sample"){
    df <- data.frame(
      index    = seq(1,x$n.ind),
      ind.name = x$ind.name
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
      mar.name = x$mar.name
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
  
  x
}