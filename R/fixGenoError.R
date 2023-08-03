#' Fix potential genotype error
#'
#' @param x binmapr object
#' @param fix.size define the neighbor size to fix the error
#'
#' @return binmapr object
#' @examples
#' x <- fixGenoError(x, fix.size = 15)
#'
#' @export
#' @author Zhougeng Xu, Guangwei Li
fixGenoError <- function(x, fix.size = 10 ){
  
  
  if (!inherits(binm, "binmapr")) {
    stop("binm is not a binmapr object or its subclass")
  }

  geno <- x$geno
  chroms <- x$CHROM
  samples <- x$ind.name
  
  for (i in samples) {
    for (j in unique(chroms)){
      tmp <- geno[ chroms == j , i ]
      res <- fixGenoErrorHelper(tmp, fix.size = fix.size)
      geno[ chroms == j , i ] <- res
    }
  }
  x$geno <- geno
  x
  
}

fixGenoErrorHelper <- function(geno, fix.size){
  
  # fix potential error with rle
  geno_rle <- rle(geno)
  error_id <- which(geno_rle$lengths < fix.size)
  for(i in error_id){
    left_id <- sum(geno_rle$lengths[1:i]) - geno_rle$lengths[i]
    right_id <- sum(geno_rle$lengths[1:i])
    
    if( i == 1 ){
      geno[(left_id+1):right_id] <- geno[right_id+1]
    }else{
      geno[(left_id+1):right_id] <- geno[left_id]
    }
  }
  return(geno)
}
