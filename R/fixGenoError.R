#' Fix potential genotype error
#'
#' @param geno a vector object, storing genotyping information
#' @param fix.size the size of short genotpe error
#'
#' @return vector contains error-fixed genotype
#' @examples
#' genos <- c(1,1,1,1,1,0,1,1,1,1,1,0)
#' fixGenoError(genos, fix.size = 2)
#'
#' @export
#' @author Zhougeng Xu, Guangwei Li
fixGenoError <- function(geno, fix.size = 10 ){

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
