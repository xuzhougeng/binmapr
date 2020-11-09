#' call genotype from Allele depth by Frequency
#'
#' @param x a binmapr object
#' @param  min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param  max.depth maximum depth to infer the genotype, if
#' depth larger than it, it will be conside as NA
#' @param low  threshold to infer one parent, encoded as 0
#' @param high threshold to infer another parent, encoded as 2
#'
#' @details The original GT in VCF may be wrong, if the variant is sequencing error.
#' To avoid this kind of error, we used the ALT/(REF + ALT) to infer
#' the genotype. If the ratio <= 0.2, it will be encoded as 0. If
#' the ratio >= 0.8, it will be encoded as 2. Otherwise, it will be
#' encoded as 1.
#'
#' @return binmapr object, containing  genotype slot
#'
#' @examples
#' AD <- matrix(data = c("30,1","1,30","0,0","15,15"), nrow = 2)
#' row.names(AD) <- c("chr1_1","chr1_100")
#' colnames(AD) <- c("A","B")
#' geno <- callGtFromAd(AD)
#'
#' @export
#' @author Zhougeng Xu
callGtFromAdByFreq <- function(x, min.depth = 10, max.depth = 200,
						 low = 0.2, high = 0.8){

  freq_mt <- calcFreqFromAd(x, min.depth = min.depth,
                               max.depth = max.depth)

  x$geno <- ifelse(freq_mt <= low, 0, ifelse(freq_mt >= high, 2 , 1) )

  return(x)
}

#' Calculate the frequcy of Allele depth
#'
#' @param x a binmapr object
#' @param min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param  max.depth maximum depth to infer the genotype, if
#' depth larger than it, it will be conside as NA
#'
#' @return matrix contains. The alt/(alt+ref) will be caculted
#' from the allele depth
#'
#' @examples
#' AD <- matrix(data = c("30,1","1,30","0,0","15,15"), nrow = 2)
#' row.names(AD) <- c("chr1_1","chr1_100")
#' colnames(AD) <- c("A","B")
#' freq <- calcFreqFromAd(AD)
#'
#' @export
#' @author Zhou-geng Xu
calcFreqFromAd <- function(x, min.depth = 10, max.depth = 200){


  Ref_COUNT <- x$refmat
  Alt_COUNT <- x$altmat
  AD_Depth <- Ref_COUNT + Alt_COUNT

  NA_pos <- which(AD_Depth < min.depth | AD_Depth > max.depth)
  Ref_COUNT[NA_pos] <- NA
  Alt_COUNT[NA_pos] <- NA

  AD_freq <- round(Alt_COUNT / (Ref_COUNT + Alt_COUNT), 2)


  return(AD_freq)

}



#' call genotype from Allele depth with updog
#' 
#' @importFrom updog multidog
#' @importFrom updog format_multidog
#' 
#' @param x a binmapr object
#' @param ploidy The ploidy of the species. Assumed to be the same for each individual.
#' @param model What form should the prior (genotype distribution) take?
#' More detail in help document of updog::multidog
#' @param n.cpus 	The number of computing cores to use
#'
#' @details call genotype with updog
#'
#' @return binmapr object, containing  genotype slot
#'
#' @examples
#' AD <- matrix(data = c("30,1","1,30","0,0","15,15"), nrow = 2)
#' row.names(AD) <- c("chr1_1","chr1_100")
#' colnames(AD) <- c("A","B")
#' geno <- callGtFromAd(AD)
#'
#' @export
#' @author Zhougeng Xu
callGtFromAd <- function(x, 
                         ploidy  = 2,
                         model = "norm",
                         n.cpus = 1){
  
  refmat  <- x$refmat
  sizemat <- x$refmat + x$altmat
  
  mout <- multidog(refmat = refmat,
                   sizemat = sizemat,
                   ploidy = ploidy,
                   model = model,
                   nc = n.cpus)
  
  x$geno <- format_multidog(mout)
  
  return(x)
}
