#' get the Allele depth from vcf
#'
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract.gt
#' @importFrom vcfR getCHROM
#' @importFrom vcfR getPOS
#' @param file vcf file path
#' @param keep a vector store high confidence site, format should be "chr_1"
#' @param chromosome vector, which chromosome to use
#' if it is NULL, all chromosome will be include in analysis
#' @return a binmapr object with AD and related information
#'
#' @examples
#' library(vcfR)
#' data(vcfR_test)
#' orig_dir <- getwd()
#' temp_dir <- tempdir()
#' setwd( temp_dir )
#' write.vcf( vcfR_test, file = "test.vcf.gz" )
#' ad <- getAdFromVcf("test.vcf.gz")
#' ad
#' # return is full of NA, because the origin vcf don't have INFO/AD
#' setwd( orig_dir )
#'
#' @export
#' @author Zhougeng Xu
getAdFromVcf <- function(file, keep = NULL, chromosome = NULL){

  vcf <- read.vcfR(file)

  AD <- extract.gt(vcf, "AD")
  
  CHROM <- getCHROM(vcf)
  POS <- getPOS(vcf)

  # filter AD by chromosome
  if ( ! is.null(chromosome)){
    idx <- CHROM %in% chromosome
    AD <- AD[ idx, ]
    CHROM <- CHROM[ idx,  ]
    POS  <- POS[ idx, ]
  }
  
  #  only use given SNP
  if ( !is.null(keep) ){
    idx <- row.names(AD) %in% keep
    AD <- AD[ idx,  ]
    CHROM <- CHROM[ idx,  ]
    POS  <- POS[ idx, ]
  }
  
  # fix the bug that AD is NA
  AD[which(is.na(AD))] <- "0,0"
  
  AD_COUNT_LIST <-  strsplit(AD, ",", fixed = TRUE, useBytes = TRUE)
  
  AD_COUNT <- as.integer(unlist(AD_COUNT_LIST))
  Ref_COUNT <- AD_COUNT[seq.int(1,length(AD_COUNT),2)]
  Alt_COUNT <- AD_COUNT[seq.int(2,length(AD_COUNT),2)]
  
  ref_mat <- matrix(Ref_COUNT, nrow = nrow(AD))
  row.names(ref_mat) <- row.names(AD)
  colnames(ref_mat) <- colnames(AD)
  alt_mat <- matrix(Alt_COUNT, nrow = nrow(AD))
  row.names(alt_mat) <- row.names(AD)
  colnames(alt_mat) <- colnames(AD)
  
  structure(list(refmat = ref_mat,
                 altmat = alt_mat,
                 n.ind = dim(AD)[2], 
                 n.mar = dim(AD)[1], 
                 CHROM = CHROM,
                 POS = POS,
                 mar.name = row.names(AD),
                 ind.name = colnames(AD),
                 geno = NULL
                 ), 
            class = "binmapr")

}


# Print method for object class 'binmapr'
#' @export
#' @method print binmapr
print.binmapr <- function (x, ...) {
  ## Print a brief summary of the data
  
  cat("  This is an object of class 'binmapr'\n")
  
  cat("    No. individuals:   ", x$n.ind, "\n")
  cat("    No. markers:       ", x$n.mar, "\n")
  cat("    CHROM information: ", ifelse(is.null(x$CHROM), "no", "yes"), "\n")
  cat("    POS information:   ", ifelse(is.null(x$POS), "no", "yes"), "\n")
  
}
## end of file
