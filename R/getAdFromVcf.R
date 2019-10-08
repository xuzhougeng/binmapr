#' get the AD from vcf
#'
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract.gt
#' @importFrom vcfR getCHROM
#' @param file vcf file path
#' @param keep a vector store high confidence site, format should be "chr_1"
#' @param chromosome vector, which chromosome to use
#' if it is NULL, all chromosome will be include in analysis
#' @return a list with AD and CHROM
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

  # filter AD by chromosome
  if ( ! is.null(chromosome)){

    AD <- AD[CHROM %in% chromosome, ]

  }

  if ( !is.null(keep) ){
    AD <- AD[row.names(AD) %in% keep, ]
  }

  return(list(AD = AD, CHROM = unique(CHROM)))

}
