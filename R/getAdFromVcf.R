#' get the AD from vcf
#'
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract.gt
#' @importFrom vcfR getCHROM
#' @param file vcf file path
#' @param keep a vector store high confidence site, format should be "chr_1"
#' @param chromosome vector, which chromosome to use
#' if it is NULL, all chromosome will be include in analysis
#' @export
#' @author Zhou-geng Xu
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
