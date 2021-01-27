#' convet vcf to binmap automatically
#'
#' @inheritParams getAdFromVcf
#' @inheritParams batchCallGeno
#' @param  min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param  max.depth maximum depth to infer the genotype, if
#' depth larger than it, it will be conside as NA
#' @param low.ratio  threshold to infer one parent, encoded as 0
#' @param high.ratio threshold to infer another parent, encoded as 2
#' @param miss.ratio only keep site with miss ratio lower than it
#'
#' @return a list
#' @export
#' @author Zhougeng xu
vcf2binmap <- function(file, outdir = ".",
                       keep = NULL, chromosome = NULL,
                       min.depth = 10, max.depth = 100,
                       low.ratio = 0.2, high.ratio = 0.8,
                       miss.ratio = 0.05,
                       window.size = 15, low.count = 6, high.count = 24,
                       fix.size = 5
                       ){

  # read vcf file and extract AD
  message("Reading VCF")
  info <- getAdFromVcf(file, keep = keep, chromosome = chromosome)

  AD <- info$AD
  CHROM <- unique(info$CHROM)

  # call GT from AD
  message("Genotyping From AD")
  GT <- callGtFromAd(AD, min.depth = min.depth,
                     low = low.ratio, high = high.ratio)

  # filter snp site with high missing ratio
  message("Filter marker with high missing ratio")
  miss_ratio <- rowSums(is.na(GT)) / ncol(GT)
  GT_flt <- GT[miss_ratio < miss.ratio, ]
  geno_list <- batchCallGeno(x = GT_flt,
                             CHROM = CHROM,
                             outdir = outdir)

  return(geno_list)
}
