#' convet vcf to binmap automatically
#'
#' @param vcf vcf file
#' @param keep high confidence snp site
#' @param chromosome vector, which chromosome to use
#' if it is NULL, all chromosome will be include in analysis
#' @param min.depth parameter of callGtFromAd
#' @param low.ratio parameter of callGtFromAd
#' @param high.ratio parameter of callGtFromAd
#' @param miss.ratio maximum missing ratio
#' @param window.size parameter of callWindowGeno
#' @param high.count parameter of callWindowGeno
#' @param low.count parameter of callWindowGeno
#'
#' @rdname vcf2binmap
#' @export
#' @author Zhou-geng xu
vcf2binmap <- function(vcf, outdir = ".",
                       keep = NULL, chromosome = NULL,
                       min.depth = 10, low.ratio = 0.2, high.ratio = 0.8,
                       miss.ratio = 0.05,
                       window.size = 15, low.count = 6, high.count = 24,
                       fix.size = 5,
                       ...
                       ){

  # read vcf and extract AD
  message("Reading VCF")
  info <- getAdFromVcf(vcf, keep = keep, chromosome = chromosome)

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
  geno_list <- bacthCallGeno(x = GT_flt,
                             CHROM = CHROM,
                             outdir = outdir)

  return(geno_list)
}
