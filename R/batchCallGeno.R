#' Genotyping all sample by chromsome
#'
#' @description The genotype of each sample will be called according to
#' their allele depth, and then the potential error will be fixed.

#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#' @importFrom utils write.csv
#'
#' @param x GT matrix
#' @param CHROM chromosome vector
#' @param min.marker.nums the minimum marker number on the chromosome
#' @param collapsed merge the neighbor SNPs into single site
#' @param outdir outdir
#' @param window.len window length to merge the genotype
#' @param delim delimeter to split the position, e.g the delim of
#' "chr_01" should be "_"
#' @inheritParams callWindowGeno
#' @inheritParams fixGenoError
#' @param plot.geno whether to plot the genotype or not
#' @param pos.start position start index,
#' for exmaple, the pos.start of chr1_1234 is 6, chr01_1234 is 7
#' @param pdf.height pdf width
#' @param pdf.width pdf width
#'
#' @return A list object including genotype of each chromosome. It will
#' also create PDF and CSV of these genotype.
#'
#' @examples
#' data(geno)
#' GT <- as.matrix(geno[1:100,5:6])
#' row.names(GT) <- paste0(geno$CHR[1:100], "_", geno$POS[1:100])
#' gt <- batchCallGeno(GT, "chr01",pos.start = 7)
#'
#' @export
#'
#' @author Zhougeng xu
batchCallGeno <- function(x, CHROM,
                          min.marker.nums = 10,
                          collapsed = FALSE,
                          outdir = ".",
                          window.len = 1000,
                          delim = "_",
                          window.size = 15,
                          low.count = 6, high.count = 24,
                          fix.size = 5,
                          pos.start = 6,
                          plot.geno = TRUE,
                          pdf.height = 4,
                          pdf.width = 8){

  GT_flt <- x
  geno_list <- list()

  # flag for show the progress
  flag <- 0
  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  for ( chr in CHROM){

    # flag
    flag <- flag + 1
    setTxtProgressBar(pb, round(flag * 100 / length(CHROM), 0))

    # subset the chromosome
    GT_chr <- GT_flt[grepl(chr , row.names(GT_flt)),]

    if (nrow(GT_chr) < min.marker.nums) { next }

    pdf_path <- file.path(outdir, paste0(chr, ".pdf"))
    csv_path <- file.path(outdir, paste0(chr, ".csv"))

    if (plot.geno){
      pdf(pdf_path, height = pdf.height, width = pdf.width)
      par(mfrow = c(1,2))
    }


    result <- lapply(1:ncol(GT_chr), function(x){

      sample_name <- colnames(GT_chr)[x]
      x <- GT_chr[,x]

      # collpased
      if (collapsed){
        x <- amalgamate(x, window.len = window.len, delim = delim)
      }

      # Call Geno by window
      geno <- callWindowGeno(x, window.size = window.size,
                             low = low.count, high = high.count)

      # Fix the potential error
      geno_fix <- fixGenoError(geno, fix.size = fix.size)

      if (plot.geno){
      plotGeno(geno, pos.start = pos.start,
               ylab = "before fix", title = sample_name)
      plotGeno(geno_fix, pos.start = pos.start,
               ylab = "after fix", title = sample_name)
      }
      return(geno_fix)
    })

    if (plot.geno){
      dev.off()
    }

    geno_mt <- Reduce(cbind, result)
    colnames(geno_mt) <- colnames(GT_chr)

    geno_list[[chr]] <- geno_mt
    write.csv(geno_mt, csv_path)
  }

  return(geno_list)
}
