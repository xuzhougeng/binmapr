#' Genotyping all sample by chromsome
#'
#' @description The genotype of each sample will be called according to
#' their allele depth, and then the potential error will be fixed.
#'
#'
#' @param x GT matrix
#' @param CHROM chromosome vector
#' @param outdir outdir
#' @param window.size parameter of callWindowGeno
#' @param high.count parameter of callWindowGeno
#' @param low.count parameter of callWindowGeno
#' @param fix.size parameter of fixGenoError
#' @param pos.start position start index,
#' for exmaple, the pos.start of chr1_1234 is 6, chr01_1234 is 7
#' @param pdf.height pdf width
#' @param pdf.width pdf width
#' @export
#'
#' @author Zhou-Geng xu
batchCallGeno <- function(x, CHROM, outdir = ".",
                          window.size = 15,
                          low.count = 6, high.count = 24,
                          fix.size = 5,
                          pos.start = 6,
                          pdf.height = 4,
                          pdf.width = 8){

  GT_flt <- x
  geno_list <- list()

  for ( chr in CHROM){

    # subset the chromosome
    GT_chr <- GT_flt[grepl(chr , row.names(GT_flt)),]

    pdf_path <- file.path(outdir, paste0(chr, ".pdf"))
    csv_path <- file.path(outdir, paste0(chr, ".csv"))

    pdf(pdf_path, height = pdf.height, width = pdf.width)
    par(mfrow = c(1,2))

    result <- lapply(1:ncol(GT_chr), function(x){
      sample_name <- colnames(GT_chr)[x]
      x <- GT_chr[,x]
      geno <- callWindowGeno(x, window.size = window.size,
                             low = low.count, high = high.count)
      geno_fix <- fixGenoError(geno, fix.size = fix.size)

      plotGeno(geno, pos.start = pos.start,
               ylab = "before fix", title = sample_name)
      plotGeno(geno_fix, pos.start = pos.start,
               ylab = "after fix", title = sample_name)
      return(geno_fix)
    })

    dev.off()

    geno_mt <- Reduce(cbind, result)
    colnames(geno_mt) <- colnames(GT_chr)

    geno_list[[chr]] <- geno_mt
    write.csv(geno_mt, csv_path)
  }

  return(geno_list)
}
