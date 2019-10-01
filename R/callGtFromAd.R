#' call genotype from AD
#'
#' @param x AD matrix
#' @param  min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param low  threshold to infer one parent, 0
#' @param high threshold to infer another parent, 2
#'
#' @rdname callGtFromAd
#' @export
#' @author Zhou-geng Xu
callGtFromAd <- function(x, min.depth = 10, low = 0.2, high = 0.8){

  AD_COUNT_LIST <-  strsplit(x, ",", fixed = TRUE, useBytes = TRUE)
  AD_COUNT <- as.integer(unlist(AD_COUNT_LIST))
  Ref_COUNT <- AD_COUNT[seq.int(1,length(AD_COUNT),2)]
  Alt_COUNT <- AD_COUNT[seq.int(2,length(AD_COUNT),2)]
  AD_Depth <- Ref_COUNT + Alt_COUNT

  NA_pos <- which(AD_Depth < min.depth)
  Ref_COUNT[NA_pos] <- NA
  Alt_COUNT[NA_pos] <- NA

  AD_freq <- round(Alt_COUNT / (Ref_COUNT + Alt_COUNT), 2)

  GT <- ifelse(AD_freq < low, 0, ifelse(AD_freq > high, 2 , 1) )

  GT_mt <- matrix(GT, nrow = nrow(x))
  row.names(GT_mt) <- row.names(x)
  colnames(GT_mt) <- colnames(x)
  return(GT_mt)
}
