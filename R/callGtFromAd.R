#' call genotype from AD
#'
#' @param x AD matrix
#' @param  min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param  max.depth maximum depth to infer the genotype, if
#' depth larger than it, it will be conside as NA
#' @param low  threshold to infer one parent, 0
#' @param high threshold to infer another parent, 2
#'
#' @export
#' @author Zhou-geng Xu
callGtFromAd <- function(x, min.depth = 10, max.depth = 200, 
						 low = 0.2, high = 0.8){

  freq_mt <- calcFreqFromAd(x, min.depth = min.depth,
                               max.depth = max.depth)

  GT_mt <- ifelse(freq_mt < low, 0, ifelse(freq_mt > high, 2 , 1) )

  return(GT_mt)
}

#' Calculate the frequcy of AD
#' 
#' @param x AD matrix
#' @param min.depth minimum depth to infer the genotype, if
#' depth lower than it, it will be conside as NA
#' @param  max.depth maximum depth to infer the genotype, if
#' depth larger than it, it will be conside as NA
#'
#' @export
#' @author Zhou-geng Xu
calcFreqFromAd <- function(x, min.depth = 10, max.depth = 200){

  # fix the bug that AD is NA
  x[which(is.na(x))] <- "0,0"

  AD_COUNT_LIST <-  strsplit(x, ",", fixed = TRUE, useBytes = TRUE)

  AD_COUNT <- as.integer(unlist(AD_COUNT_LIST))
  Ref_COUNT <- AD_COUNT[seq.int(1,length(AD_COUNT),2)]
  Alt_COUNT <- AD_COUNT[seq.int(2,length(AD_COUNT),2)]
  AD_Depth <- Ref_COUNT + Alt_COUNT

  NA_pos <- which(AD_Depth < min.depth | AD_Depth > max.depth)
  Ref_COUNT[NA_pos] <- NA
  Alt_COUNT[NA_pos] <- NA

  AD_freq <- round(Alt_COUNT / (Ref_COUNT + Alt_COUNT), 2)
  freq_mt <- matrix(AD_freq, nrow = nrow(x))
  row.names(freq_mt) <- row.names(x)
  colnames(freq_mt) <- colnames(x)

  return(freq_mt)

}
