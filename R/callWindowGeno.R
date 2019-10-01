#' call window genotype
#'
#' @param x a vector object, storing genotype information
#' @param window.size default is 15
#' @param low default is 6
#' @param high default is 24
#'
#' @rdname callWindowGeno
#' @export
#' @author Zhougeng Xu
callWindowGeno <- function(x, window.size = 15,
                           low = 6, high = 24){

  splitSNP <- split(x, ceiling(seq_along(x)/window.size))

  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchSNP <- sapply(seq_along(splitSNP), function(x){

    setTxtProgressBar(pb, round(x * 100 / length(splitSNP), 0))

    snp_num <- length(splitSNP[[x]])
    start_name <- names(splitSNP[[x]][1])
    end_name <- names(splitSNP[[x]][])

    chr_name <- strsplit(start_name,"_",fixed = T)[[1]][1]
    start_pos <- as.numeric(strsplit(start_name,"_",fixed = T)[[1]][2])
    end_pos <- as.numeric(strsplit(end_name,"_",fixed = T)[[1]][2])

    total <- sum(splitSNP[[x]], na.rm = T)
    names(total) <- paste(chr_name, floor((start_pos + end_pos) / 2) , sep = "_")
    return(total)
  })

  geno <- ifelse(batchSNP < low, 0, ifelse(batchSNP > high ,2 , 1))

  return(geno)
}
