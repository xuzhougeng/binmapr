#' Generate input for mstmap from callWindowGeno result
#'
#' @param geno gentype matrix or list from callWindowGeno
#'
#' @return data.frame
#' @examples
#' data(geno)
#' GT <- as.matrix(geno[1:5,5:10])
#' rownames(GT) <- paste0(geno$CHR[1:5], "_", geno$POS[1:5])
#' df <- buildMstmapInput(GT)
#'
#' @export
#' @author Zhougeng Xu
buildMstmapInput <- function(geno){

  if ( inherits(geno, "matrix")){
    mt <- geno
  } else if( inherits(geno, "list") ) {
    mt <- Reduce(rbind, geno)
  } else{
    stop("input should be matrix or output from batchCallGeno")
  }

  mt[mt == 0] <- "A"
  mt[mt == 2] <- "B"
  mt[mt == 1] <- "X"
  mt[is.na(mt)] <- "U"
  df <- cbind.data.frame(mt, stringsAsFactors = FALSE)
  return(df)

}
