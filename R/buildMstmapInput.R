#' Generate input for mstmap from callWindowGeno result
#'
#' @param geno gentype matrix or list from callWindowGeno
#' @param pop.type Character string specifying the population
#' type of the data frame object. Accepted values are "DH"
#' (doubled haploid), "BC" (backcross), "RILn" (non-advanced
#' RIL population with n generations of selfing) and "ARIL"
#' (advanced RIL) (see Details). Default is "DH".
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
buildMstmapInput <- function(geno, pop.type = "DH"){

  if ( inherits(geno, "matrix")){
    mt <- geno
  } else if( inherits(geno, "list") ) {
    mt <- Reduce(rbind, geno)
  } else{
    stop("input should be matrix or output from batchCallGeno")
  }

  # encode the genotype
  if (pop.type == "DH"){
    mt[mt == 0] <- "A"
    mt[mt == 2] <- "B"
    mt[mt == 1] <- "U"
    mt[is.na(mt)] <- "U"
  } else if (pop.type == "BC"){
    mt[mt == 0] <- "A"
    mt[mt == 2] <- "U"
    mt[mt == 1] <- "B"
    mt[is.na(mt)] <- "U"
  } else if (pop.type == "RILn"){
    mt[mt == 0] <- "A"
    mt[mt == 2] <- "B"
    mt[mt == 1] <- "X"
    mt[is.na(mt)] <- "U"
  } else if(pop.type == "ARIL"){
    mt[mt == 0] <- "A"
    mt[mt == 2] <- "B"
    mt[mt == 1] <- "U"
    mt[is.na(mt)] <- "U"
  } else{
    stop("unknown popultaion type, only DH, BC, RILn, ARIL is allowed")
  }

  df <- cbind.data.frame(mt, stringsAsFactors = FALSE)
  return(df)

}
