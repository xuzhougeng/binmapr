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
#'
#' @examples
#' data(geno)
#' GT <- as.matrix(geno[1:5,5:10])
#' rownames(GT) <- paste0(geno$CHR[1:5], "_", geno$POS[1:5])
#' df <- export2asmap(GT)
#'
#' @export
#' @author Zhougeng Xu
export2asmap <- function(geno, pop.type = "DH"){
  
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


#' Generate cross object for R/qtl
#'
#' @param x 
#' @param BC.gen Used only for cross type "bcsft"
#' @param F.gen  Used only for cross type "bcsft"
#'
#' @return data.frame
#'
#' @examples
#' data(geno)
#' GT <- as.matrix(geno[1:5,5:10])
#' rownames(GT) <- paste0(geno$CHR[1:5], "_", geno$POS[1:5])
#' df <- export2asmap(GT)
#'
#' @export
#' @author Zhougeng Xu
export2cross <- function(x, BC.gen=0, F.gen=0, 
                         alleles = c("A","B"),
                         parent1 = NULL, parent2 = NULL){
  
  if (is.null(parent1) | is.null((parent2))){
    stop("please provide parent1 and parent 2")
  }
  type = "f2"
  if (BC.gen != 0 | F.gen != 0){
    type = "bcsft"
  }
  
  mt <- x$geno
  
  p1_vec <- mt[, parent1]
  p1_idx <- match(parent1, colnames(mt))
  p2_vec <- mt[, parent2]
  p2_idx <- match(parent2, colnames(mt))
  
  # filtering the row which parents is NA
  # and only keep the  Homozygote allele
  na.pos <- c(which(is.na(p1_vec)), which( is.na(p1_vec) ))
  heter.pos <- c(which(p1_vec == 1),  which(p2_vec==1))
  same.pos <- which(p1_vec == p2_vec)
  combine.pos <- unique( c(na.pos, heter.pos, same.pos) )
  
  # filter parent 
  if ( length(combine.pos) > 0 ){
    p1_vec <- p1_vec[ -combine.pos ]
    p2_vec <- p2_vec[ -combine.pos ]
    # filter genotype matrix
    tmp_mt <- mt[-combine.pos, c(-p1_idx, -p2_idx), drop=FALSE]
  } else{
    tmp_mt <- mt
  }

  n.row <- nrow(tmp_mt)
  n.col <- ncol(tmp_mt)
  
  data <- matrix(data = NA, nrow = n.row, ncol=n.col)
  
  # setting genotyep base on parent
  for (i in seq.int(1, n.row)){
    for (j in seq.int(1, n.col)) {
      
      if (is.na(tmp_mt[i,j])){
        data[i,j] <- NA
      } else if (tmp_mt[i,j] == 1){
        data[i,j] <- 2
      } else if( tmp_mt[i,j] == p1_vec[j] ){
        data[i,j] <- 1
      } else if( tmp_mt[i,j] == p2_vec[j]){
        data[i,j] <- 3
      } else{
        data[i,j] <- NA
      }
    }
  }
  row.names(data) <- row.names(tmp_mt)
  map <- seq(0, nrow(tmp_mt)-1)
  names(map) <- row.names(tmp_mt)
  
  geno <- vector("list", 1)
  names(geno) <- 1
  geno[[1]] <- list(data = t(data), map = map)
  class(geno[[1]]) <- "A"
  
  pheno <- data.frame(id = colnames(tmp_mt))
  
  cross <- list(geno = geno, pheno = pheno)
  class(cross) <- c(type, "cross")
  
  if (type == "bcsft"){
    attr(cross, "scheme") <- c(BC.gen, F.gen)
  }
  attr(cross, "alleles") <- alleles
  cross
  
}