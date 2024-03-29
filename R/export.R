
#' Generate cross object for R/qtl
#'
#' @param x 
#' @param BC.gen Used only for cross type "bcsft"
#' @param F.gen  Used only for cross type "bcsft"
#' @param alleles A vector of two one-letter character strings 
#' (or four, for the four-way cross), to be used as labels for the two alleles.
#' @param parents a vector with the position or name of two parents
#' @return a binmapr object
#'
#' @export
#' @author Zhougeng Xu
export2cross <- function(x, BC.gen=0, F.gen=0, 
                         alleles = c("A","B"),
                         parents = NULL){
  
  # cross type
  type = "f2"
  if (BC.gen != 0 | F.gen != 0){
    type = "bcsft"
  }
  
  if (!is.null(parents) ) {
    if (length(parents) < 2){
      stop("please give the position of two parent")
    }
  }
  
  mt <- x$geno
  
  if ( is.null(parents)){
    tmp_mt <- mt
    n.row <- nrow(tmp_mt)
    n.col <- ncol(tmp_mt)
    
    data <- matrix(data = NA, nrow = n.row, ncol=n.col)
    
    # setting genotyep base on following rule
    # ref: 0 -> 1
    # het: 1 -> 2
    # alt: 2 -> 3
    for (i in seq.int(1, n.row)){
      for (j in seq.int(1, n.col)) {
        
        if (is.na(tmp_mt[i,j])){
          data[i,j] <- NA
        } else if (tmp_mt[i,j] == 1){
          data[i,j] <- 2
        } else if( tmp_mt[i,j] == 0 ){
          data[i,j] <- 1
        } else if( tmp_mt[i,j] == 2){
          data[i,j] <- 3
        } else{
          data[i,j] <- NA
        }
      }
    }
  } else{
    parent1 <- parents[1]
    parent2 <- parents[2]
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

