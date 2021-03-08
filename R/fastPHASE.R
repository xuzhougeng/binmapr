#' write fastPHASE format
#' 
#' @param obj a binmapr object
#' @param file output file name
#' 
#' @export
write.fastphase <- function(obj, file=""){
  
  out_file <- vector(mode = "character", length = obj$n.ind * 3 + 2)
  out_file[1] <- obj$n.ind
  out_file[2] <- obj$n.mar
  genos <- obj$geno
  for (i in seq(1, obj$n.ind)) {
    out_file[2+3*(i-1)+1] <- obj$ind.name[i]
    ref_vector <- vector("character", length = obj$n.mar)
    alt_vector <- vector("character", length = obj$n.mar)
    for (j in seq(1, obj$n.mar)){
      
      if ( is.na(genos[j,i])){
        a <- "?"
        b <- "?"
      } else if( genos[j,i] == 0){
        a <- "0"
        b <- "0"
      } else if (genos[j,i] == 2){
        a <- "1"
        b <- "1"
      } else{
        a <- "0"
        b <- "1"
      }
      
      ref_vector[j] <- a
      alt_vector[j] <- b
    }
    
    out_file[2+3*(i-1)+2] <- paste(ref_vector,collapse = "")
    out_file[2+3*(i-1)+3] <- paste(alt_vector,collapse = "")
  }
  writeLines(out_file, con = file)
}

#' read fastPHASE format
#' 
#' @param file fastaPHASE output file name
#' 
#' @export
read.fastphase <- function(file=""){
  
  in_file <- readLines(file)
  
  # get the genotype  vector
  begin_index <- which(in_file == "BEGIN GENOTYPES")
  end_index <- which(in_file == "END GENOTYPES")
  genotypes <- in_file[(begin_index+1):(end_index-1)]
  n.ind <- length(genotypes) / 3
  
  # sample
  samples <- genotypes[seq(1, length(genotypes), 3)]
  ref_vector <- genotypes[seq(2, length(genotypes), 3)]
  alt_vector <- genotypes[seq(3, length(genotypes), 3)]
  
  ref_code <- gsub(" ","",ref_vector)
  alt_code <- gsub(" ","",alt_vector)
  
  n.mar <- nchar(ref_code[1])
  
  # 
  genos <- matrix(nrow = n.mar, ncol = n.ind)
  colnames(genos) <- samples
  for (i in seq(1, n.ind)) {
    
    for (j in seq(1, n.mar)){
      a <- substr(ref_code[i],j,j)
      b <- substr(alt_code[i],j,j)
      if ( a == "0" && b == "0" ){
        genos[j,i] <- 0
      } else if ( a == "1" && b == "1" ){
        genos[j,i] <- 2
      } else{
        genos[j,i] <- 1
      }
    }
    
  }
  genos
}





