#' Call genotype by fix-size window
#'
#' @param binm, a binmapr object
#' @param window.size default is 15
#' @param low default is 6
#' @param high default is 24
#'
#' @return genotype matrix list
#'
#' @export
#' @author Zhougeng Xu
callWindowGeno <- function(binm, window.size = 15,
                           low = 6, high = 24){
  
  
  if (!inherits(binm, "binmapr")) {
    stop("binm is not a binmapr object or its subclass")
  }
  
  
  chroms <- unique(binm$CHROM)
  
  warn_chrom <- grepl("_", chroms)
  # the gene name with "_" is not allowed in this method,
  if (sum(warn_chrom) > 1){
    warning("there are '_' in  chromosome name, filter for following analysis")
    
    
  }
  
  chroms <- chroms[! grepl("_", chroms)]
  
  geno_list <- list()
  
  # iterate the chromosomes
  for (chrom in chroms){
    
    # subset the marker by chromosome 
    binm_sub <- subset(binm, by = "mar", CHROM==chrom)
    
    # iterate the sample
    samples <- unique(binm_sub$ind.name)
    sample_list <- list()
    for (sample in samples){
      # 构建x
      x <- binm_sub$geno[, sample]
      
      splitSNP <- split(x, ceiling(seq_along(x)/window.size))
      
      batchSNP <- sapply(seq_along(splitSNP), function(x){
        
        snp_num <- length(splitSNP[[x]])
        start_name <- names(splitSNP[[x]][1])
        end_name <- names(splitSNP[[x]][snp_num])
        
        chr_name <- strsplit(start_name,"_",fixed = TRUE)[[1]][1]
        start_pos <- as.numeric(strsplit(start_name,"_",fixed = TRUE)[[1]][2])
        end_pos <- as.numeric(strsplit(end_name,"_",fixed = TRUE)[[1]][2])
        
        total <- sum(splitSNP[[x]], na.rm = TRUE)
        names(total) <- paste(chr_name, floor((start_pos + end_pos) / 2) , sep = "_")
        return(total)
      })
      
      sample_list[[sample]] <- ifelse(batchSNP < low, 0, ifelse(batchSNP > high ,2 , 1))
      
    }
    geno_list[[chrom]] <-  do.call(cbind, sample_list )
    
  }
  
  return(geno_list)
  
  
}


