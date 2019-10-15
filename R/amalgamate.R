#' amalgamate the neighborhood site
#'
#' @param geno genotype vecotor, with only 0,1,2
#' @param window.len window length to merge the genotype
#' @param delim delimeter to split the position, e.g the delim of
#' "chr_01" should be "_"
#' @return a amalgamated genotype vector
#'
#' @examples
#' geno <- rep(c(rep(1, 10),0, rep(1,3),2,rep(2,3)), 2)
#' names(geno) <- c(paste0("chr1_", c(1,10:15, 20, 30:38, 51)),
#'               paste0("chr2_", c(1,10:15, 20, 30:38, 51)))
#' geno
#' amalgamate(geno, 8)
#'
#' @export
#' @author Zhougeng Xu
amalgamate <- function(geno, window.len = 500, delim = "_"){


  # splite names to chrname and position
  name_splited <- unlist(strsplit(names(geno), delim,
                                  fixed = TRUE,
                                  useBytes = TRUE))

  chr <- name_splited[seq(1,length(name_splited),2)]
  pos <- as.numeric(name_splited[seq(2,length(name_splited),2)])



  l <- lapply(unique(chr), function(x){

    # for debug
    # message(paste0(x))
    pos_chr <- pos[chr == x]
    x_chr <- geno[chr == x]

    # get the window matrix of the position
    index_mt <- findWindow(pos_chr, windowLen = window.len)

    new_x <- findRepresent(index_mt, x_chr)
    return(new_x)

  })

  return(unlist(l))

}

# Depreated function, I have rewrited this function with Rcpp
findWindowDepreated <- function(pos, window.len = 500){

  # The vector only have one window
  last_index <- length(pos)
  if (pos[last_index] - pos[1] < window.len){
    index_mt <-  matrix(c(1, last_index), byrow = T, ncol = 2)
    return(index_mt)
  }

  # index matrix stores the start and end of each window
  index <- 1
  start_index <- 1
  end_index <- 0
  index_mt <- matrix(data = 0, nrow = length(pos), ncol = 2)

  # phyiscal position of genotype
  start_pos <- pos[1]
  end_pos <- 0

  # While loop
  ptr <- 1

  while(ptr <= length(pos)){

    end_pos <- pos[ptr]

    if ((end_pos - start_pos) > window.len){

      # If middle position is isolated, the index of start
      # and end is same
      if (end_index < start_index){
        end_index <- start_index
      }

      # fill in the matrix with the index of start and end
      index_mt[index,] <- c(start_index, end_index)

      # increase the matriw row index
      index <- index + 1

      # set the start position
      start_index <- ptr
      start_pos <- end_pos

      # increase the ptr
      ptr <- ptr + 1

    } else{
      end_index <- ptr
      ptr <- ptr + 1

    }
  }

  # For the case that the last position is still in the window
  if (pos[last_index] - pos[start_index] <  window.len){
    index_mt[index,] <- c(start_index, last_index)
  }

  # For the case that the last position is isolated
  if (pos[last_index] - pos[last_index - 1] > window.len){
    index_mt[index,] <- c(last_index, last_index)
  }

  index_mt <- index_mt[rowSums(index_mt) > 0, ,drop = FALSE]
  return(index_mt)

}
