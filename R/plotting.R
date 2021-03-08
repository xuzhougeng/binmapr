#' genotype heatmap
#' 
#' @param samples samples
#' @param markers markers
#' @param chrom chromosome/contig/scaffold name, vector
#' @param start start position, vector
#' @param end end position, vector
#' @param order chromosome order, vector equal length to chrom
#' @param mark.sample mark sample name in heatmap
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap rowAnnotation
#' @importFrom  RColorBrewer brewer.pal
#' 
#' @export
plotGenoHeatmap <- function(obj,
                            samples = NULL,
                            markers = NULL,
                            chrom = NULL,
                            start = NULL,
                            end = NULL,
                            order = NULL,
                            mark.sample = NULL,
                            cols = c("red","green","blue")
                            ){
  
  mat <- obj$geno
  
  # filter matrix by  samples
  if ( ! is.null(samples)){
    if (sum(samples  %in% ad$ind.name)  < length(samples)){
      stop("some sample is not in your data")
    }
    
    mat <- mat[, samples, drop=FALSE]
  }
  
  # filter matrix by position
  marker_pos <- c()
  if ( !is.null(chrom) && is.null(start) && is.null(end) ){
    # filter by chromosome
    for ( i in seq(1, length(chrom))){
      new_pos <-  which(obj$CHROM == chrom[i])
      if (!is.null(order) && order[i] == "-"){
        new_pos <- rev(new_pos)
      }
      marker_pos <- c(marker_pos, new_pos) 
      }   
  } else if ( ! is.null(chrom) && ! is.null(start) && is.null(end)){
    # filter by chromosome and start
    for ( i in seq(1, length(chrom))){
      new_pos <-which(obj$CHROM == chrom[i] & 
                        obj$POS > start[i])
      if (!is.null(order) && order[i] == "-"){
        new_pos <- rev(new_pos)
      }
      marker_pos <- c(marker_pos, new_pos)     
    }
  } else if ( ! is.null(chrom) && is.null(start) && !is.null(end)){
    # filter by chromosome and end
    for ( i in seq(1, length(chrom))){
      new_pos <- which(obj$CHROM == chrom[i] & 
                         obj$POS < end[i])     
      if (!is.null(order) && order[i] == "-"){
        new_pos <- rev(new_pos)
      }
      marker_pos <- c(marker_pos, new_pos) 
    }
  } else if ( ! is.null(chrom) && !is.null(start) && !is.null(end)){
    # filter by chromosome start and end
    for ( i in seq(1, length(chrom))){
      new_pos <- which(obj$CHROM == chrom[i] & 
                         obj$POS > start[i] & 
                         obj$POS < end[i])
      if (!is.null(order) && order[i] == "-"){
        new_pos <- rev(new_pos)
      }
      marker_pos <- c(marker_pos, new_pos)   
    }
  } else if ( is.null(chrom) && is.null(start) && is.null(end)) {
    # use all
    marker_pos <- seq(1, obj$n.mar)
  } else{
    stop("chrom should not be null")
  }
  mat <- mat[marker_pos, , drop=FALSE]
  
  # add annotation
  # set enough color for contig
  # if the contig is too many, repeat the color
  chrom <- unique(obj$CHROM[marker_pos])
  col_brewer <-  RColorBrewer::brewer.pal(n=9,"Set1")
  if ( length(chrom) <= 9 ){
    chrom_cols <- col_brewer[1:length(chrom)]
  } else{
    reps <- length(chrom) %/% 9
    remainder <- length(chrom) %% 9 
    chrom_cols <- c(rep(col_brewer, times = reps) ,
                     col_brewer[1:remainder])
  }
  names(chrom_cols) <- chrom
  rowAnno <- rowAnnotation(
    chrom = obj$CHROM[marker_pos],
    col = list( chrom = chrom_cols)
  )
  
  ht <- Heatmap(mat,
                col = cols,
                name = "genotype",
                left_annotation = rowAnno,
                cluster_rows = FALSE,
                cluster_columns =  FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE)  

}

