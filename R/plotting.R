#' genotype heatmap
#' 
#' @param obj binmapr object
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




#' Recombination fractions heatmap
#' 
#' @param obj binmapr object
#' @param ctg a vector contain contig name
#' @param ort a vector contain contig orientation
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap rowAnnotation HeatmapAnnotation
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' 
#' @export
plotRfHeatmap <- function(obj, ctg, ori ){
  
  corona <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", 
              "#9467bd", "#8c564b", "#e377c2","#7f7f7f", 
              "#bcbd22", "#17becf", "#ad494a", "#e7ba52", 
              "#8ca252", "#756bb1", "#636363", "#aec7e8", 
              "#ff9896", "#98df8a", "#ffbb78", "#c5b0d5", 
              "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d",
              "#9edae5", "#e7969c", "#e7cb94", "#c7e9c0",
              "#de9ed6", "#d9d9d9")
  
  rf <- obj$rf
  num <- length(ctg)
  idx <- c()
  mns <- substr(colnames(rf), 1,11)
  for (i in seq(1, num)){
    idx_i <- which(mns %in% ctg[i])
    if (ori[i] == "-"){
      idx_i <- rev(idx_i)
    }
    idx <- c(idx, idx_i)
    
  }
  
  rf_idx <- rf[idx,idx]
  
  contig_cols <- corona[1:num]
  names(contig_cols) <- ctg
  ha <- HeatmapAnnotation(contig = mns[idx],
                          col = list(contig=contig_cols))
  row_ha <- rowAnnotation(contig = mns[idx],
                          col = list(contig=contig_cols))
  
  ht <- Heatmap(rf_idx,
                name = "rf",
                col = colorRampPalette(viridis::viridis(7))(100),
                top_annotation = ha,
                left_annotation = row_ha,
                cluster_rows = FALSE,
                cluster_columns =  FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE)
  
  
  
}