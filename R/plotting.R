#' genotype heatmap
#' 
#' @param obj binmapr object
#' @param samples samples
#' @param markers markers
#' @param chrom chromosome/contig/scaffold name, vector
#' @param start start position, vector
#' @param end end position, vector
#' @param order chromosome order, vector equal length to chrom
#' @param show.name show sample name in heatmap
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
                            show.name = FALSE,
                            cols = c("red","green","blue")
                            ){
  
  mat <- obj$geno
  
  # filter matrix by  samples
  if ( ! is.null(samples)){
    if (sum(samples  %in% obj$ind.name)  < length(samples)){
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
 
  storage.mode(mat) <- 'character'

  
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
  
  colors = structure(cols, names = c("0", "1", "2"))
  
  if ( show.name ){
    ht <- Heatmap(mat,
                  col = colors,
                  name = "genotype",
                  left_annotation = rowAnno,
                  cluster_rows = FALSE,
                  cluster_columns =  FALSE,
                  show_row_names = FALSE,
                  show_column_names = TRUE)  
  } else{
    ht <- Heatmap(mat,
                  col = colors,
                  name = "genotype",
                  left_annotation = rowAnno,
                  cluster_rows = FALSE,
                  cluster_columns =  FALSE,
                  show_row_names = FALSE,
                  show_column_names = FALSE)  
  }
  ht

}




#' Recombination fractions heatmap
#' 
#' @param obj binmapr object
#' @param ctg a vector contain contig name
#' @param ort a vector contain contig orientation
#' @param delim delimiter of marker name
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap rowAnnotation HeatmapAnnotation
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis viridis
#' 
#' @export
plotRfHeatmap <- function(obj, ctg, ori,
                          delim="_"){
  
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
  mns <- gsub(paste0(delim,".*"),"",colnames(rf))
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


#' marker distribution in scaffold/contigs
#' @param obj binmapr object
#' @param ctg a vector contain contig name
#' @param ort a vector contain contig orientation
#' @param delim delimiter of marker name
#' @importFrom ggplot2 ggplot scale_fill_gradient ylab
#' @importFrom ggplot2 geom_rect theme scale_x_discrete
#' 
#' @export
plotMarker <- function(obj, min_size =100000, min_num = 10,
                       step_size = 1000000){
  df.list <- split(obj$POS, obj$CHROM)
  df.list.stat <- lapply(df.list, function(x){
    stat = data.frame(max_len = max(x), num = length(x))
  })
  
  stat.df <- do.call(rbind, df.list.stat)
  stat.df$chrom <- rownames(stat.df)
  
  # record the maxmium length
  stat.df <- stat.df[stat.df$max_len > min_size & stat.df$num > min_num,  ]
  stat.df$chrom <- factor(stat.df$chrom, levels = stat.df$chrom)
  
  # create sliding window
  slide_window.list <- list()
  for (i in seq(1,nrow(stat.df))){
    seq_id <- stat.df$chrom[i]
    start <- seq(1,stat.df$max_len[i],step_size)
    end <- c(start[-1], stat.df$max_len[i])
    slide_window <- data.frame(
      chrom = seq_id,
      start = start,
      end = end)
    slide_window.list[[seq_id]] <- slide_window
  }
  
  slide_window_df <- do.call(rbind, slide_window.list)
  slide_window_df$count <- 0
  
  
  for (i in seq(1, length(obj$POS))) {
    seq_id <- obj$CHROM[i]
    pos <- obj$POS[i]
    
    locate <- which(slide_window_df$chrom == seq_id & slide_window_df$start <= pos &
                      slide_window_df$end > pos)
    slide_window_df$count[locate] <- slide_window_df$count[locate] + 1
    
    
  }
  
  slide_window_df$chrom <- factor(slide_window_df$chrom,
                                  levels = unique(slide_window_df$chrom))
  
  # modified by https://www.biostars.org/p/269857/
  ggplot(data = stat.df) + 
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    geom_rect(aes(xmin = as.numeric(chrom) - 0.2, 
                  xmax = as.numeric(chrom) + 0.2, 
                  ymax = max_len, ymin = 0), 
              colour="black", fill = "white") +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) +
    # give the appearance of a discrete axis with chrom labels
    scale_x_discrete(name = "chromosome", limits = as.character(stat.df$chrom)) +
    # add marker count
    geom_rect(data = slide_window_df, aes(xmin = as.numeric(chrom) - 0.2, 
                                          xmax = as.numeric(chrom) + 0.2, 
                                          ymax = end, ymin = start, fill = count))  +
    
    scale_fill_gradient(low='#fffdbf', high = "#a51526") +
    # supress scientific notation on the y-axis
    scale_y_continuous(labels = scales::comma) +
    ylab("region (bp)")
  
}