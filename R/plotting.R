#' genotype plotting
#'
#' @param geno genotype vector
#' @param pos.start position start index,
#' for exmaple, the pos.start of chr1_1234 is 6
#' @param xlab x lab
#' @param ylab y lab
#' @param title title
#'
#' @rdname plotGeno
#' @export
#' @author Zhou-geng Xu
plotGeno <- function(geno,
                     pos.start = 6,
                     xlab = "position",
                     ylab = NULL,
                     title = NULL){

  pos <- as.numeric(substring(names(geno), pos.start))

  plot(x=pos, y = geno, pch=20 , cex = 0.2,
       xlab = xlab, ylab = ylab, ylim = c(0,2),
       main = title)
}

#' QTL mapping plot
#'
#' @param pos position of each p value
#' @param pvalue p value
#' @param xlab x lab
#' @param ylab y lab
#' @param chr.name chromosome name
#' @param threshold QTL threshold
#'
#' @rdname plotQtl
#' @export
#' @author Zhou-geng Xu
plotQtl <- function(pos, pvalue, chr.name = NULL,
                    threshold = 2,
                    ...){


  # plot according to position
  df <- data.frame(pos = pos,
                   pvalue = pvalue,
                   lod = -log10(pvalue))

  plot(df$pos, df$lod,
       ylab = "LOD",
       xlab = chr.name,
       axes = FALSE, xaxs = "i", yaxs = "i",
       ylim=c(1,8), pch=20 , cex = 0.2, type="l")

  max_pos <- ceiling(max(df$pos) / 1000000) + 1
  axis(1, at=seq(0,max_pos)* 1000000, labels = seq(0,max_pos))
  axis(2, at=seq(1,10))
  lines(x=c(0,max_pos*1000000),y=c(threshold,threshold), col = "red")

  # add start position and end position
  #start_index <- min(which(df$lod > threshold))
  #end_index <-  min(which(rev(df$lod) > threshold))

  #if ( length(start_index) > 0){
  #  text(grconvertX(0.2,"npc"), grconvertY(0.5, "npc"),
  #       labels = paste0("Start: ", df$pos[start_index]))
  #  text(grconvertX(0.2,"npc"), grconvertY(0.6, "npc"),
  #       labels = paste0("End: ", rev(df$pos)[end_index]))
  #}

}

