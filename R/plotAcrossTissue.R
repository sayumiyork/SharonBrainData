#' plotAcrossTissue
#'
#' Creates a figure plotting the gene expression using filtered TD-data from the Sharon et al. (2019) dataset across the prefrontal cortex and striatum tissues.
#'
#'
#' @param x A string containing a valid Mouse Genome Informatics (MGI) ID that is found in the dataset
#' @return A figure plotting the gene expression of filtered TD-data from the Sharon et al. (2019) dataset across the prefrontal cortex and striatum tissues.
#' @export




plotAcrossTissue <- function( x ) {
  df <- data.frame(counts(brain, normalized=TRUE)[x,], brain_metadata$tissue, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "tissue")
  df$tissue<- fct_relevel( factor( df$tissue), "STR", "PFC")
  ggplot( df ) +
    geom_bar( aes( tissue, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}

