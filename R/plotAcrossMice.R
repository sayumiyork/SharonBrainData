
#' plotAcrossMice
#'
#' Creates a figure plotting the gene expression using filtered TD-data from the Sharon et al. (2019) dataset across mice.
#'
#'
#' @param x A string containing a valid Mouse Genome Informatics (MGI) ID that is found in the dataset
#' @return A figure plotting the gene expression of filtered TD-data from the Sharon et al. (2019) dataset across mice
#' @export


plotAcrossMice <- function( x ) {
  df <- data.frame(counts(brain, normalized=TRUE)[x,], brain_metadata$bio_sample, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "bio_sample")
  df$bio_sample <- fct_relevel( factor( df$bio_sample), "Striatum_TD_1_1","Striatum_TD_1_2","Striatum_TD_1_3","Striatum_TD_1_4",
                                "Striatum_TD_2_1","Striatum_TD_2_2","Striatum_TD_2_3","Striatum_TD_2_4",
                                "Striatum_TD_3_1","Striatum_TD_3_2","Striatum_TD_3_3","Striatum_TD_3_4","Striatum_TD_3_6","Striatum_TD_3_7",
                                "PFC_TD_1_1", "PFC_TD_1_2", "PFC_TD_1_3","PFC_TD_1_4",
                                "PFC_TD_2_1", "PFC_TD_2_2", "PFC_TD_2_3","PFC_TD_2_4",
                                "PFC_TD_3_1","PFC_TD_3_2","PFC_TD_3_3","PFC_TD_3_4","PFC_TD_3_5","PFC_TD_3_6","PFC_TD_3_7")
  ggplot( df ) +
    geom_bar( aes( bio_sample, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}
