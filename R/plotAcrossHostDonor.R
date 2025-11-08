
#' plotAcrossHostDonor
#'
#' Creates a figure plotting the gene expression using filtered TD-data from the Sharon et al. (2019) dataset across samples from a given host donor
#'
#'
#' @param x A string containing a valid Mouse Genome Informatics (MGI) ID that is found in the dataset
#' @return A figure plotting the gene expression of filtered TD-data from the Sharon et al. (2019) dataset across samples from a given host donor





plotAcrossHostDonor<- function( x ) {
  df <- data.frame(counts(brain)[x,], brain_metadata$host_donor_tissue, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "host_donor_tissue")
  df$host_donor_tissue <- fct_relevel( factor( df$host_donor_tissue), "TD_1_STR","TD_2_STR","TD_3_STR","TD_1_PFC", "TD_2_PFC", "TD_3_PFC")
  ggplot( df ) +
    geom_bar( aes( host_donor_tissue, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() +
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}
