

#' getClusterProfilerGenesMouse
#'
#' Reports the gene names of genes found during a clusterProfiler analysis
#'
#' @param x The results of an over representation test - what genes are significantly up/down regulated - and their correlating KEGG categories. Typically the output of runClusterProfilerMouse
#' @param i A KEGG group from the dot plot. Case sensitive.
#' @return Reports the gene names of genes found during a clusterProfiler analysis


getClusterProfilerGenesMouse <- function (x, i) {
  data.frame( x ) %>%
    filter( str_detect( Description, i ) ) %>%
    pull( geneID ) %>%
    strsplit( "/" ) %>%
    unlist() %>%
    bitr( fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db") %>%
    pull( SYMBOL )
}
