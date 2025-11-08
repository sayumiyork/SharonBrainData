#' runClusterProfilerMouse
#'
#' Performs a clusterProfiler analysis on mouse data using KEGG categories from org.Mm.eg.db.
#'
#' @param x A data frame containing the results of a DESeq2 analysis
#' @return The results of an over representation test - what genes are significantly up/down regulated - and their correlating KEGG categories
#' @export

runClusterProfilerMouse <- function (x) {
  ids <- bitr( x$GeneID, "ENSEMBL", "ENTREZID", "org.Mm.eg.db" )
  kegg <- enrichKEGG(ids$ENTREZID, "mmu", keyType="ncbi-geneid")
  kegg@result$Description <- sub( " - Mus musculus \\(house mouse\\)", "", kegg@result$Description )
  return(kegg)
}
