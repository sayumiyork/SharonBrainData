#Load other libraries
library("DESeq2")
library("clusterProfiler")
library("org.Dm.eg.db")
library("tidyverse")
library("org.Mm.eg.db")
library("MarianesMidgutData")

# Generic functions

formatDESeq2Results <- function( x ) {
  df <- as.data.frame(x)
  df <- data.frame(rownames(df), df)
  colnames(df) <- c("GeneID", colnames(df)[-1])
  rownames(df) <- c()
  return(df)
}

# Fly functions

runClusterProfiler <- function (x) {
  ids <- bitr( x$GeneID, "ENSEMBL", "ENTREZID", "org.Dm.eg.db" )
  kegg <- enrichKEGG(ids$ENTREZID, "dme", keyType="ncbi-geneid")
  kegg@result$Description <- sub( " - Drosophila melanogaster \\(fruit fly\\)", "", kegg@result$Description )
  return(kegg)
}

getClusterProfilerGenes <- function (x, i) {
  data.frame( x ) %>%
    filter( str_detect( Description, i ) ) %>%
    pull( geneID ) %>%
    strsplit( "/" ) %>%
    unlist() %>%
    bitr( fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Dm.eg.db") %>%
    pull( SYMBOL )
}

plotAcrossRegions <- function( x ) {
  df <- data.frame(counts(midgut)[x,], midgut_tsv$condition)
  colnames(df) <- c("counts", "region")
  df <- df[10:30,]
  df$region <- fct_relevel( factor( df$region ), "a1", "a2_3", "Cu", "LFCFe", "Fe", "p1", "p2_4" )
  ggplot( df ) +
    geom_bar( aes( region, counts ), stat="identity" ) +    
    ggtitle(x)+
    theme_bw() + 
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}

#Mouse functions
runClusterProfilerMouse <- function (x) {
  ids <- bitr( x$GeneID, "ENSEMBL", "ENTREZID", "org.Mm.eg.db" )
  kegg <- enrichKEGG(ids$ENTREZID, "mmu", keyType="ncbi-geneid")
  kegg@result$Description <- sub( " - Mus musculus \\(house mouse\\)", "", kegg@result$Description )
  return(kegg)
}

getClusterProfilerGenesMouse <- function (x, i) {
  data.frame( x ) %>%
    filter( str_detect( Description, i ) ) %>%
    pull( geneID ) %>%
    strsplit( "/" ) %>%
    unlist() %>%
    bitr( fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db") %>%
    pull( SYMBOL )
}

plotAcrossTissue <- function( x ) {
  df <- data.frame(counts(brain)[x,], brain_metadata$tissue, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "tissue")  
  df$tissue<- fct_relevel( factor( df$tissue), "STR", "PFC")
  ggplot( df ) +
    geom_bar( aes( tissue, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() + 
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}

plotAcrossMice <- function( x ) {
  df <- data.frame(counts(brain)[x,], brain_metadata$bio_sample, row.names = row.names(brain_metadata))
  colnames(df) <- c("counts", "bio_sample")  
  df$bio_sample <- fct_relevel( factor( df$bio_sample), "PFC_TD_1_1", "PFC_TD_1_2", "PFC_TD_1_3","PFC_TD_1_4",
                               "PFC_TD_2_1", "PFC_TD_2_2", "PFC_TD_2_3","PFC_TD_2_4",
                               "PFC_TD_3_1","PFC_TD_3_2","PFC_TD_3_3","PFC_TD_3_4","PFC_TD_3_5","PFC_TD_3_6","PFC_TD_3_7",
                               "Striatum_TD_1_1","Striatum_TD_1_2","Striatum_TD_1_3","Striatum_TD_1_4",
                               "Striatum_TD_2_1","Striatum_TD_2_2","Striatum_TD_2_3","Striatum_TD_2_4",
                               "Striatum_TD_3_1","Striatum_TD_3_2","Striatum_TD_3_3","Striatum_TD_3_4","Striatum_TD_3_6","Striatum_TD_3_7")
  ggplot( df ) +
    geom_bar( aes( bio_sample, counts ), stat="identity" ) +
    ggtitle(x)+
    theme_bw() + 
    theme( axis.text.x = element_text( angle=90, vjust=0.5 ) )
}


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

###########################################
#Create objects
###########################################
mouse_counts<-read.csv("GEMS_counts_TD_filt.csv", row.names = 1)
mouse_metadata<-read.csv("GEMS_metadata_TD.csv", row.names = 1)

mouse_brain <- DESeqDataSetFromMatrix(countData = brain,
                                      colData = brain_metadata,
                                      design = ~ tissue)
mouse_brain_collapse <- collapseReplicates(mouse_brain, mouse_brain$bio_sample)
mouse_brain_collapse <- DESeq(mouse_brain_collapse)
mouse_brain_collapse
saveRDS(mouse_brain_collapse, file="GEMS_counts_TD_filt_TISSUE_BIOSAMPLE_COLLAPSE.RDS", compress="xz")


###########################################
#Test objects
###########################################

brain<-readRDS("brain.RDS")
brain_metadata<-read.csv("brain_metadata.csv", row.names = 1)

###########################################
#DESeq2 and ClusterProfiler
###########################################

temp <- results(brain, contrast = c("tissue", "STR", "PFC") )
R1_vs_R2 <- formatDESeq2Results(temp)
sig_genes <- filter(R1_vs_R2, padj <= 0.05)
dim(sig_genes)
R1_vs_R2_clusters <- runClusterProfilerMouse(sig_genes)

# Look at the dotplot to see the groups of genes
dim(R1_vs_R2_clusters)
dotplot(R1_vs_R2_clusters, showCategory=34, title="YOUR TITLE HERE", font.size=10, label_format = 50)

# Get the gene symbols for a group of genes
#Replace YOUR CATEGORY HERE with a category from your dotplot
getClusterProfilerGenesMouse(R1_vs_R2_clusters, "Endocytosis")


###########################################
#plotAcrossRegionsMouse
###########################################

plotAcrossTissue("ENSMUSG00000000184")
plotAcrossMice("ENSMUSG00000000184")
plotAcrossHostDonor("ENSMUSG00000000184")

