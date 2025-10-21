#' Differential expression analysis using Sharon et al. (2020) data.
#'
#'  Expression data for regions of the mouse brain after inoculation of gut microbiota from typically developing humans. Data has been filtered to contain only sequences that were found in at least half of samples with a count of at least 10 reads.
#'
#' @format A DESeq2 object containing results of differential expression analysis:
#' \describe{
#'   \item{assays}{List of matrices with raw and normalized count data}
#'   \item{rowRanges}{Genomic ranges for each gene}
#'   \item{colData}{Data frame with sample metadata}
#'   \item{design}{Formula representing the experimental design}
#'   \item{dispersions}{Estimated dispersion values}
#'   \item{results}{Data frame with differential expression results}
#' }
#' @source \url{doi:10.1016/j.cell.2019.05.004}
"brain"

#' Sample meta data for the gutbrain
#'
#'  Sample meta data for the gutbrain data
#'
#' @format A data frame  with 30 rows and 3 variables:
#' \describe{
#'   \item{SRA_id}{The SRA identifier for that specific sample}
#'   \item{geo_accession}{The GEO identifier for the concatenated sample}
#'   \item{source_name}{Human host identifier for human who donated microbiome sample}
#'   \item{title}{Human host identifier + mouse identifier}
#'   \item{bio_sample}{Human host identifier + mouse identifier, codified}
#'   \item{host_donor}{Human host identifier for human who donated microbiome sample based on renamed host metadata from original study}
#'   \item{host_donor_tissue}{Human host identifier + tissue identifier, codified}
#'   \item{development}{ASD or typically developing. This data contains only typically developing (TD) samples}
#'   \item{tissue}{Prefrontal cortex or striatum tissue sample}
#'   \item{condition}{Combination identifier of development and tissue}
#'   \item{t_replicate}{Technical replicate. Either A or B, to denote between specific samples}
#' }
#' @source \url{doi:10.1016/j.cell.2019.05.004}
"brain_metadata"

