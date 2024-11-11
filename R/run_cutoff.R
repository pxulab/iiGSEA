#' run_cutoff: Determine Enrichment Score Cutoffs
#'
#' Calculates optimal cutoff values within enrichment scores to define significant cells based on
#' changes in enrichment score trends. This helps identify cells with prominent gene set activity.
#'
#' @param es_out A data frame of enrichment scores (cells x gene sets) generated from \code{run_gsea}.
#' @param start Numeric, indicating the starting proportion of enrichment scores to consider when defining cutoffs.
#' @param end Numeric, indicating the ending proportion of enrichment scores to consider when defining cutoffs.
#' @return A data frame containing calculated cutoff values for each gene set. Each cutoff is identified based on
#'         score trends, highlighting the point where the score derivative changes significantly.
#' @export
#' @examples
#' # Example enrichment scores
#' es_out <- data.frame(GeneSet1 = runif(10), GeneSet2 = runif(10))
#' run_cutoff(es_out, start = 0.01, end = 0.6)
#' 
run_cutoff <- function(es_out, start = 0.01, end = 0.6) {
  library(pracma)
  
  es_cutoff <- data.frame(matrix(nrow = 1, ncol = ncol(es_out)), row.names = "cutoff")
  colnames(es_cutoff) <- colnames(es_out)
  
  for (i in colnames(es_out)) {
    gsea <- es_out[, i, drop = FALSE]
    gsea <- gsea[order(gsea[, 1], decreasing = TRUE), , drop = FALSE]
    gsea$n <- 1:nrow(gsea)
    
    derivative <- diff(gsea[, 1]) / diff(gsea[, 2])
    derivative[floor(length(derivative) * end):length(derivative)] <- 0
    derivative[0:floor(length(derivative) * start)] <- 0
    
    peaks <- findpeaks(-derivative, threshold = 0.5, minpeakdistance = nrow(gsea) * 0.02, sortstr = TRUE)
    cutoff <- peaks[1, 2]
    es_cutoff[, i] <- cutoff
  }
  return(es_cutoff)
}
