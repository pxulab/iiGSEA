#' gs_generate: Refine Gene Sets Based on ssGSEA Scores
#'
#' This function refines input gene sets based on ssGSEA scores and differential expression analysis.
#' Gene sets are iteratively adjusted by identifying marker genes associated with significant cells
#' in each enrichment score set.
#'
#' @param es_out A data frame of enrichment scores (cells x gene sets) from \code{run_gsea}.
#' @param es_cutoff A data frame containing cutoffs for enrichment scores, used to define significant cells.
#' @param normalized_data A matrix of normalized expression data (genes x cells) for calculating marker genes.
#' @param n_cutoff Integer, the maximum number of top markers to select when refining each gene set.
#' @param start Numeric, the starting proportion for defining significant enrichment score cells.
#' @param end Numeric, the ending proportion for defining significant enrichment score cells.
#' @return A list where the first element is the refined gene sets, and the second element is a data frame
#'         of differentially expressed genes for each gene set.
#' @export
#' @examples
#' # Example data
#' es_out <- data.frame(GeneSet1 = runif(10), GeneSet2 = runif(10))
#' es_cutoff <- data.frame(GeneSet1 = 5, GeneSet2 = 4)
#' normalized_data <- matrix(rnorm(1000), nrow=100, ncol=10)
#' gs_generate(es_out, es_cutoff, normalized_data)
#' 
gs_generate <- function(es_out, es_cutoff, normalized_data, n_cutoff = 50, start = 0.01, end = 0.6) {
  library(dplyr)
  library(Seurat)
  
  obj <- CreateSeuratObject(counts = normalized_data, project = "iiGSEA", min.cells = 0, min.features = 0)
  obj[["RNA"]]$data <- obj[["RNA"]]$counts
  gs <- list()
  DEG_out <- data.frame()
  
  for (i in colnames(es_out)) {
    gsea <- es_out[, i, drop = FALSE]
    gsea <- gsea[order(gsea[, 1], decreasing = TRUE), , drop = FALSE]
    
    cutoff <- es_cutoff[, i]
    obj@meta.data$iigsea <- ifelse(colnames(obj) %in% row.names(gsea)[1:(cutoff - 1)], 0, 1)
    
    DEG <- FindMarkers(obj, ident.1 = 0, group.by = "iigsea", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
    DEG$gene <- row.names(DEG)
    DEG$type <- i
    DEG_out <- rbind(DEG_out, DEG)
    
    DEG_sig <- DEG[DEG$avg_log2FC > 0.25 & DEG$pct.1 > 0.5 & DEG$p_val_adj < 0.05, ]
    target_gs <- list(row.names(DEG_sig[1:min(nrow(DEG_sig), n_cutoff), ]))
    names(target_gs) <- i
    gs <- c(gs, target_gs)
  }
  
  return(list(gs, DEG_out))
}
