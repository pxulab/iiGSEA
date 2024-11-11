#' iiGSEA: Iterative Gene Set Enrichment Analysis for scRNA-seq
#'
#' This function performs an iterative gene set enrichment analysis on single-cell RNA-seq data.
#' The analysis starts with an initial imputation (optional) followed by single-sample GSEA (ssGSEA),
#' and subsequently refines gene sets by identifying differentially expressed genes (DEGs) in significant cells.
#'
#' @param normalized_data A matrix of normalized expression data (genes x cells).
#' @param gs A list of gene sets for the initial analysis.
#' @param imputed_data Optional matrix of imputed expression data (genes x cells); if not provided, MAGIC imputation is applied.
#' @param iteration Logical, whether iterative refinement of gene sets should be performed.
#' @param n_cutoff Integer, the maximum number of top markers to refine.
#' @param start Numeric, the start position for defining significant cells.
#' @param end Numeric, the end position for defining significant cells.
#' @param magic_seed Integer, a seed for reproducibility of MAGIC imputation.
#' @param jobs Integer, the number of jobs/threads for parallel processing.
#' @return A list containing:
#'         - EnrichmentScores: The final enrichment scores for each cell and gene set.
#'         - Cutoff: Calculated cutoffs for significant cells.
#'         - Markers: A data frame of differentially expressed genes if iteration is TRUE.
#' @export
#' @examples
#' # Example usage
#' iiGSEA(normalized_data = matrix(rnorm(1000), nrow=100, ncol=10), gs = list(GeneSet1 = c("GeneA", "GeneB")))
iiGSEA <- function(
    normalized_data,
    gs,
    imputed_data = NA,
    iteration = TRUE,
    n_cutoff = 50,
    start = 0.01,
    end = 0.6,
    magic_seed = 123,
    jobs = 20
) {
  # Load necessary libraries
  library(Rmagic)
  library(reticulate)

  # Step 1: Impute missing expression data (if imputed data is not provided)
  if (is.null(imputed_data) || (length(imputed_data) == 1 && is.na(imputed_data))) {
    message("Applying MAGIC imputation...")
    datExpr_magic <- Rmagic::magic(t(normalized_data), t = 3, genes = "all_genes", seed = magic_seed, n.jobs = jobs)
    datExpr <- t(datExpr_magic$result)
  } else {
    datExpr <- as.matrix(imputed_data)
  }
  
  # Step 2: Run initial ssGSEA on expression data
  es_out <- run_gsea(datExpr, gs)
  
  # Step 3: Calculate cutoffs for significant cells in the initial ssGSEA results
  es_cutoff <- run_cutoff(es_out, start, end)
  
  # Step 4: Perform gene set refinement if iteration is TRUE
  if (iteration) {
    message("Refining gene sets based on initial ssGSEA results...")
    
    # Step 4a: Refine gene sets based on the first iteration of ssGSEA
    gs2 <- gs_generate(es_out, es_cutoff, normalized_data, n_cutoff, start, end)
    
    # Step 4b: Run a second iteration of ssGSEA with the refined gene sets
    iigsea_es <- run_gsea(datExpr, gs2[[1]])
    iigsea_es_cutoff <- run_cutoff(iigsea_es)
    
    # Step 5: Return results
    return(list(
      "EnrichmentScores" = iigsea_es,
      "Cutoff" = iigsea_es_cutoff,
      "Markers" = gs2[[2]]
    ))
  } else {
    # Return results without refinement
    return(list(
      "EnrichmentScores" = es_out,
      "Cutoff" = es_cutoff,
      "Markers" = ""
    ))
  }
}
