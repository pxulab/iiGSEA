#' run_gsea: Single-sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' This function performs ssGSEA to calculate enrichment scores for specified gene sets in single-cell RNA-seq data.
#' The method is particularly useful for identifying gene set activity within individual cells.
#'
#' @param datExpr A matrix of expression data (genes x cells) to analyze.
#' @param gs A list containing gene sets, where each element is a vector of gene names representing a single gene set.
#' @return A data frame with enrichment scores (cells x gene sets), where each cell's enrichment score
#'         corresponds to the activity of each gene set in that cell.
#' @export
#' @examples
#' # Example data and gene sets
#' datExpr <- matrix(rnorm(1000), nrow=100, ncol=10)
#' gs <- list(GeneSet1 = c("GeneA", "GeneB", "GeneC"), GeneSet2 = c("GeneD", "GeneE"))
#' run_gsea(datExpr, gs)
#' 
run_gsea <- function(datExpr, gs) {
    library(GSVA)
  
	  ssgseaPar <- ssgseaParam(datExpr, gs, normalize = FALSE)
    ssgsea.es <- gsva(ssgseaPar)
    return(as.data.frame(t(ssgsea.es)))
}


