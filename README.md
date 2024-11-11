```markdown
# iiGSEA

[![Build Status](https://github.com/yourusername/iiGSEA/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/iiGSEA/actions)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/iiGSEA)](https://CRAN.R-project.org/package=iiGSEA)
[![Coverage Status](https://coveralls.io/repos/github/yourusername/iiGSEA/badge.svg?branch=main)](https://coveralls.io/github/yourusername/iiGSEA?branch=main)

The iiGSEA package is designed for performing iteractive imputed gene set enrichment analysis (GSEA) on normalized gene expression data. It includes functions to compute enrichment scores and identify marker genes associated with specific gene sets.

## Installation

You can install the latest stable release from CRAN:

```r
devtools::install("iiGSEA")
```

Or the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("iiGSEA")
```

## Usage Example

Below is a quick demonstration of how to use the `iiGSEA` package:

```r
library(iiGSEA)

# Load normalized gene expression data
normalized_data <- readRDS("iiGSEA/data/pancreas_normalized.RDS")

# Read gene set information
gene_set <- read.delim("iiGSEA/data/input_gene_set.tsv")

# Prepare gene sets for analysis
gs <- list()
for (i in unique(gene_set$marker_group)) {
  target_gs <- as.list(gene_set[gene_set$marker_group == i, "gene"])
  names(target_gs) <- i
  gs <- c(gs, target_gs)
}

# Perform GSEA
iiGSEA_result <- iiGSEA(normalized_data, gs)
iiGSEA_es <- iiGSEA_result[[1]]
iiGSEA_markers <- iiGSEA_result[[2]]
```

(Optional) Run iiGSEA with imputed expression matrix, without iteration.

```r
iiGSEA_result = iiGSEA(normalized_data,gs,imputed_data = imputed_data,iteration = FALSE)
```


## Contributing

Contributions are welcome! If you find bugs or want to contribute new features, please submit an issue or pull request on GitHub.

## Citation

If you use this package in your research, please cite it as follows:

```
BibTeX entry:
@Manual{,
  title = {iiGSEA: Gene Set Enrichment Analysis},
  author = {Your Name},
  year = {2024},
  note = {URL: https://github.com/yourusername/iiGSEA},
}
```

## License

This project is licensed under the [MIT License](LICENSE).
```
