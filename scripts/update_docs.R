#!/usr/bin/env Rscript

# Regenerate package interfaces and docs from source comments.
# Run from repository root:
#   Rscript scripts/update_docs.R

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required.")
}

Rcpp::compileAttributes()
devtools::document()

cat("\nDocumentation refresh complete.\n")
cat("Remember to commit:\n")
cat("  - NAMESPACE\n")
cat("  - man/*.Rd\n")
cat("  - R/RcppExports.R\n")
cat("  - src/RcppExports.cpp\n\n")
