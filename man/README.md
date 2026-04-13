# man directory

This directory is intentionally versioned and should contain generated
`*.Rd` help files created by:

```r
Rcpp::compileAttributes()
devtools::document()
```

Before release, ensure this directory includes up-to-date help pages and commit
the generated files together with `NAMESPACE`, `R/RcppExports.R`, and
`src/RcppExports.cpp`.
