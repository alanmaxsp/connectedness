# Release checklist

Use this checklist before tagging or publishing a new release.

## 1) Regenerate bindings and docs

```r
Rcpp::compileAttributes()
devtools::document()
```

Or run:

```bash
Rscript scripts/update_docs.R
```

## 2) Run checks and tests

```r
devtools::test()
devtools::check()
```

Or via command line:

```bash
R CMD build .
R CMD check connectedness_*.tar.gz
```

## 3) Verify package documentation artifacts

- `NAMESPACE` is up to date.
- `man/*.Rd` exists and reflects current function signatures.
- `help(package = "connectedness")` lists expected user-facing pages.

## 4) Verify metadata and changelog

- `DESCRIPTION` version updated.
- `NEWS.md` updated with release highlights and fixes.
- `URL` and `BugReports` are valid.

## 5) Install and smoke test

```r
devtools::install()
```

Quick smoke checks:

- `build_Ainv()`
- `build_Ginv()`
- `build_Hinv()`
- `compute_connectedness()`

## 6) Tag and publish

- Commit generated files (`NAMESPACE`, `man/`, `RcppExports.*`, docs updates).
- Tag release (`vX.Y.Z`).
- Publish GitHub release notes using `NEWS.md` as base.
