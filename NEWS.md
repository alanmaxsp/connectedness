# connectedness 0.1.0

## Features

- Added core user workflows for pedigree, genomic, and H-kernel connectedness:
  `build_Ainv()`, `build_Ginv()`, `build_Hinv()`, and `compute_connectedness()`.
- Added support for contrast-based connectedness metrics:
  CD contrast and PEVD contrast.
- Added temporal overlap summaries in `compute_connectedness()`.

## Improvements

- Improved robustness of genomic inverse construction by excluding SNPs with
  zero observed genotypes from MAF filtering and downstream use.
- Removed LDLT fallback in `compute_Ginv_cpp()` and now stop with a clear
  regularization message if LLT fails.
- Harmonized package terminology to describe `H` as a combined
  pedigree-genomic kernel and removed ssGBLUP wording from user-facing docs.
- Marked low-level C++ entry points as internal in roxygen
  (`@keywords internal`, `@noRd`) to keep the public API focused.

## Documentation

- Refreshed README wording and examples.
- Added references for H-kernel context (including Legarra et al., 2009).
