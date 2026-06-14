# connectedness development version

## New features

- Added an alternative Schur-complement backend for `compute_connectedness()` via
  `mme_backend = "schur"`. The backend avoids direct factorization of the full
  MME by factorizing `Cuu = Z'Z + lambda Kinv` and absorbing fixed effects
  through a Schur complement.
- Added `schur_solver` selection for the Schur backend. `"auto"` routes sparse
  kernels such as `Ainv` through CHOLMOD via Matrix, switches to
  `"cholmod_lowmem"` when the dense Schur working matrix `Cuu^{-1} Z'X`
  would be too large, routes dense kernels such as `Ginv` through a dense
  compiled solver, and keeps `"eigen_sparse"` as a diagnostic solver for small
  comparisons.
- The Schur backend is now the default `compute_connectedness()` backend; the
  full-MME backend remains available with `mme_backend = "full_mme"`.
- Added `target_scope` to define which animals from temporally selected MUs
  receive non-zero contrast weights. The MME is built with all available data;
  `target_scope = "window"` targets animals inside `year_window`, while
  `target_scope = "selected_mus"` targets all animals from selected MUs.

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
