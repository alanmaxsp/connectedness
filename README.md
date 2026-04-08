# connectedness

`connectedness` is an R package for computing genetic connectedness between
management units (MUs) in animal breeding evaluations.

The package focuses on **contrast-based** connectedness metrics computed from
mixed model equations, and supports connectedness analyses based on:

- **A-inverse** (pedigree-based connectedness),
- **G-inverse** (genomic connectedness),
- **H-inverse** (single-step connectedness), and
- **custom inverse kernels** supplied by the user.

## What problem does it solve?

Animals are often compared across herds, flocks, regions, years, or other
management units. Those comparisons are not equally reliable in all data sets.
When management units are weakly linked genetically, differences in breeding
values across units become less precise.

This package quantifies that problem by computing pairwise connectedness between
management units under the same relationship structure used to define the
analysis.

## Metrics

The package currently reports two complementary metrics:

- **CD contrast**: Coefficient of Determination of contrasts between MUs.
  Higher values indicate stronger connectedness.
- **PEVD contrast**: Prediction Error Variance of Differences for the same
  contrasts. Lower values indicate stronger connectedness.

Both metrics are computed under the **contrast approach** via Mixed Model
Equations (MME) following Laloë (1993) and Laloë et al. (1996).

If desired, `PEVD` can also be returned scaled by the additive genetic variance
through `scale_pevd = TRUE` in `compute_connectedness()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

A working C++ toolchain is required (`Rtools` on Windows, Xcode command line
tools on macOS, or a standard compiler toolchain on Linux).

## Quick start

### Pedigree-based connectedness (A-inverse)

```r
library(connectedness)

res_A <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex + birth_year,
  sigma2a       = 5.66,
  sigma2e       = 10.24,
  relationship  = "Ainv",
  pedigree      = my_pedigree
)

print(res_A)
plot.connectedness(res_A, which = "all")
```

If you want `PEVD / sigma2a` instead of `PEVD`, use:

```r
res_A_scaled <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex + birth_year,
  sigma2a       = 5.66,
  sigma2e       = 10.24,
  relationship  = "Ainv",
  pedigree      = my_pedigree,
  scale_pevd    = TRUE
)
```

### Genomic connectedness (G-inverse)

```r
res_G <- compute_connectedness(
  data          = my_genotyped_data,
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex,
  sigma2a       = 5.66,
  sigma2e       = 10.24,
  relationship  = "Ginv",
  X             = my_genotypes,
  animal_index  = my_index
)
```

### Single-step connectedness (H-inverse)

```r
res_H <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex,
  sigma2a       = 5.66,
  sigma2e       = 10.24,
  relationship  = "Hinv",
  pedigree      = my_pedigree,
  X             = my_genotypes,
  genotyped_idx = my_genotyped_idx
)
```

## Main functions

- `compute_connectedness()`
- `renum_pedigree()`
- `build_Ainv()`
- `build_Ginv()`
- `build_Hinv()`

## Output

`compute_connectedness()` returns an object of class `"connectedness"` with
components such as:

- `CD`
- `PEVD`
- `qK`
- `qC`
- `n_target`
- `relationship`
- `overlap` (when temporal restrictions are used)

## Temporal overlap

The package can optionally restrict analyses to a common time window and report
which MU pairs overlap across years. This is useful when connectedness is not
stable over time.

```r
res_time <- compute_connectedness(
  data                 = my_data,
  animal_col           = "animal_id",
  mu_col               = "region",
  fixed_formula        = ~ 1 + sex,
  sigma2a              = 5.66,
  sigma2e              = 10.24,
  relationship         = "Ainv",
  pedigree             = my_pedigree,
  year_col             = "birth_year",
  year_window          = c(2018, 2022),
  min_records_per_year = 10
)

plot.connectedness(res_time, which = "overlap")
```

## References

Kennedy, B. W., & Trus, D. (1993). Considerations on genetic connectedness
between management units under an animal model. *Journal of Animal Science*,
71, 2341–2352.

Laloë, D. (1993). Precision and information in linear models of genetic
evaluation. *Genetics Selection Evolution*, 25, 557–576.

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on measures of
precision and connectedness in mixed linear models of genetic evaluation.
*Genetics Selection Evolution*, 28, 359–378.
