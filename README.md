# connectedness

`connectedness` is an R package for computing genetic connectedness between
management units (MUs) in animal genetic evaluations.

The package focuses on **contrast-based connectedness metrics** derived from the
mixed model equations (MME), and supports analyses based on:

* **pedigree relationships** through (A^{-1}),
* **genomic relationships** through (G^{-1}),
* **single-step relationships** through (H^{-1}), and
* **custom inverse kernels** supplied by the user.

## What problem does it solve?

In animal breeding, animals are often compared across herds, flocks, regions, countries, or other management units. However, those comparisons are not equally
reliable in all data sets. When management units are weakly linked genetically,
contrasts between units become less accurate, and comparisons of breeding values
across units are less informative.

`connectedness` quantifies pairwise connectedness between management units using
the same relationship structure assumed in the analysis, allowing users to
assess how strongly units are genetically linked under pedigree-based, genomic,
or single-step settings.

## Metrics

The package currently reports two core contrast-based metrics:

* **CD contrast**: Coefficient of Determination of contrasts between
  management units. Higher values indicate stronger connectedness.
* **PEVD contrast**: Prediction Error Variance of Differences for the same
  contrasts. Lower values indicate stronger connectedness.

Both metrics are computed under the **contrast approach** from the mixed model
equations, following Laloë (1993) and Laloë et al. (1996).

## Installation

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

A working C++ toolchain is required (`Rtools` on Windows, Xcode command line
tools on macOS, or a standard compiler toolchain on Linux).

## Quick start

The following examples are schematic and illustrate the main arguments required
by the package.

### Pedigree-based connectedness ((A^{-1}))

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

If you want `PEVD / sigma2a` instead of unscaled `PEVD`, use:

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

### Genomic connectedness ((G^{-1}))

```r
res_G <- compute_connectedness(
  data          = my_genotyped_data,
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex,
  sigma2a       = 5.66,
  sigma2e       = 10.24,
  relationship  = "Ginv",
  X             = my_genotypes_matrix,
  animal_index  = my_index
)
```

### Single-step connectedness ((H^{-1}))

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
  X             = my_genotypes_matrix,
  genotyped_idx = my_genotyped_idx
)
```

## Main functions

* `compute_connectedness()`
* `build_Ainv()`
* `build_Ginv()`
* `build_Hinv()`

## Output

`compute_connectedness()` returns an object of class `"connectedness"` with
components such as:

* `CD`: matrix or summary of CD-based connectedness values between management
  units.
* `PEVD`: matrix or summary of PEVD-based connectedness values between
  management units.
* `n_target`: number of target animals used in the analysis, when applicable.
* `relationship`: relationship structure used in the analysis (`"Ainv"`,
  `"Ginv"`, `"Hinv"`, or custom).

## Temporal overlap

The package can optionally restrict analyses to a common time window and report
which MU pairs overlap across years. This is useful when genetic links between
management units vary over time and connectedness needs to be evaluated within a
restricted temporal window.

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

Yu, H., & Morota, G. (2021). GCA: An R package for genetic connectedness analysis using pedigree and genomic data. *BMC Genomics*, 22, 119.
