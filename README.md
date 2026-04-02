# connectedness

An R package for computing connectedness metrics between management units (MUs)
in animal breeding genetic evaluations.

## Metrics

- **CD** (Coefficient of Determination): ranges from 0 (disconnected) to 1 (fully connected).
- **PEVD** (Prediction Error Variance of Differences): in units of residual variance.

Both metrics are based on the contrast approach via Mixed Model Equations
(Laloë, 1993; Laloë et al., 1996), implemented with sparse matrix arithmetic
via RcppEigen for scalability to large pedigrees.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yourusername/connectedness")
```

Requires a C++ compiler (Rtools on Windows, Xcode CLT on macOS).

## Quick start

```r
library(connectedness)

res <- compute_connectedness(
  data          = my_data,
  pedigree      = my_pedigree,   # complete pedigree including ancestors
  animal_col    = "animal_id",
  mu_col        = "region",
  fixed_formula = ~ 1 + sex + birth_year,
  sigma2a       = 5.66,          # from prior variance components analysis
  sigma2e       = 10.24,
  year_col      = "birth_year",
  year_window   = c(2003, 2022)
)

print(res)
plot(res, which = "all")
```

## References

Laloë, D. (1993). Precision and information in linear models of genetic
evaluation. *Genetics Selection Evolution*, 25, 557–576.

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on measures
of precision and connectedness in mixed linear models of genetic evaluation.
*Genetics Selection Evolution*, 28, 359–378.
