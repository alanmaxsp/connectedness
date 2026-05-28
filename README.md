# connectedness

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

`connectedness` es un paquete de R para calcular conectividad genética entre
unidades de manejo (MUs, por sus siglas en inglés) en evaluaciones genéticas
animales.

Implementa **métricas de conectividad basadas en contrastes** a partir de las
ecuaciones de modelos mixtos (MME) y permite análisis basados en:

* relaciones de pedigree mediante **A⁻¹**
* relaciones genómicas mediante **G⁻¹**
* relaciones combinadas pedigree-genómicas mediante **H⁻¹**
* kernels inversos provistos por el usuario

### ¿Qué calcula?

El paquete actualmente provee dos métricas de conectividad entre
unidades de manejo:

* **Contraste CD**: coeficiente de determinación de contrastes entre MUs
* **Contraste PEVD**: varianza del error de predicción de diferencias entre MUs

Valores más altos de CD y más bajos de PEVD indican mayor conectividad.

### Instalación

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

Se requiere una herramienta de compilación C++ funcional:

* **Windows**: Rtools
* **macOS**: herramientas de línea de comandos de Xcode
* **Linux**: toolchain estándar de compilación

Para una introducción más extensa con ejemplos desarrollados, ver la
[vignette introductoria](https://alanmaxsp.github.io/connectedness/intro.html).

### Antes de empezar

Para ejecutar `compute_connectedness()`, los datos deben incluir:

* una columna identificadora del animal (`animal_col`)
* una columna de unidad de manejo (`mu_col`)
* todas las variables de efectos fijos usadas en `fixed_formula`

Dependiendo de la estructura de relaciones, también se necesita:

* **Ainv**: un pedigree con animal, padre y madre
* **Ginv**: una matriz de genotipos `X` y un `animal_index`
* **Hinv**: un pedigree, una matriz de genotipos `X` y `genotyped_idx`
  (índices de los animales genotipados en el pedigree renumerado)

Cada animal debe pertenecer a **una y solo una** unidad de manejo.

### Ejemplos mínimos

#### Conectividad basada en pedigree (Ainv)

```r
library(connectedness)

res_A <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree
)
```

#### Conectividad genómica (Ginv)

```r
res_G <- compute_connectedness(
  data          = my_genotyped_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ginv",
  X             = my_genotypes_matrix,
  animal_index  = my_index
)
```

#### Conectividad combinada pedigree-genómica (Hinv)

```r
res_H <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Hinv",
  pedigree      = my_pedigree,
  X             = my_genotypes_matrix,
  genotyped_idx = my_genotyped_idx
)
```

### Salida

`compute_connectedness()` devuelve un objeto de clase `"connectedness"` con
componentes como:

* `CD`: matriz o resumen de valores de conectividad basados en CD entre unidades
  de manejo.
* `PEVD`: matriz o resumen de valores de conectividad basados en PEVD entre
  unidades de manejo.
* `n_target`: número de animales objetivo usados en el análisis, cuando
  corresponde.
* `relationship`: estructura de relaciones usada en el análisis (`"Ainv"`,
  `"Ginv"`, `"Hinv"` o custom).

Los resultados se pueden inspeccionar o visualizar con:

```r
print(res_A)
plot(res_A, which = "all")
```

### Restricción temporal opcional

La conectividad también puede evaluarse dentro de una ventana temporal
restringida:

```r
res_time <- compute_connectedness(
  data                 = my_data,
  animal_col           = "animal_id",
  mu_col               = "herd",
  fixed_formula        = ~ 1 + herd + sex,
  sigma2a              = 2.0,
  sigma2e              = 5.0,
  relationship         = "Ainv",
  pedigree             = my_pedigree,
  year_col             = "birth_year",
  year_window          = c(2018, 2022),
  min_records_per_year = 10
)

plot(res_time, which = "overlap")
```

### Funciones principales

* `compute_connectedness()`
* `build_Ainv()`
* `build_Ginv()`
* `build_Hinv()`

### Referencias

Kennedy, B. W., & Trus, D. (1993). Considerations on genetic connectedness
between management units under an animal model. *Journal of Animal Science*,
71, 2341–2352.

Laloë, D. (1993). Precision and information in linear models of genetic
evaluation. *Genetics Selection Evolution*, 25, 557–576.

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on measures of
precision and connectedness in mixed linear models of genetic evaluation.
*Genetics Selection Evolution*, 28, 359–378.

Yu, H., & Morota, G. (2021). GCA: An R package for genetic connectedness analysis using pedigree and genomic data. *BMC Genomics*, 22, 119.

---

<a id="english"></a>

## English

`connectedness` is an R package for computing genetic connectedness between
management units (MUs) in animal genetic evaluations.

It implements **contrast-based connectedness metrics** from the mixed model equations (MME) and supports analyses based on:

* pedigree relationships through **A⁻¹**
* genomic relationships through **G⁻¹**
* combined pedigree-genomic relationships through **H⁻¹**
* user-supplied inverse kernels

### What does it compute?

The package currently provides two pairwise connectedness metrics between management units:

* **CD contrast**: Coefficient of Determination of contrasts between MUs
* **PEVD contrast**: Prediction Error Variance of Differences between MUs

Higher CD and lower PEVD indicate stronger connectedness.

### Installation

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

A working C++ toolchain is required:

* **Windows**: Rtools
* **macOS**: Xcode command line tools
* **Linux**: standard compiler toolchain

For a longer introduction with worked examples, see the [intro vignette](https://alanmaxsp.github.io/connectedness/intro.html).

### Before you start

To run `compute_connectedness()`, your data should include:

* an animal identifier column
* a management-unit column (`mu_col`)
* all fixed-effect variables used in `fixed_formula`

Depending on the relationship structure, you will also need:

* **Ainv**: a pedigree with animal, sire, and dam
* **Ginv**: a genotype matrix `X` and an `animal_index`
* **Hinv**: a pedigree, a genotype matrix `X`, and `genotyped_idx`

Each animal must belong to **one and only one** management unit.

### Minimal examples

#### Pedigree-based connectedness (Ainv)

```r
library(connectedness)

res_A <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree
)
```

#### Genomic connectedness (Ginv)

```r
res_G <- compute_connectedness(
  data          = my_genotyped_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ginv",
  X             = my_genotypes_matrix,
  animal_index  = my_index
)
```

#### Combined pedigree-genomic connectedness (Hinv)

```r
res_H <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Hinv",
  pedigree      = my_pedigree,
  X             = my_genotypes_matrix,
  genotyped_idx = my_genotyped_idx
)
```

### Output

`compute_connectedness()` returns an object of class `"connectedness"` with
components such as:

* `CD`: matrix or summary of CD-based connectedness values between management
  units.
* `PEVD`: matrix or summary of PEVD-based connectedness values between
  management units.
* `n_target`: number of target animals used in the analysis, when applicable.
* `relationship`: relationship structure used in the analysis (`"Ainv"`,
  `"Ginv"`, `"Hinv"`, or custom).

You can inspect or visualize results with:

```r
print(res_A)
plot(res_A, which = "all")
```

### Optional temporal restriction

Connectedness can also be evaluated within a restricted time window:

```r
res_time <- compute_connectedness(
  data                 = my_data,
  animal_col           = "animal_id",
  mu_col               = "herd",
  fixed_formula        = ~ 1 + herd + sex,
  sigma2a              = 2.0,
  sigma2e              = 5.0,
  relationship         = "Ainv",
  pedigree             = my_pedigree,
  year_col             = "birth_year",
  year_window          = c(2018, 2022),
  min_records_per_year = 10
)

plot(res_time, which = "overlap")
```

### Main functions

* `compute_connectedness()`
* `build_Ainv()`
* `build_Ginv()`
* `build_Hinv()`

### References

Kennedy, B. W., & Trus, D. (1993). Considerations on genetic connectedness
between management units under an animal model. *Journal of Animal Science*,
71, 2341–2352.

Laloë, D. (1993). Precision and information in linear models of genetic
evaluation. *Genetics Selection Evolution*, 25, 557–576.

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on measures of
precision and connectedness in mixed linear models of genetic evaluation.
*Genetics Selection Evolution*, 28, 359–378.

Yu, H., & Morota, G. (2021). GCA: An R package for genetic connectedness analysis using pedigree and genomic data. *BMC Genomics*, 22, 119.
