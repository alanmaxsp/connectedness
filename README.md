# connectedness

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

`connectedness` es un paquete de R para calcular conectividad genética entre
unidades de manejo (MUs, por sus siglas en inglés) en evaluaciones genéticas
animales.

Implementa métricas de conectividad basadas en contrastes a partir de las
ecuaciones de modelos mixtos (MME) y permite usar relaciones de pedigree
(**A⁻¹**), genómicas (**G⁻¹**), combinadas pedigree-genómicas (**H⁻¹**) o kernels
inversos definidos por el usuario.

### Instalación

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

Se requiere una herramienta de compilación C++ funcional:

* **Windows**: Rtools
* **macOS**: herramientas de línea de comandos de Xcode
* **Linux**: toolchain estándar de compilación

### Versión estable, desarrollo y versiones numeradas

La instalación anterior usa `main`, la rama reservada para cambios revisados.
`develop` es la rama para cambios y pruebas de la próxima versión. Sus cambios
no modifican `main` hasta que se revisan y se incorporan explícitamente.

Para probar la versión de desarrollo, elíjala expresamente:

```r
remotes::install_github("alanmaxsp/connectedness@develop")
```

Para reproducir un análisis con una versión numerada, use su etiqueta. Este
ejemplo será válido **cuando publiquemos la etiqueta `v0.1.0`**:

```r
remotes::install_github("alanmaxsp/connectedness@v0.1.0")
```

Todavía no hay versiones numeradas publicadas. Una rama puede avanzar; una
etiqueta publicada se conserva sin cambios y las correcciones reciben una nueva
versión. Registre `packageVersion("connectedness")` y `sessionInfo()`; si usa una rama,
registre además el commit instalado (`packageDescription("connectedness")$RemoteSha`).

En un repositorio público, `develop` también es visible y descargable, pero sólo
se instala si el usuario la selecciona. Instalar estable y desarrollo en la misma
biblioteca de R reemplaza la instalación anterior; use bibliotecas separadas si
necesita conservar ambas simultáneamente.

La [guía de trabajo](CONTRIBUTING.md#espanol) explica cómo cambiar de rama en GitHub Desktop
y revisar cambios antes de incorporarlos a `main`.

### ¿Qué calcula?

`compute_connectedness()` devuelve dos métricas entre pares de MUs:

* **Contraste CD**: coeficiente de determinación de contrastes entre MUs.
* **Contraste PEVD**: varianza del error de predicción de diferencias entre MUs.

Valores más altos de CD y más bajos de PEVD indican mayor conectividad.

### Ejemplo mínimo con pedigree (Ainv)

```r
library(connectedness)

res <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree
)

print(res)
plot(res, which = "all")
```

El paquete también soporta `relationship = "Ginv"`, `"Hinv"` y `"custom"`.
Para ejemplos desarrollados, ver la
[vignette introductoria](https://alanmaxsp.github.io/connectedness/articles/intro.html).

### Selección temporal de animales target

Una ventana temporal puede usarse para seleccionar MUs activas y definir los
animales target de los contrastes. Por defecto (`target_scope = "window"`), el
MME se ajusta usando todos los registros disponibles en `data`, pero CD/PEVD se
reportan para animales de MUs activas dentro de la ventana.

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
  min_records_per_year = 30
)
```

### Diagnóstico rápido

Para bases grandes, `dry_run = TRUE` permite inspeccionar el tamaño esperado del
sistema antes de resolver las MME, incluido si conviene usar el método
`schur_solver = "cholmod_lowmem"` para evitar matrices densas intermedias grandes:

```r
diag <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree,
  dry_run       = TRUE
)
```

### Salida principal

El objeto `connectedness` incluye, entre otros componentes:

* `CD` y `PEVD`: matrices de conectividad entre MUs reportadas.
* `n_target`: número de animales target por MU reportada; estos animales reciben
  pesos distintos de cero en los contrastes pareados.
* `report_mus`: MUs incluidas en las matrices CD/PEVD.
* `target_scope`: definición de los animales target usados en los contrastes.

### Más información

La [vignette introductoria](https://alanmaxsp.github.io/connectedness/articles/intro.html)
es el documento principal para la explicación metodológica, ejemplos con
`Ginv`, `Hinv` y kernels custom, diagnóstico computacional y referencias.

### Funciones principales

* `compute_connectedness()`
* `build_Ainv()`
* `build_Ginv()`
* `build_Hinv()`

---

<a id="english"></a>

## English

`connectedness` is an R package for computing genetic connectedness between
management units (MUs) in animal genetic evaluations.

It implements contrast-based connectedness metrics from the mixed model
equations (MME) and supports pedigree relationships (**A⁻¹**), genomic
relationships (**G⁻¹**), combined pedigree-genomic relationships (**H⁻¹**), and
user-supplied inverse kernels.

### Installation

```r
# install.packages("remotes")
remotes::install_github("alanmaxsp/connectedness")
```

A working C++ toolchain is required:

* **Windows**: Rtools
* **macOS**: Xcode command line tools
* **Linux**: standard compiler toolchain

### Stable, development and numbered versions

The installation above uses `main`, the branch reserved for reviewed changes.
`develop` is the branch for changes and tests for the next version. Its changes
do not modify `main` until they are reviewed and explicitly merged.

To try the development version, select it explicitly:

```r
remotes::install_github("alanmaxsp/connectedness@develop")
```

To reproduce an analysis with a numbered version, use its tag. This example will
work **once we publish the `v0.1.0` tag**:

```r
remotes::install_github("alanmaxsp/connectedness@v0.1.0")
```

No numbered versions have been published yet. A branch may advance; a published
tag is retained unchanged and fixes receive a new version. Record
`packageVersion("connectedness")` and `sessionInfo()`; when using a branch, also record the
installed commit (`packageDescription("connectedness")$RemoteSha`).

In a public repository, `develop` is also visible and downloadable, but is only
installed when users select it. Installing stable and development in the same R
library replaces the previous installation; use separate libraries if you need
both simultaneously.

The [working guide](CONTRIBUTING.md#english) explains switching branches in GitHub Desktop
and reviewing changes before merging them into `main`.

### What does it compute?

`compute_connectedness()` returns two pairwise metrics between MUs:

* **CD contrast**: coefficient of determination of contrasts between MUs.
* **PEVD contrast**: prediction error variance of differences between MUs.

Higher CD and lower PEVD indicate stronger connectedness.

### Minimal pedigree example (Ainv)

```r
library(connectedness)

res <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree
)

print(res)
plot(res, which = "all")
```

The package also supports `relationship = "Ginv"`, `"Hinv"`, and `"custom"`.
See the [intro vignette](https://alanmaxsp.github.io/connectedness/articles/intro.html)
for worked examples.

### Temporal definition of target animals

A time window can be used to select active MUs and define the target animals for
the contrasts. By default (`target_scope = "window"`), the MME is fitted using
all records in `data`, but CD/PEVD are reported for animals from active MUs
inside the time window.

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
  min_records_per_year = 30
)
```

### Quick diagnostics

For large datasets, use `dry_run = TRUE` to inspect the expected MME system size
before solving, including whether `schur_solver = "cholmod_lowmem"` should be
used to avoid large intermediate dense matrices:

```r
diag <- compute_connectedness(
  data          = my_data,
  animal_col    = "animal_id",
  mu_col        = "herd",
  fixed_formula = ~ 1 + herd + sex,
  sigma2a       = 2.0,
  sigma2e       = 5.0,
  relationship  = "Ainv",
  pedigree      = my_pedigree,
  dry_run       = TRUE
)
```

### Main output

The `connectedness` object includes, among other components:

* `CD` and `PEVD`: connectedness matrices among reported MUs.
* `n_target`: number of target animals per reported MU; these animals receive
  non-zero weights in the pairwise contrasts.
* `report_mus`: MUs included in the reported CD/PEVD matrices.
* `target_scope`: definition of the target animals used in the contrasts.

### More information

The [intro vignette](https://alanmaxsp.github.io/connectedness/articles/intro.html) is
the main document for methodological background, examples with `Ginv`, `Hinv`
and custom kernels, computational diagnostics, and references.

### Main functions

* `compute_connectedness()`
* `build_Ainv()`
* `build_Ginv()`
* `build_Hinv()`
