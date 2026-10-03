# Desarrollo / Development

<p align="center">
  <a href="#espanol">Español</a> · <a href="#english">English</a>
</p>

<a id="espanol"></a>

## Español

### Ramas y publicación de versiones

| Referencia | Uso |
|---|---|
| `main` | Línea estable para uso habitual; sólo recibe cambios revisados |
| `develop` | Cambios y pruebas para la próxima versión |
| `vX.Y.Z` | Etiqueta de una versión publicada que se conserva sin cambios |

El nombre de una rama no certifica por sí mismo la calidad del código. Esta
separación es una política de trabajo; no se han configurado reglas de protección
que impidan escribir directamente en `main`. Mantenga `main` como rama
predeterminada del repositorio.

Para el trabajo habitual en GitHub Desktop:

1. Abra o clone el repositorio y use **Fetch origin** para consultar novedades.
   Aplique **Pull origin** si hay cambios pendientes. Antes de cambiar de rama,
   guarde sus cambios mediante un commit en su rama de trabajo.
2. En **Current branch**, seleccione `develop`. Para una tarea aislada puede crear
   una rama nueva desde `develop`, por ejemplo `feature/nueva-funcion`, e
   incorporarla posteriormente a `develop` mediante una pull request.
3. Edite y pruebe allí. Revise **Changes**, haga el commit y use **Push origin**.
   Esto actualiza desarrollo; no actualiza `main` ni versiones ya publicadas.
4. Cuando el cambio esté listo, abra una **pull request** con destino `main` y
   origen `develop`. Revise el diff y los resultados de las comprobaciones.
   Las pruebas aprobadas no sustituyen la revisión metodológica.
5. Incorpore la pull request sólo después de revisar esos resultados. Esto cambia
   lo que reciben quienes instalan desde `main`. **No elimine `develop`** después
   de incorporarla: es la rama de desarrollo permanente.
6. Para publicar una versión, actualice `DESCRIPTION`, `NEWS.md` y las referencias
   de instalación/cita necesarias dentro del cambio revisado. En **Releases**,
   cree una publicación con etiqueta `vX.Y.Z` sobre el commit validado de `main`.
   No mueva ni reemplace etiquetas anteriores; publique una nueva para corregirlas.
7. Después de publicar, sincronice `develop` incorporando los cambios de `main`
   sin sobrescribir trabajos pendientes. Durante desarrollo puede usar una versión
   como `X.Y.Z.9000` en `DESCRIPTION`; la versión publicada debe volver a una
   numeración de publicación apropiada.

Una release identifica y describe una versión del código; no implica publicación
en CRAN ni creación automática de un instalador binario de R. Crear `develop`
tampoco crea una release. Antes de una primera publicación numerada, las
instalaciones por etiqueta que muestra el README son ejemplos pendientes.

Las ramas y etiquetas pertenecen al mismo repositorio. Cuando éste es público,
desarrollo también lo es. La separación controla qué versión se instala por
defecto, no oculta el código en desarrollo. Instalar otra rama en la misma
biblioteca de R reemplaza el paquete instalado; para comparar instalaciones
simultáneas use bibliotecas separadas.

Fuentes: [ramas de GitHub](https://docs.github.com/en/pull-requests/reference/branches),
[releases de GitHub](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
y [selección de referencias en remotes](https://remotes.r-lib.org/reference/install_github.html).

### Comprobaciones y web de documentación

R-CMD-check se ejecuta en pushes y pull requests dirigidas a `main` o
`develop`, y mediante ejecución manual. Los errores del check hacen
fallar el trabajo; las advertencias y notas deben revisarse antes de
incorporar o publicar. El sitio pkgdown se actualiza desde `main`; los
pushes a `develop` no publican la web estable. Este sitio puede generar
un commit adicional de documentación en `main`: incorpórelo después
a `develop` sin sobrescribir cambios pendientes.

---

<a id="english"></a>

## English

### Branches and version releases

| Reference | Purpose |
|---|---|
| `main` | Stable line for normal use; receives reviewed changes only |
| `develop` | Changes and tests for the next version |
| `vX.Y.Z` | Published version tag retained unchanged |

A branch name does not certify code quality by itself. This separation is a
working policy; no protection rules have been configured to prevent direct writes
to `main`. Keep `main` as the repository's default branch.

For routine work in GitHub Desktop:

1. Open or clone the repository and use **Fetch origin** to check for updates.
   Use **Pull origin** if updates are pending. Before switching branches, save
   your changes in a commit on your working branch.
2. Under **Current branch**, select `develop`. For an isolated task, create a new
   branch from `develop`, such as `feature/new-function`, and later merge it into
   `develop` through a pull request.
3. Edit and test there. Review **Changes**, commit and use **Push origin**. This
   updates development; it does not update `main` or published versions.
4. When ready, open a **pull request** with `main` as base and `develop` as head.
   Review the diff and check results. Passing tests do not replace methodological
   review.
5. Merge only after reviewing those results. This changes what users installing
   from `main` receive. **Do not delete `develop`** after merging: it is the
   permanent development branch.
6. To release a version, update `DESCRIPTION`, `NEWS.md` and any required
   installation/citation references within the reviewed change. Under
   **Releases**, create a release with tag `vX.Y.Z` on the validated `main` commit.
   Do not move or replace previous tags; publish a new version to fix them.
7. After release, synchronize `develop` by merging changes from `main` without
   overwriting pending work. During development, a version such as `X.Y.Z.9000`
   may be used in `DESCRIPTION`; the published version must return to an
   appropriate release version number.

A release identifies and describes a code version; it does not imply CRAN
publication or automatically create an R binary installer. Creating `develop`
does not create a release either. Until the first numbered release, the README's
tag installation commands are pending examples.

Branches and tags belong to the same repository. When it is public, development
is public too. This separation controls the default installation version; it
does not hide development code. Installing another branch in the same R library
replaces the installed package; use separate libraries to compare simultaneous
installations.

Sources: [GitHub branches](https://docs.github.com/en/pull-requests/reference/branches),
[GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
and [reference selection in remotes](https://remotes.r-lib.org/reference/install_github.html).

### Checks and documentation website

R-CMD-check runs on pushes and pull requests targeting `main` or
`develop`, and on manual dispatch. Check errors fail the job; warnings
and notes must be reviewed before merging or releasing. The pkgdown
website updates from `main`; pushes to `develop` do not publish the
stable website. Site generation may create an additional documentation
commit on `main`: merge it back into `develop` without overwriting
pending changes.
