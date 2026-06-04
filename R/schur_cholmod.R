.cd_force_symmetric_sparse <- function(x, label = "Kinv", tol = 1e-8) {
  if (inherits(x, "symmetricMatrix")) {
    return(methods::as(x, "dsCMatrix"))
  }

  x <- .as_dgCMatrix(x)
  asym <- Matrix::drop0(x - t(x))
  max_abs_x <- if (length(x@x)) max(abs(x@x)) else 0
  max_abs_asym <- if (length(asym@x)) max(abs(asym@x)) else 0
  scale <- max(1, max_abs_x)

  if (max_abs_asym > tol * scale) {
    stop(sprintf(
      paste0(
        "%s must be symmetric. If you supplied only one triangle of a custom ",
        "sparse matrix, pass it as forceSymmetric(Kinv, uplo = 'L' or 'U') ",
        "before calling compute_connectedness()."
      ),
      label
    ), call. = FALSE)
  }

  x <- Matrix::drop0((x + t(x)) * 0.5)
  methods::as(Matrix::forceSymmetric(x, uplo = "L"), "dsCMatrix")
}

.cd_cholesky <- function(x, label) {
  tryCatch(
    Matrix::Cholesky(x, LDL = FALSE, perm = TRUE, super = TRUE),
    error = function(e) {
      stop(sprintf(
        "Schur-CHOLMOD: Cholesky factorization of %s failed: %s",
        label,
        conditionMessage(e)
      ), call. = FALSE)
    }
  )
}

.cd_solve_cholesky <- function(factor, rhs, label) {
  tryCatch(
    solve(factor, rhs, system = "A"),
    error = function(e) {
      stop(sprintf(
        "Schur-CHOLMOD: solve for %s failed: %s",
        label,
        conditionMessage(e)
      ), call. = FALSE)
    }
  )
}

.cd_contrast_mu_schur_cholmod_R <- function(Kinv,
                                            id_rec,
                                            Xsp,
                                            mu_animal,
                                            target,
                                            sigma2a,
                                            sigma2e,
                                            mu_names,
                                            block_size = 16L,
                                            verbose = TRUE) {

  if (!.is_sparse_matrix(Kinv)) {
    stop("'Kinv' must be a sparse Matrix object for schur_solver = 'cholmod'.")
  }
  if (sigma2a <= 0) stop("'sigma2a' must be > 0.")
  if (sigma2e < 0) stop("'sigma2e' must be >= 0.")

  block_size <- as.integer(block_size)
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L) {
    stop("'block_size' must be a positive integer.")
  }

  lambda <- sigma2e / sigma2a
  N <- nrow(Kinv)
  nrec <- length(id_rec)
  U <- length(mu_names)
  p <- ncol(Xsp)

  if (ncol(Kinv) != N) stop("'Kinv' must be square.")
  if (nrow(Xsp) != nrec) stop("'Xsp' must have length(id_rec) rows.")
  if (length(mu_animal) != N) stop("'mu_animal' must have length N.")
  if (length(target) != N) stop("'target' must have length N.")
  if (any(is.na(id_rec)) || any(id_rec < 1L) || any(id_rec > N)) {
    stop("'id_rec' values must be in 1..N.")
  }

  D <- tabulate(id_rec, nbins = N)

  Kinv_sym <- .cd_force_symmetric_sparse(Kinv, label = "Kinv")
  Cuu_sym <- Matrix::drop0(lambda * Kinv_sym + Matrix::Diagonal(n = N, x = as.numeric(D)))
  Cuu_sym <- Matrix::forceSymmetric(Cuu_sym, uplo = "L")

  if (verbose) {
    message(sprintf(
      "Schur-CHOLMOD: Cuu dim = %d x %d, nnz stored = %d",
      nrow(Cuu_sym), ncol(Cuu_sym), Matrix::nnzero(Cuu_sym)
    ))
    message(sprintf(
      "Schur-CHOLMOD: Kinv dim = %d x %d, nnz stored = %d",
      nrow(Kinv_sym), ncol(Kinv_sym), Matrix::nnzero(Kinv_sym)
    ))
  }

  if (verbose) message("Schur-CHOLMOD: factorizing Cuu...")
  fac_Cuu <- .cd_cholesky(Cuu_sym, "Cuu")

  if (verbose) message("Schur-CHOLMOD: factorizing Kinv...")
  fac_Kinv <- .cd_cholesky(Kinv_sym, "Kinv")

  if (verbose) message("Schur-CHOLMOD: building XtX and Z'X...")
  XtX <- crossprod(Xsp)
  Z <- Matrix::sparseMatrix(
    i = seq_along(id_rec),
    j = id_rec,
    x = 1,
    dims = c(nrec, N)
  )
  ZtX <- crossprod(Z, Xsp)

  if (verbose) {
    message(sprintf(
      "Schur-CHOLMOD: solving W = Cuu^{-1}Z'X; W dim = %d x %d",
      N, p
    ))
  }
  W <- as.matrix(.cd_solve_cholesky(fac_Cuu, as.matrix(ZtX), "W = Cuu^{-1}Z'X"))

  if (verbose) message("Schur-CHOLMOD: building and factorizing fixed-effect Schur complement...")
  S <- as.matrix(XtX) - as.matrix(crossprod(ZtX, W))
  S <- 0.5 * (S + t(S))
  chol_S <- tryCatch(
    chol(S),
    error = function(e) {
      stop(sprintf(
        "Schur-CHOLMOD: factorization of fixed-effect Schur complement failed. Check collinearity in fixed effects: %s",
        conditionMessage(e)
      ), call. = FALSE)
    }
  )

  idx <- vector("list", U)
  target <- as.logical(target)
  for (k in seq_len(U)) {
    idx[[k]] <- which(mu_animal == k & target)
  }
  nk <- lengths(idx)

  G_num <- matrix(0, U, U)
  G_den <- matrix(0, U, U)

  for (j0 in seq(1L, U, by = block_size)) {
    j1 <- min(j0 + block_size - 1L, U)
    cols <- j0:j1
    bs <- length(cols)

    if (verbose) {
      message(sprintf("Schur-CHOLMOD: processing MU block %d-%d / %d", j0, j1, U))
    }

    B_blk <- matrix(0, nrow = N, ncol = bs)
    for (kk in seq_along(cols)) {
      B_blk[idx[[cols[kk]]], kk] <- 1
    }

    R_blk <- as.matrix(.cd_solve_cholesky(fac_Cuu, B_blk, "Cuu block"))
    rhs_b <- -as.matrix(crossprod(ZtX, R_blk))
    b_blk <- backsolve(chol_S, forwardsolve(t(chol_S), rhs_b))
    U_blk <- R_blk - W %*% b_blk
    YK_blk <- as.matrix(.cd_solve_cholesky(fac_Kinv, B_blk, "Kinv block"))

    for (i in seq_len(U)) {
      ii <- idx[[i]]
      if (!length(ii)) next
      for (kk in seq_along(cols)) {
        j <- cols[kk]
        G_num[i, j] <- sum(U_blk[ii, kk])
        G_den[i, j] <- sum(YK_blk[ii, kk])
      }
    }

    rm(B_blk, R_blk, rhs_b, b_blk, U_blk, YK_blk)
    gc(FALSE)
  }

  G_num <- 0.5 * (G_num + t(G_num))
  G_den <- 0.5 * (G_den + t(G_den))

  CD <- PEVD <- qK <- qC <- matrix(NA_real_, U, U)

  for (i in seq_len(U - 1L)) {
    ni <- nk[i]
    if (ni <= 0L) next
    for (j in seq.int(i + 1L, U)) {
      nj <- nk[j]
      if (nj <= 0L) next

      qK_ij <- G_den[i, i] / (ni * ni) +
        G_den[j, j] / (nj * nj) -
        2 * G_den[i, j] / (ni * nj)

      qC_ij <- G_num[i, i] / (ni * ni) +
        G_num[j, j] / (nj * nj) -
        2 * G_num[i, j] / (ni * nj)

      qK[i, j] <- qK[j, i] <- qK_ij
      qC[i, j] <- qC[j, i] <- qC_ij

      if (qK_ij > 0 && qC_ij >= 0) {
        CD[i, j] <- CD[j, i] <- 1 - lambda * (qC_ij / qK_ij)
        PEVD[i, j] <- PEVD[j, i] <- sigma2e * qC_ij
      }
    }
  }

  dn <- list(mu_names, mu_names)
  dimnames(CD) <- dimnames(PEVD) <- dimnames(qK) <- dimnames(qC) <- dn
  names(nk) <- mu_names

  list(
    CD = CD,
    PEVD = PEVD,
    qK = qK,
    qC = qC,
    n_target_by_MU = nk
  )
}
