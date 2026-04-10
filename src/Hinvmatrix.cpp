#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <unordered_set>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// ===========================================================================
// Helpers para tuning de G (estilo blupf90)
// ===========================================================================

static inline double mean_diag_cpp(const MatrixXd& M) {
  return M.diagonal().mean();
}

static inline double mean_all_cpp(const MatrixXd& M) {
  return M.sum() / static_cast<double>(M.rows() * M.cols());
}

static inline double mean_offdiag_cpp(const MatrixXd& M) {
  const int n = M.rows();
  if (n < 2) return 0.0;
  const double sdiag = M.diagonal().sum();
  return (M.sum() - sdiag) / static_cast<double>(n * n - n);
}

static inline MatrixXd tune_G_affine_cpp(const MatrixXd& G,
                                         int tunedG,
                                         const MatrixXd* A22ptr,
                                         double& a,
                                         double& b) {
  if (tunedG == 0) {
    a = 0.0;
    b = 1.0;
    return G;
  }

  const double gd = mean_diag_cpp(G);
  const double go = mean_offdiag_cpp(G);
  const double denom = gd - go;

  if ((tunedG == 1 || tunedG == 2) && std::abs(denom) < 1e-14)
    stop("Cannot tune G: mean(diag(G)) and mean(offdiag(G)) are too similar.");

  MatrixXd Gt = G;

  if (tunedG == 1) {
    b = 1.0 / denom;
    a = -b * go;
    Gt.array() = a + b * G.array();

  } else if (tunedG == 2) {
    if (A22ptr == nullptr)
      stop("A22 must be supplied when tunedG = 2.");

    const MatrixXd& A22 = *A22ptr;
    if (A22.rows() != G.rows() || A22.cols() != G.cols())
      stop("A22 must have the same dimensions as G when tunedG = 2.");

    const double ad = mean_diag_cpp(A22);
    const double ao = mean_offdiag_cpp(A22);

    b = (ad - ao) / denom;
    a = ao - b * go;
    Gt.array() = a + b * G.array();

  } else if (tunedG == 3) {
    if (A22ptr == nullptr)
      stop("A22 must be supplied when tunedG = 3.");

    const MatrixXd& A22 = *A22ptr;
    if (A22.rows() != G.rows() || A22.cols() != G.cols())
      stop("A22 must have the same dimensions as G when tunedG = 3.");

    b = 1.0;
    a = mean_all_cpp(A22) - mean_all_cpp(G);
    Gt.array() = a + G.array();

  } else {
    stop("Unsupported tunedG value. Use 0, 1, 2 or 3.");
  }

  Gt = (Gt + Gt.transpose()) * 0.5;
  return Gt;
}

// ===========================================================================
// 1. compute_F_ML92
// ===========================================================================

//' Compute inbreeding coefficients via an ML92-style traversal
 //'
 //' Computes inbreeding coefficients using an ancestor-traversal algorithm
 //' inspired by the Meuwissen & Luo (1992) strategy and adapted from the
 //' implementation logic used in MCMCglmm.
 //'
 //' Animals must be numbered 1..N in chronological order (parents before
 //' offspring). Unknown parents coded as 0.
 //'
 //' @param sire Integer vector of sire indices (0 = unknown), length N.
 //' @param dam  Integer vector of dam indices (0 = unknown), length N.
 //' @return Numeric vector of inbreeding coefficients F, length N.
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 NumericVector compute_F_ML92(const IntegerVector& sire,
                              const IntegerVector& dam) {

   if (sire.size() != dam.size())
     stop("sire and dam must have the same length.");
   const int n = sire.size();
   if (n < 1)
     stop("Pedigree is empty.");

   for (int i = 0; i < n; ++i) {
     int s = sire[i], d = dam[i];
     if (s < 0 || s > n || d < 0 || d > n)
       stop("Parent IDs must be in [0, N]. Out-of-range at i = %d.", i + 1);
     int id = i + 1;
     if (s > id || d > id)
       stop("Pedigree not ordered: parent ID > offspring ID at animal %d.", id);
   }

   std::vector<int> S(n + 1), D(n + 1);
   for (int i = 1; i <= n; ++i) {
     S[i] = sire[i - 1];
     D[i] = dam[i - 1];
   }

   std::vector<double> Fint(n + 1, 0.0);
   Fint[0] = -1.0;

   std::vector<double> dii(n + 1, 1.0);
   std::vector<double> li(n + 1, 0.0);
   std::vector<int> AN(2 * n + 10, -1);

   for (int i = 1; i <= n; ++i) {
     dii[i] = 0.5 - 0.25 * (Fint[S[i]] + Fint[D[i]]);
     if (dii[i] <= 0.0)
       stop("Non-positive d_ii at animal %d while computing F (d_ii = %g).", i, dii[i]);

     std::fill(li.begin(), li.end(), 0.0);
     std::fill(AN.begin(), AN.end(), -1);

     li[i] = 1.0;
     double aii = 0.0;

     int j = i;
     int cnt = 0;

     while (j >= 1) {
       const int sj = S[j];
       const int dj = D[j];

       if (sj > 0) {
         if (cnt >= static_cast<int>(AN.size()))
           stop("Internal buffer overflow in compute_F_ML92. Increase AN size.");
         AN[cnt] = sj;
         li[sj] += 0.5 * li[j];
         ++cnt;
       }

       if (dj > 0) {
         if (cnt >= static_cast<int>(AN.size()))
           stop("Internal buffer overflow in compute_F_ML92. Increase AN size.");
         AN[cnt] = dj;
         li[dj] += 0.5 * li[j];
         ++cnt;
       }

       aii += li[j] * li[j] * dii[j];

       int nextj = -1;
       for (int k = 0; k < cnt; ++k) {
         if (AN[k] > nextj)
           nextj = AN[k];
       }

       if (nextj >= 1) {
         for (int k = 0; k < cnt; ++k) {
           if (AN[k] == nextj)
             AN[k] = -1;
         }
       }

       j = nextj;
     }

     Fint[i] = aii - 1.0;
   }

   NumericVector Fout(n);
   for (int i = 1; i <= n; ++i)
     Fout[i - 1] = Fint[i];

   return Fout;
 }

// ===========================================================================
// 2. build_Ainv_sparse_RA
// ===========================================================================

//' Build sparse A-inverse from a renumbered pedigree
 //'
 //' Implements the Henderson (1976) direct method. Animals must be numbered
 //' 1..N in chronological order (parents before offspring). Unknown parents
 //' coded as 0.
 //'
 //' Inbreeding coefficients are computed via \code{compute_F_ML92()}.
 //'
 //' @param sire Integer vector of sire indices (0 = unknown), length N.
 //' @param dam  Integer vector of dam  indices (0 = unknown), length N.
 //' @return List with \code{Ainv} (dgCMatrix, N x N) and \code{F}
 //'   (numeric vector of inbreeding coefficients, length N).
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 List build_Ainv_sparse_RA(const IntegerVector& sire,
                           const IntegerVector& dam) {

   if (sire.size() != dam.size())
     stop("sire and dam must have the same length.");
   const int n = sire.size();
   if (n < 1)
     stop("Pedigree is empty.");

   for (int i = 0; i < n; ++i) {
     int s = sire[i], d = dam[i];
     if (s < 0 || s > n || d < 0 || d > n)
       stop("Parent IDs must be in [0, N]. Out-of-range at i = %d.", i + 1);
     int id = i + 1;
     if (s > id || d > id)
       stop("Pedigree not ordered: parent ID > offspring ID at animal %d.", id);
   }

   std::vector<int> S(n + 1), D(n + 1);
   for (int i = 1; i <= n; ++i) {
     S[i] = sire[i - 1];
     D[i] = dam[i - 1];
   }

   NumericVector Fout = compute_F_ML92(sire, dam);

   std::vector<double> F(n + 1, 0.0);
   F[0] = -1.0;
   for (int i = 1; i <= n; ++i)
     F[i] = Fout[i - 1];

   std::vector<Triplet<double>> tri;
   tri.reserve(static_cast<size_t>(9) * n);

   auto add = [&](int r, int c, double v) {
     tri.emplace_back(r - 1, c - 1, v);
   };

   for (int i = 1; i <= n; ++i) {
     int s = S[i], d = D[i];

     double dii = 0.5 - 0.25 * (F[s] + F[d]);
     if (dii <= 0.0)
       stop("Non-positive d_ii at animal %d (d_ii = %g). Check pedigree.", i, dii);

     double w = 1.0 / dii;

     add(i, i, w);
     if (s) {
       add(i, s, -0.5 * w);
       add(s, i, -0.5 * w);
       add(s, s,  0.25 * w);
     }
     if (d) {
       add(i, d, -0.5 * w);
       add(d, i, -0.5 * w);
       add(d, d,  0.25 * w);
     }
     if (s && d) {
       add(s, d, 0.25 * w);
       add(d, s, 0.25 * w);
     }
   }

   SparseMatrix<double> Ainv(n, n);
   Ainv.setFromTriplets(tri.begin(), tri.end());
   Ainv.makeCompressed();

   return List::create(
     _["Ainv"] = Ainv,
     _["F"]    = Fout
   );
 }

// ===========================================================================
// 3. build_A22 internals
// ===========================================================================

static MatrixXd build_A22_eigen(const IntegerVector& sire,
                                const IntegerVector& dam,
                                const IntegerVector& genotyped_idx,
                                const NumericVector& F) {

  const int N     = sire.size();
  const int n_gen = genotyped_idx.size();

  if (dam.size() != N)
    stop("sire and dam must have the same length N.");
  if (F.size() != N)
    stop("F must have length N (one inbreeding coefficient per animal).");
  if (n_gen == 0)
    stop("genotyped_idx is empty.");

  std::unordered_set<int> seen;
  seen.reserve(static_cast<size_t>(n_gen) * 2);

  for (int i = 0; i < n_gen; ++i) {
    int idx = genotyped_idx[i];
    if (idx < 1 || idx > N)
      stop("genotyped_idx[%d] = %d is out of range [1, N=%d].", i + 1, idx, N);
    if (!seen.insert(idx).second)
      stop("Duplicate entry in genotyped_idx at position %d: %d.", i + 1, idx);
  }

  std::vector<int> S(N), D(N);
  for (int i = 0; i < N; ++i) {
    S[i] = sire[i] - 1;
    D[i] = dam[i]  - 1;
  }

  std::vector<int> gidx(n_gen);
  for (int i = 0; i < n_gen; ++i) gidx[i] = genotyped_idx[i] - 1;

  std::vector<double> Dvec(N);
  for (int i = 0; i < N; ++i) {
    int s = S[i], d = D[i];
    if      (s >= 0 && d >= 0) Dvec[i] = 0.5  - 0.25 * (F[s] + F[d]);
    else if (s >= 0)           Dvec[i] = 0.75 - 0.25 * F[s];
    else if (d >= 0)           Dvec[i] = 0.75 - 0.25 * F[d];
    else                       Dvec[i] = 1.0;

    if (Dvec[i] <= 0.0)
      stop("Non-positive D value at animal %d while building A22.", i + 1);
  }

  MatrixXd A22 = MatrixXd::Zero(n_gen, n_gen);
  std::vector<double> q(N), u(N), a(N);

  Rcout << "Building A22 (" << n_gen << " x " << n_gen
        << ") via Colleau algorithm (N = " << N << ")..." << std::endl;

  for (int jj = 0; jj < n_gen; ++jj) {
    const int j = gidx[jj];

    std::fill(q.begin(), q.end(), 0.0);
    q[j] = 1.0;

    for (int i = j; i >= 0; --i) {
      if (q[i] == 0.0) continue;
      if (S[i] >= 0) q[S[i]] += 0.5 * q[i];
      if (D[i] >= 0) q[D[i]] += 0.5 * q[i];
    }

    for (int i = 0; i < N; ++i)
      u[i] = Dvec[i] * q[i];

    std::fill(a.begin(), a.end(), 0.0);
    for (int i = 0; i < N; ++i) {
      double ps = (S[i] >= 0) ? a[S[i]] : 0.0;
      double pd = (D[i] >= 0) ? a[D[i]] : 0.0;
      a[i] = u[i] + 0.5 * (ps + pd);
    }

    for (int ii = 0; ii < n_gen; ++ii)
      A22(ii, jj) = a[gidx[ii]];

    if ((jj + 1) % 500 == 0 || jj == n_gen - 1)
      Rcout << "  " << (jj + 1) << " / " << n_gen << "\r" << std::flush;
  }
  Rcout << std::endl;

  A22 = (A22 + A22.transpose()) * 0.5;
  return A22;
}

//' Build the dense A22 submatrix for genotyped animals
 //'
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 NumericMatrix build_A22(const IntegerVector& sire,
                         const IntegerVector& dam,
                         const IntegerVector& genotyped_idx,
                         const NumericVector& F) {
   MatrixXd A22 = build_A22_eigen(sire, dam, genotyped_idx, F);

   const int n_gen = A22.rows();
   NumericMatrix out(n_gen, n_gen);
   for (int i = 0; i < n_gen; ++i)
     for (int j = 0; j < n_gen; ++j)
       out(i, j) = A22(i, j);

   return out;
 }

// ===========================================================================
// 4. compute_Ginv internals
// ===========================================================================

struct GinvResultCpp {
  MatrixXd Ginv;
  int n_snps_used;
  int n_snps_removed;
  NumericVector allele_freqs;
  int tunedG;
  double tune_a;
  double tune_b;
  double mean_diag_before;
  double mean_offdiag_before;
  double mean_diag_after;
  double mean_offdiag_after;
};

static GinvResultCpp compute_Ginv_cpp(const MatrixXi& X,
                                      double maf_threshold,
                                      int    missing_code,
                                      double blend,
                                      int    chunk_size,
                                      int    n_threads,
                                      int    tunedG,
                                      const MatrixXd* A22ptr) {

  const int n = X.rows();
  const int m = X.cols();

  if (n < 2) stop("At least 2 animals required.");
  if (m < 2) stop("At least 2 SNPs required.");
  if (maf_threshold < 0.0 || maf_threshold >= 0.5)
    stop("maf_threshold must be in [0, 0.5).");
  if (blend < 0.0 || blend >= 1.0)
    stop("blend must be in [0, 1).");
  if (chunk_size < 1)
    stop("chunk_size must be >= 1.");
  if (n_threads < 1)
    stop("n_threads must be >= 1.");
  if (!(tunedG == 0 || tunedG == 1 || tunedG == 2 || tunedG == 3))
    stop("Unsupported tunedG value. Use 0, 1, 2 or 3.");
  if ((tunedG == 2 || tunedG == 3) && A22ptr == nullptr)
    stop("A22 must be supplied when tunedG = 2 or 3.");

  Eigen::setNbThreads(n_threads);
  Eigen::initParallel();

  Rcout << "Computing allele frequencies (" << m << " SNPs, "
        << n << " animals)..." << std::endl;

  VectorXd p(m);
  std::vector<int> obs_count(m, 0);

#pragma omp parallel for num_threads(n_threads) schedule(static)
  for (int j = 0; j < m; ++j) {
    double sum = 0.0;
    int count  = 0;
    for (int i = 0; i < n; ++i) {
      int g = X(i, j);
      if (g != missing_code) {
        if (g < 0 || g > 2)
          continue;
        sum += g;
        ++count;
      }
    }
    obs_count[j] = count;
    p(j) = (count > 0) ? sum / (2.0 * count)
                       : std::numeric_limits<double>::quiet_NaN();
  }

  Rcout << "Filtering SNPs by MAF >= " << maf_threshold << "..." << std::endl;

  std::vector<unsigned char> keep(m, 0);

#pragma omp parallel for num_threads(n_threads) schedule(static)
  for (int j = 0; j < m; ++j) {
    double f = p(j);
    keep[j] = static_cast<unsigned char>(
      (obs_count[j] > 0) &&
      std::isfinite(f) &&
      (f >= maf_threshold) &&
      (f <= 1.0 - maf_threshold)
    );
  }

  std::vector<int> valid_snps;
  valid_snps.reserve(m);
  for (int j = 0; j < m; ++j)
    if (keep[j]) valid_snps.push_back(j);

    const int m_valid = static_cast<int>(valid_snps.size());
    if (m_valid < 2)
      stop("Fewer than 2 SNPs passed the MAF filter.");

    Rcout << "Retained " << m_valid << " / " << m << " SNPs." << std::endl;

    double denom = 0.0;
#pragma omp parallel for reduction(+:denom) num_threads(n_threads) schedule(static)
    for (int jj = 0; jj < m_valid; ++jj) {
      double f = p(valid_snps[jj]);
      denom += f * (1.0 - f);
    }
    denom *= 2.0;

    if (denom <= 0.0 || !std::isfinite(denom))
      stop("Non-positive or non-finite denominator in G construction. Check SNP filtering.");

    Rcout << "Denominator = " << denom << std::endl;
    Rcout << "Building G matrix (chunks of " << chunk_size
          << " SNPs, using BLAS dsyrk)..." << std::endl;

    const int n_chunks = (m_valid + chunk_size - 1) / chunk_size;
    MatrixXd G = MatrixXd::Zero(n, n);

    MatrixXd Zbuf(n, chunk_size);

    for (int chunk = 0; chunk < n_chunks; ++chunk) {
      const int start    = chunk * chunk_size;
      const int end      = std::min(start + chunk_size, m_valid);
      const int chunk_sz = end - start;

      auto Z = Zbuf.leftCols(chunk_sz);

      for (int j = 0; j < chunk_sz; ++j) {
        const int col = valid_snps[start + j];
        const double c = 2.0 * p(col);

        for (int i = 0; i < n; ++i) {
          int g = X(i, col);
          if (g == missing_code) {
            Z(i, j) = 0.0;
          } else {
            if (g < 0 || g > 2)
              stop("Unexpected genotype code at row %d, col %d: %d", i + 1, col + 1, g);
            Z(i, j) = static_cast<double>(g) - c;
          }
        }
      }

      G.selfadjointView<Lower>().rankUpdate(Z);

      if ((chunk + 1) % 10 == 0 || chunk == n_chunks - 1)
        Rcout << "  chunk " << (chunk + 1) << "/" << n_chunks << "\r" << std::flush;
    }
    Rcout << std::endl;

    G = G.selfadjointView<Lower>();

    Rcout << "Normalizing and blending G..." << std::endl;
    G /= denom;
    G *= (1.0 - blend);
    G.diagonal().array() += blend;
    G = (G + G.transpose()) * 0.5;

    const double mean_diag_before    = mean_diag_cpp(G);
    const double mean_offdiag_before = mean_offdiag_cpp(G);

    double tune_a = 0.0, tune_b = 1.0;
    if (tunedG != 0) {
      Rcout << "Applying G tuning (tunedG = " << tunedG << ")..." << std::endl;
      G = tune_G_affine_cpp(G, tunedG, A22ptr, tune_a, tune_b);
    }

    const double mean_diag_after    = mean_diag_cpp(G);
    const double mean_offdiag_after = mean_offdiag_cpp(G);

    Rcout << "Cholesky factorization (LLT)..." << std::endl;
    LLT<MatrixXd> llt(G);

    if (llt.info() != Success) {
      stop(
        "G is not numerically positive definite and cannot be inverted without regularization. "
        "Set a positive 'blend' value or revise marker filtering/tuning. "
        "No LDLT fallback is used because it may return an unstable inverse."
      );
    }

    Rcout << "Inverting via triangular solve..." << std::endl;
    MatrixXd Linv = llt.matrixL().solve(MatrixXd::Identity(n, n));
    MatrixXd Ginv = Linv.transpose() * Linv;
    Ginv = (Ginv + Ginv.transpose()) * 0.5;
    Rcout << "LLT successful." << std::endl;

    NumericVector freqs(m_valid);
    for (int i = 0; i < m_valid; ++i)
      freqs(i) = p(valid_snps[i]);

    Rcout << "Done (Ginv " << n << " x " << n << ")." << std::endl;

    GinvResultCpp out;
    out.Ginv                = std::move(Ginv);
    out.n_snps_used         = m_valid;
    out.n_snps_removed      = m - m_valid;
    out.allele_freqs        = freqs;
    out.tunedG              = tunedG;
    out.tune_a              = tune_a;
    out.tune_b              = tune_b;
    out.mean_diag_before    = mean_diag_before;
    out.mean_offdiag_before = mean_offdiag_before;
    out.mean_diag_after     = mean_diag_after;
    out.mean_offdiag_after  = mean_offdiag_after;
    return out;
}

//' Compute dense G-inverse for genotyped animals (VanRaden method 1)
 //'
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 List compute_Ginv(const Eigen::MatrixXi& X,
                   double maf_threshold = 0.05,
                   int    missing_code  = 5,
                   double blend         = 0.05,
                   int    chunk_size    = 2000,
                   int    n_threads     = 1,
                   int    tunedG        = 0,
                   Rcpp::Nullable<Rcpp::NumericMatrix> A22 = R_NilValue) {

   MatrixXd A22mat;
   MatrixXd* A22ptr = nullptr;
   if (A22.isNotNull()) {
     NumericMatrix A22r(A22);
     if (A22r.nrow() != X.rows() || A22r.ncol() != X.rows())
       stop("A22 must be n_gen x n_gen with n_gen = number of rows in X.");
     A22mat = as<MatrixXd>(A22r);
     A22mat = (A22mat + A22mat.transpose()) * 0.5;
     A22ptr = &A22mat;
   }

   GinvResultCpp res = compute_Ginv_cpp(X, maf_threshold, missing_code, blend,
                                        chunk_size, n_threads, tunedG, A22ptr);

   return List::create(
     _["Ginv"]                = res.Ginv,
     _["n_snps_used"]         = res.n_snps_used,
     _["n_snps_removed"]      = res.n_snps_removed,
     _["allele_freqs"]        = res.allele_freqs,
     _["tunedG"]              = res.tunedG,
     _["tune_a"]              = res.tune_a,
     _["tune_b"]              = res.tune_b,
     _["mean_diag_before"]    = res.mean_diag_before,
     _["mean_offdiag_before"] = res.mean_offdiag_before,
     _["mean_diag_after"]     = res.mean_diag_after,
     _["mean_offdiag_after"]  = res.mean_offdiag_after
   );
 }

// ===========================================================================
// 5. compute_Hinv
// ===========================================================================

//' Compute sparse H-inverse for a combined pedigree-genomic relationship model
 //'
 //' Implements the generalized form:
 //' \deqn{H^{-1} = A^{-1} + \begin{bmatrix} 0 & 0 \\ 0 & \tau G^{-1} - \omega A_{22}^{-1} \end{bmatrix}}
 //'
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 Eigen::SparseMatrix<double> compute_Hinv(
     const Eigen::SparseMatrix<double>& Ainv,
     const Eigen::MatrixXd&             Ginv,
     const Eigen::MatrixXd&             A22,
     const Rcpp::IntegerVector&         genotyped_idx,
     double tau   = 1.0,
     double omega = 1.0) {

   const int N     = Ainv.rows();
   const int n_gen = genotyped_idx.size();

   if (Ainv.cols() != N)
     stop("Ainv must be square.");
   if (Ginv.rows() != n_gen || Ginv.cols() != n_gen)
     stop("Ginv must be n_gen x n_gen where n_gen = length(genotyped_idx) = %d.", n_gen);
   if (A22.rows() != n_gen || A22.cols() != n_gen)
     stop("A22 must be n_gen x n_gen where n_gen = length(genotyped_idx) = %d.", n_gen);
   if (tau <= 0.0)
     stop("tau must be > 0.");
   if (omega < 0.0)
     stop("omega must be >= 0.");

   std::unordered_set<int> seen;
   seen.reserve(static_cast<size_t>(n_gen) * 2);

   std::vector<int> gidx(n_gen);
   for (int i = 0; i < n_gen; ++i) {
     int idx = genotyped_idx[i];
     if (idx < 1 || idx > N)
       stop("genotyped_idx[%d] = %d is out of range [1, N=%d].", i + 1, idx, N);
     if (!seen.insert(idx).second)
       stop("Duplicate entry in genotyped_idx at position %d: %d.", i + 1, idx);
     gidx[i] = idx - 1;
   }

   MatrixXd A22sym  = (A22  + A22.transpose())  * 0.5;
   MatrixXd Ginvsym = (Ginv + Ginv.transpose()) * 0.5;

   Rcout << "Inverting A22 (" << n_gen << " x " << n_gen << ")..." << std::endl;

   MatrixXd A22inv;
   LLT<MatrixXd> llt_a22(A22sym);
   if (llt_a22.info() == Success) {
     MatrixXd Linv = llt_a22.matrixL().solve(MatrixXd::Identity(n_gen, n_gen));
     A22inv = Linv.transpose() * Linv;
     A22inv = (A22inv + A22inv.transpose()) * 0.5;
     Rcout << "A22 LLT successful." << std::endl;
   } else {
     Rcout << "A22 LLT failed. Trying LDLT..." << std::endl;
     LDLT<MatrixXd> ldlt_a22(A22sym);
     if (ldlt_a22.info() != Success)
       stop("A22 is not positive definite. Check pedigree and genotyped_idx.");
     A22inv = ldlt_a22.solve(MatrixXd::Identity(n_gen, n_gen));
     A22inv = (A22inv + A22inv.transpose()) * 0.5;
     Rcout << "A22 LDLT successful." << std::endl;
   }

   Rcout << "Computing D = tau * Ginv - omega * A22inv..." << std::endl;
   MatrixXd D = tau * Ginvsym - omega * A22inv;
   D = (D + D.transpose()) * 0.5;

   Rcout << "Scattering D onto Ainv (N = " << N
         << ", dense block = " << n_gen << " x " << n_gen << ")..." << std::endl;

   const size_t n_dense = static_cast<size_t>(n_gen) * n_gen;
   std::vector<Triplet<double>> triplets;
   triplets.reserve(static_cast<size_t>(Ainv.nonZeros()) + n_dense);

   for (int k = 0; k < Ainv.outerSize(); ++k)
     for (SparseMatrix<double>::InnerIterator it(Ainv, k); it; ++it)
       triplets.emplace_back(it.row(), it.col(), it.value());

   for (int i = 0; i < n_gen; ++i) {
     for (int j = 0; j < n_gen; ++j) {
       double dval = D(i, j);
       if (dval != 0.0)
         triplets.emplace_back(gidx[i], gidx[j], dval);
     }
   }

   SparseMatrix<double> Hinv(N, N);
   Hinv.setFromTriplets(triplets.begin(), triplets.end());
   Hinv.makeCompressed();

   Rcout << "H-inverse built."
         << "\n  Total nonzeros : " << Hinv.nonZeros()
         << "\n  From Ainv      : " << Ainv.nonZeros()
         << "\n  Dense block    : " << n_dense
         << " (overlap absorbed by setFromTriplets)"
         << std::endl;

   return Hinv;
 }

// ===========================================================================
// 6. compute_Hinv_from_X
// ===========================================================================

//' Compute full H-inverse pipeline from pedigree and genotypes
 //'
 //' Convenience wrapper that performs the full workflow.
 //'
 //' @param sire Integer vector of sire indices (0 = unknown), length N.
 //' @param dam Integer vector of dam indices (0 = unknown), length N.
 //' @param genotyped_idx Integer vector (1-based pedigree indices) for the
 //'   genotyped animals, in the same order as the rows of \code{X}.
 //' @param X Genotype matrix (n_gen x m), coded 0/1/2.
 //' @param maf_threshold Minor allele frequency threshold.
 //' @param missing_code Integer value indicating missing genotypes.
 //' @param blend Blending factor applied to G before optional tuning.
 //' @param chunk_size Number of SNP columns per chunk in G construction.
 //' @param n_threads Number of OpenMP threads.
 //' @param tunedG Integer tuning option for G.
 //' @param tau Scaling factor multiplying \eqn{G^{-1}} in H-inverse.
 //' @param omega Scaling factor multiplying \eqn{A_{22}^{-1}} in H-inverse.
 //' @param return_Ainv Return Ainv in output list. Default TRUE.
 //' @param return_F Return F in output list. Default TRUE.
 //' @param return_A22 Return A22 in output list. Default FALSE.
 //' @param return_Ginv Return Ginv in output list. Default FALSE.
 //' @param return_allele_freqs Return allele frequencies in output list. Default FALSE.
 //' @return List containing Hinv and selected optional objects.
 //' @keywords internal
 //' @noRd
 // [[Rcpp::export]]
 List compute_Hinv_from_X(const Rcpp::IntegerVector& sire,
                          const Rcpp::IntegerVector& dam,
                          const Rcpp::IntegerVector& genotyped_idx,
                          const Eigen::MatrixXi&    X,
                          double maf_threshold      = 0.05,
                          int    missing_code       = 5,
                          double blend              = 0.05,
                          int    chunk_size         = 2000,
                          int    n_threads          = 1,
                          int    tunedG             = 0,
                          double tau                = 1.0,
                          double omega              = 1.0,
                          bool   return_Ainv        = false,
                          bool   return_F           = false,
                          bool   return_A22         = false,
                          bool   return_Ginv        = false,
                          bool   return_allele_freqs= false) {

   const int N     = sire.size();
   const int n_gen = genotyped_idx.size();

   if (dam.size() != N)
     stop("sire and dam must have the same length.");
   if (X.rows() != n_gen)
     stop("Number of rows in X (%d) must equal length(genotyped_idx) (%d).", X.rows(), n_gen);
   if (n_gen < 1)
     stop("genotyped_idx is empty.");
   if (X.cols() < 2)
     stop("X must have at least 2 SNP columns.");

   Rcout << "Step 1/5: Building Ainv and computing F..." << std::endl;
   List Ares = build_Ainv_sparse_RA(sire, dam);
   SparseMatrix<double> Ainv = as<SparseMatrix<double>>(Ares["Ainv"]);
   NumericVector F           = Ares["F"];

   Rcout << "Step 2/5: Building A22..." << std::endl;
   MatrixXd A22 = build_A22_eigen(sire, dam, genotyped_idx, F);

   Rcout << "Step 3/5: Computing Ginv..." << std::endl;
   GinvResultCpp Gres = compute_Ginv_cpp(X, maf_threshold, missing_code, blend,
                                         chunk_size, n_threads, tunedG, &A22);

   Rcout << "Step 4/5: Computing Hinv..." << std::endl;
   SparseMatrix<double> Hinv = compute_Hinv(Ainv, Gres.Ginv, A22, genotyped_idx,
                                            tau, omega);

   Rcout << "Step 5/5: Assembling output..." << std::endl;
   List out;
   out["Hinv"] = Hinv;

   if (return_Ainv)         out["Ainv"] = Ainv;
   if (return_F)            out["F"] = F;
   if (return_A22)          out["A22"] = A22;
   if (return_Ginv)         out["Ginv"] = Gres.Ginv;
   if (return_allele_freqs) out["allele_freqs"] = Gres.allele_freqs;

   out["n_snps_used"]         = Gres.n_snps_used;
   out["n_snps_removed"]      = Gres.n_snps_removed;
   out["tunedG"]              = Gres.tunedG;
   out["tune_a"]              = Gres.tune_a;
   out["tune_b"]              = Gres.tune_b;
   out["mean_diag_before"]    = Gres.mean_diag_before;
   out["mean_offdiag_before"] = Gres.mean_offdiag_before;
   out["mean_diag_after"]     = Gres.mean_diag_after;
   out["mean_offdiag_after"]  = Gres.mean_offdiag_after;
   out["tau"]                 = tau;
   out["omega"]               = omega;

   return out;
 }
