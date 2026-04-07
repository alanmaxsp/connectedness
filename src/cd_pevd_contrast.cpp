#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::Triplet;

static inline void add_trip(std::vector<Triplet<double>>& tr, int i, int j, double v) {
  if (v != 0.0) tr.emplace_back(i, j, v);
}

//' Compute CD and PEVD via MME contrast using a generic inverse relationship matrix
//'
//' This function computes connectedness statistics using a generic inverse
//' relationship/covariance matrix for the animal effect, denoted here as
//' \eqn{K^{-1}}. Therefore, it can be used with:
//' \itemize{
//'   \item \code{Ainv} for pedigree-based connectedness,
//'   \item \code{Hinv} for single-step connectedness,
//'   \item or any other compatible inverse kernel matrix.
//' }
//'
//' The model is:
//' \deqn{y = Xb + Zu + e}
//' with
//' \deqn{u \sim N(0, K \sigma_a^2), \quad e \sim N(0, I \sigma_e^2)}
//'
//' so the MME animal block is built as:
//' \deqn{Z'Z + \lambda K^{-1}}, \quad \lambda = \sigma_e^2 / \sigma_a^2
//'
//' @param Kinv Sparse N x N inverse relationship/covariance matrix
//'   (e.g. Ainv or Hinv; dgCMatrix).
//' @param id_rec Integer vector (length nrec) with animal index (1..N) for
//'   each phenotypic record.
//' @param X Sparse model matrix for fixed effects (nrec x p, dgCMatrix).
//' @param mu_animal Integer vector (length N): MU index (1..U) per animal;
//'   0 = not a target animal.
//' @param target_nullable Logical vector (length N) indicating target animals,
//'   or NULL to include all animals with mu_animal > 0.
//' @param sigma2a Additive genetic variance (or more generally, variance of u).
//' @param sigma2e Residual variance.
//' @param mu_names_nullable Optional character vector of MU names (length U).
//' @return List with CD, PEVD, qK, qC matrices (U x U) and n_target_by_MU.
//' @export
// [[Rcpp::export]]
Rcpp::List cd_contrast_mu_mme_sparse(
    const Eigen::SparseMatrix<double>& Kinv,
    const Rcpp::IntegerVector& id_rec,
    const Eigen::SparseMatrix<double>& X,
    const Rcpp::IntegerVector& mu_animal,
    Rcpp::Nullable<Rcpp::LogicalVector> target_nullable,
    const double sigma2a,
    const double sigma2e,
    Rcpp::Nullable<Rcpp::CharacterVector> mu_names_nullable = R_NilValue
) {
  const int N = Kinv.rows();
  if (Kinv.cols() != N) Rcpp::stop("Kinv must be square.");

  const int nrec = id_rec.size();
  if (X.rows() != nrec) Rcpp::stop("X must have nrec rows.");
  if (mu_animal.size() != N) Rcpp::stop("mu_animal must have length N.");
  if (sigma2a <= 0.0) Rcpp::stop("sigma2a must be > 0.");
  if (sigma2e < 0.0)  Rcpp::stop("sigma2e must be >= 0.");

  const double lambda = sigma2e / sigma2a;
  const int p = X.cols();
  const int M = p + N;

  // target animals
  std::vector<char> target(N, 1);
  if (target_nullable.isNotNull()) {
    Rcpp::LogicalVector t = Rcpp::LogicalVector(target_nullable);
    if (t.size() != N) Rcpp::stop("target must have length N.");
    for (int i = 0; i < N; ++i) target[i] = (t[i] == TRUE);
  }

  // Number of MUs and membership
  int U = 0;
  for (int i = 0; i < N; ++i) U = std::max(U, (int)mu_animal[i]);
  if (U < 2) Rcpp::stop("At least 2 MUs required (mu_animal).");

  std::vector< std::vector<int> > idx(U);
  for (int a = 0; a < N; ++a) {
    int mu = mu_animal[a];
    if (mu > 0 && target[a]) idx[mu - 1].push_back(a);
  }

  std::vector<double> nk(U, 0.0);
  for (int k = 0; k < U; ++k) nk[k] = (double)idx[k].size();

  // B (N x U)
  MatrixXd B = MatrixXd::Zero(N, U);
  for (int k = 0; k < U; ++k)
    for (int a : idx[k]) B(a, k) = 1.0;

  // D = diag(Z'Z), where Z is incidence for animal effect
  VectorXd D = VectorXd::Zero(N);
  for (int r = 0; r < nrec; ++r) {
    int a = id_rec[r];
    if (a < 1 || a > N) Rcpp::stop("id_rec out of range 1..N at r=%d.", r + 1);
    D[a - 1] += 1.0;
  }

  // XtX
  SparseMatrix<double> XtX = (X.transpose() * X).pruned();
  XtX.makeCompressed();

  // X'Z
  std::vector< std::unordered_map<int,double> > col_acc(p);
  col_acc.reserve(p);

  for (int j = 0; j < p; ++j) {
    auto& mp = col_acc[j];
    for (SparseMatrix<double>::InnerIterator it(X, j); it; ++it) {
      int r = it.row();
      double val = it.value();
      int a = id_rec[r] - 1;
      mp[a] += val;
    }
  }

  std::vector<Triplet<double>> tr_XtZ;
  tr_XtZ.reserve((size_t)X.nonZeros());
  for (int j = 0; j < p; ++j) {
    for (auto const& kv : col_acc[j])
      add_trip(tr_XtZ, j, kv.first, kv.second);
    col_acc[j].clear();
    col_acc[j].rehash(0);
  }

  SparseMatrix<double> XtZ(p, N);
  XtZ.setFromTriplets(tr_XtZ.begin(), tr_XtZ.end());
  XtZ.makeCompressed();

  SparseMatrix<double> ZtX = XtZ.transpose();

  // Cuu = Z'Z + lambda * Kinv
  SparseMatrix<double> Cuu = Kinv;
  Cuu *= lambda;
  Cuu.makeCompressed();
  for (int i = 0; i < N; ++i) Cuu.coeffRef(i, i) += D[i];
  Cuu.makeCompressed();

  // Full MME
  std::vector<Triplet<double>> tr_M;
  tr_M.reserve((size_t)(XtX.nonZeros() + XtZ.nonZeros() + ZtX.nonZeros() + Cuu.nonZeros()));

  for (int k = 0; k < XtX.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(XtX, k); it; ++it)
      add_trip(tr_M, it.row(), it.col(), it.value());

  for (int k = 0; k < XtZ.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(XtZ, k); it; ++it)
      add_trip(tr_M, it.row(), p + it.col(), it.value());

  for (int k = 0; k < ZtX.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(ZtX, k); it; ++it)
      add_trip(tr_M, p + it.row(), it.col(), it.value());

  for (int k = 0; k < Cuu.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(Cuu, k); it; ++it)
      add_trip(tr_M, p + it.row(), p + it.col(), it.value());

  SparseMatrix<double> MME(M, M);
  MME.setFromTriplets(tr_M.begin(), tr_M.end());
  MME.makeCompressed();

  // Solve MME for [ *, W_K B ]
  MatrixXd RHS = MatrixXd::Zero(M, U);
  RHS.block(p, 0, N, U) = B;

  Eigen::SimplicialLDLT<SparseMatrix<double>> solverMME;
  solverMME.compute(MME);
  if (solverMME.info() != Eigen::Success)
    Rcpp::stop("Sparse factorization of MME failed. Check collinearity in X or definiteness of the system.");

  MatrixXd SOL = solverMME.solve(RHS);
  if (solverMME.info() != Eigen::Success)
    Rcpp::stop("MME solve failed.");

  MatrixXd WKB = SOL.block(p, 0, N, U);

  // Solve K * B indirectly from Kinv * Y = B
  Eigen::SimplicialLDLT<SparseMatrix<double>> solverKinv;
  solverKinv.compute(Kinv);
  if (solverKinv.info() != Eigen::Success)
    Rcpp::stop("Sparse factorization of Kinv failed.");

  MatrixXd Y_K = solverKinv.solve(B);
  if (solverKinv.info() != Eigen::Success)
    Rcpp::stop("Solve with Kinv failed.");

  // Aggregate by MU
  MatrixXd G_den = MatrixXd::Zero(U, U);
  MatrixXd G_num = MatrixXd::Zero(U, U);

  for (int i = 0; i < U; ++i) {
    const auto& ii = idx[i];
    if (ii.empty()) continue;

    for (int j = 0; j < U; ++j) {
      double sden = 0.0, snum = 0.0;
      for (int a : ii) {
        sden += Y_K(a, j);
        snum += WKB(a, j);
      }
      G_den(i, j) = sden;
      G_num(i, j) = snum;
    }
  }

  // CD and PEVD
  Rcpp::NumericMatrix CD(U, U), PEVD(U, U), qK_mat(U, U), qC_mat(U, U);
  for (int i = 0; i < U; ++i)
    for (int j = 0; j < U; ++j)
      CD(i, j) = PEVD(i, j) = qK_mat(i, j) = qC_mat(i, j) = NA_REAL;

  for (int i = 0; i < U - 1; ++i) {
    if (nk[i] <= 0.0) continue;

    for (int j = i + 1; j < U; ++j) {
      if (nk[j] <= 0.0) continue;

      const double ni = nk[i], nj = nk[j];

      const double qK =
        G_den(i, i)/(ni * ni) +
        G_den(j, j)/(nj * nj) -
        2.0 * G_den(i, j)/(ni * nj);

      const double qC =
        G_num(i, i)/(ni * ni) +
        G_num(j, j)/(nj * nj) -
        2.0 * G_num(i, j)/(ni * nj);

      if (!(qK > 0.0) || !(qC >= 0.0)) {
        CD(i, j) = CD(j, i) = NA_REAL;
        continue;
      }

      CD(i, j)   = CD(j, i)   = 1.0 - lambda * (qC / qK);
      PEVD(i, j) = PEVD(j, i) = sigma2e * qC;
      qK_mat(i, j) = qK_mat(j, i) = qK;
      qC_mat(i, j) = qC_mat(j, i) = qC;
    }
  }

  Rcpp::NumericVector nk_out(U);
  for (int k = 0; k < U; ++k) nk_out[k] = nk[k];

  if (mu_names_nullable.isNotNull()) {
    Rcpp::CharacterVector mu_names(mu_names_nullable);
    if (mu_names.size() != U)
      Rcpp::stop("mu_names must have length U (U=%d).", U);

    Rcpp::List dn = Rcpp::List::create(mu_names, mu_names);
    CD.attr("dimnames")   = dn;
    PEVD.attr("dimnames") = dn;
    qK_mat.attr("dimnames") = dn;
    qC_mat.attr("dimnames") = dn;
    nk_out.attr("names") = mu_names;
  }

  return Rcpp::List::create(
    Rcpp::Named("CD")             = CD,
    Rcpp::Named("PEVD")           = PEVD,
    Rcpp::Named("qK")             = qK_mat,
    Rcpp::Named("qC")             = qC_mat,
    Rcpp::Named("n_target_by_MU") = nk_out
  );
}
