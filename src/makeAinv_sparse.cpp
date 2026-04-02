#include <RcppEigen.h>
#include <unordered_map>
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]

// Key para cachear R(a,b) con a<b
static inline long long pair_key(int a, int b) {
  if (a > b) std::swap(a, b);
  return ( (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b) );
}

//' Build sparse A-inverse from a renumbered pedigree
//'
//' Internal C++ function. Use [build_Ainv()] from R instead.
//'
//' @param sire Integer vector of sire indices (0 = unknown), length N.
//' @param dam Integer vector of dam indices (0 = unknown), length N.
//' @param cache_parent_pairs Logical. Cache R(sire,dam) for repeated full-sib
//'   families. Default TRUE.
//' @return List with `Ainv` (dgCMatrix, N x N) and `F` (inbreeding coefficients).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List build_Ainv_sparse_RA(const Rcpp::IntegerVector& sire,
                                const Rcpp::IntegerVector& dam,
                                bool cache_parent_pairs = true) {

  if (sire.size() != dam.size()) Rcpp::stop("sire and dam must have same length");
  const int n = sire.size();

  for (int i = 0; i < n; ++i) {
    int s = sire[i], d = dam[i];
    if (s < 0 || s > n || d < 0 || d > n)
      Rcpp::stop("Parent IDs must be in [0,n]. Out-of-range at i=%d.", i+1);
    int id = i + 1;
    if (s > id || d > id)
      Rcpp::stop("Pedigree not ordered: parent ID > animal ID at animal %d.", id);
  }

  std::vector<int> S(n + 1), D(n + 1);
  for (int i = 1; i <= n; ++i) {
    S[i] = sire[i - 1];
    D[i] = dam[i - 1];
  }

  std::vector<double> F(n + 1, 0.0);
  F[0] = -1.0;

  std::unordered_map<long long, double> Rcache;
  if (cache_parent_pairs) Rcache.reserve(static_cast<size_t>(n) * 2);

  std::function<double(int,int)> R = [&](int a, int b) -> double {
    if (a == 0 || b == 0) return 0.0;
    if (a == b) return 1.0 + F[a];
    if (a < b) {
      int sb = S[b], db = D[b];
      return 0.5 * ( R(a, sb) + R(a, db) );
    } else {
      int sa = S[a], da = D[a];
      return 0.5 * ( R(b, sa) + R(b, da) );
    }
  };

  for (int i = 1; i <= n; ++i) {
    int s = S[i], d = D[i];
    if (s == 0 || d == 0) { F[i] = 0.0; continue; }
    if (cache_parent_pairs) {
      long long k = pair_key(s, d);
      auto it = Rcache.find(k);
      if (it != Rcache.end()) {
        F[i] = 0.5 * it->second;
      } else {
        double Rij = R(s, d);
        Rcache.emplace(k, Rij);
        F[i] = 0.5 * Rij;
      }
    } else {
      F[i] = 0.5 * R(s, d);
    }
  }

  std::vector<Eigen::Triplet<double>> tri;
  tri.reserve(static_cast<size_t>(9) * n);

  auto add = [&](int r, int c, double v) {
    tri.emplace_back(r - 1, c - 1, v);
  };

  for (int i = 1; i <= n; ++i) {
    int s = S[i], d = D[i];
    double dii = 0.5 - 0.25 * (F[s] + F[d]);
    if (dii <= 0.0)
      Rcpp::stop("Non-positive Dii at animal %d (Dii=%g). Check pedigree.", i, dii);
    double w = 1.0 / dii;

    add(i, i, w);
    if (s != 0) { add(i,s,-0.5*w); add(s,i,-0.5*w); add(s,s,0.25*w); }
    if (d != 0) { add(i,d,-0.5*w); add(d,i,-0.5*w); add(d,d,0.25*w); }
    if (s != 0 && d != 0) { add(s,d,0.25*w); add(d,s,0.25*w); }
  }

  Eigen::SparseMatrix<double> Ainv(n, n);
  Ainv.setFromTriplets(tri.begin(), tri.end());
  Ainv.makeCompressed();

  Rcpp::NumericVector Fout(n);
  for (int i = 1; i <= n; ++i) Fout[i - 1] = F[i];

  return Rcpp::List::create(
    Rcpp::_["Ainv"] = Ainv,
    Rcpp::_["F"]    = Fout
  );
}
