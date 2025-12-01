#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>

using namespace Rcpp;

// ================================
// Fenwick (Binary Indexed) Tree
// ================================

inline void bit_add(std::vector<int> &bit, int idx, int val = 1) {
  int n = static_cast<int>(bit.size()) - 1; // bit[1..n]
  while (idx <= n) {
    bit[idx] += val;
    idx += (idx & -idx);
  }
}

inline int bit_sum(const std::vector<int> &bit, int idx) {
  int s = 0;
  while (idx > 0) {
    s += bit[idx];
    idx -= (idx & -idx);
  }
  return s;
}

// Count inversions of permutation a (values in 1..n, stored in a[1..n])
int kendall_inversions(const std::vector<int> &a) {
  int n = static_cast<int>(a.size()) - 1; // a[1..n]
  std::vector<int> bit(n + 1, 0);
  int inv = 0;
  for (int k = n; k >= 1; --k) {
    int x = a[k];
    if (x > 1) inv += bit_sum(bit, x - 1);
    bit_add(bit, x, 1);
  }
  return inv;
}

// Inverse permutation: p[0..n-1] with values in 1..n -> inv[1..n] = positions (1-based)
std::vector<int> invert_perm(const std::vector<int> &p) {
  int n = static_cast<int>(p.size());
  std::vector<int> inv(n + 1);
  for (int i = 0; i < n; ++i) {
    int val = p[i];       // in 1..n
    inv[val] = i + 1;     // position (1..n)
  }
  return inv;
}

// Distance from a candidate x to a single observation inverse permutation invp
int dist_to_single_obs(const std::vector<int> &x,
                       const std::vector<int> &invp) {
  int n = static_cast<int>(x.size());
  // a[1..n]
  std::vector<int> a(n + 1);
  for (int k = 0; k < n; ++k) {
    int val = x[k];          // 1..n
    a[k + 1] = invp[val];    // position 1..n
  }
  return kendall_inversions(a);
}

// Distances from candidate x to all observations (precomputed inverse perms)
std::vector<int> dist_to_obs(const std::vector<int> &x,
                             const std::vector< std::vector<int> > &obs_inv) {
  int m = static_cast<int>(obs_inv.size());
  std::vector<int> d(m);
  for (int i = 0; i < m; ++i) {
    d[i] = dist_to_single_obs(x, obs_inv[i]);
  }
  return d;
}

// Tie-splitting accumulator increment: acc[winners] += 1/|winners|
void accumulate_winners(std::vector<double> &acc,
                        const std::vector<int> &d) {
  int m = static_cast<int>(d.size());
  int md = d[0];
  for (int i = 1; i < m; ++i) {
    if (d[i] < md) md = d[i];
  }
  // count winners
  int cnt = 0;
  for (int i = 0; i < m; ++i) {
    if (d[i] == md) ++cnt;
  }
  if (cnt == 0) return;
  double inc = 1.0 / static_cast<double>(cnt);
  for (int i = 0; i < m; ++i) {
    if (d[i] == md) acc[i] += inc;
  }
}

// factorial as 64-bit integer (n is small)
long long factorial_int(int n) {
  long long res = 1;
  for (int i = 2; i <= n; ++i) res *= i;
  return res;
}

// =========================================================
// EXACT enumeration (serial) – C++ version of weights_exact_serial()
// =========================================================
List weights_exact_serial_cpp(const IntegerMatrix &obs) {
  int m = obs.nrow();
  int n = obs.ncol();
  long long W = factorial_int(n);

  // handle m == 1, m == 2
  if (m == 1) {
    NumericVector weights(1);
    weights[0] = static_cast<double>(W);
    return List::create(
      _["weights"] = weights,
      _["total"]   = static_cast<double>(W),
      _["probs"]   = NumericVector::create(1.0),
      _["mode"]    = std::string("exact-serial")
    );
  }
  if (m == 2) {
    NumericVector weights(2);
    weights[0] = static_cast<double>(W) / 2.0;
    weights[1] = static_cast<double>(W) / 2.0;
    NumericVector probs(2, 0.5);
    return List::create(
      _["weights"] = weights,
      _["total"]   = static_cast<double>(W),
      _["probs"]   = probs,
      _["mode"]    = std::string("exact-serial")
    );
  }

  // obs_inv: list of inverse permutations (1..n)
  std::vector< std::vector<int> > obs_inv(m);
  for (int i = 0; i < m; ++i) {
    std::vector<int> row(n);
    for (int j = 0; j < n; ++j) row[j] = obs(i, j);
    obs_inv[i] = invert_perm(row);
  }

  // enumerate all permutations of 1..n
  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = i + 1;

  std::vector<double> acc(m, 0.0);

  do {
    std::vector<int> d = dist_to_obs(perm, obs_inv);
    accumulate_winners(acc, d);
  } while (std::next_permutation(perm.begin(), perm.end()));

  NumericVector weights(m);
  NumericVector probs(m);
  for (int i = 0; i < m; ++i) {
    weights[i] = acc[i];
    probs[i]   = acc[i] / static_cast<double>(W);
  }

  return List::create(
    _["weights"] = weights,
    _["total"]   = static_cast<double>(W),
    _["probs"]   = probs,
    _["mode"]    = std::string("exact-serial")
  );
}

// =========================================================
// MONTE CARLO (serial) – C++ version of weights_mc_serial()
// =========================================================
List weights_mc_serial_cpp(const IntegerMatrix &obs,
                           int K = 100000,
                           int seed = -1) {
  int m = obs.nrow();
  int n = obs.ncol();
  long long total_perm = factorial_int(n);

  if (m == 1) {
    NumericVector weights(1);
    weights[0] = static_cast<double>(K);
    return List::create(
      _["weights"] = weights,
      _["total"]   = static_cast<double>(K),
      _["probs"]   = NumericVector::create(1.0),
      _["mode"]    = std::string("mc-serial")
    );
  }
  if (m == 2) {
    NumericVector weights(2);
    weights[0] = static_cast<double>(K) / 2.0;
    weights[1] = static_cast<double>(K) / 2.0;
    NumericVector probs(2, 0.5);
    return List::create(
      _["weights"] = weights,
      _["total"]   = static_cast<double>(K),
      _["probs"]   = probs,
      _["mode"]    = std::string("mc-serial")
    );
  }

  // obs_inv
  std::vector< std::vector<int> > obs_inv(m);
  for (int i = 0; i < m; ++i) {
    std::vector<int> row(n);
    for (int j = 0; j < n; ++j) row[j] = obs(i, j);
    obs_inv[i] = invert_perm(row);
  }

  // RNG
  std::mt19937 rng;
  if (seed > 0) {
    rng.seed(static_cast<unsigned int>(seed));
  } else {
    std::random_device rd;
    rng.seed(rd());
  }

  std::vector<double> acc(m, 0.0);
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 1); // 1..n

  for (int t = 0; t < K; ++t) {
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<int> d = dist_to_obs(perm, obs_inv);
    accumulate_winners(acc, d);
  }

  NumericVector weights(m);
  NumericVector weights_K(m);
  NumericVector probs(m);
  double Kd = static_cast<double>(K);
  double total_perm_d = static_cast<double>(total_perm);

  for (int i = 0; i < m; ++i) {
    weights_K[i] = acc[i];
    probs[i]     = acc[i] / Kd;
    weights[i]   = acc[i] / Kd * total_perm_d;
  }

  return List::create(
    _["weights"]   = weights,
    _["weights_K"] = weights_K,
    _["total"]     = total_perm_d,
    _["total_K"]   = Kd,
    _["probs"]     = probs,
    _["mode"]      = std::string("mc-serial")
  );
}

// =========================================================
// Helper: automatic chooser – C++ version of closest_perm_weights()
// =========================================================
List closest_perm_weights_cpp(const IntegerMatrix &obs,
                              const std::string &method,
                              int K,
                              int workers,          // currently unused (serial only)
                              int seed,
                              long long max_exact_factorial) {
  int n = obs.ncol();

  if (method == "auto") {
    long long W = factorial_int(n);
    if (W <= max_exact_factorial) {
      return weights_exact_serial_cpp(obs);
    } else {
      return weights_mc_serial_cpp(obs, K, seed);
    }
  } else if (method == "exact-serial" || method == "exact-parallel") {
    // "exact-parallel" mapped to the same serial implementation
    return weights_exact_serial_cpp(obs);
  } else if (method == "mc-serial" || method == "mc-parallel") {
    // "mc-parallel" mapped to serial MC here
    return weights_mc_serial_cpp(obs, K, seed);
  } else {
    stop("Unknown method: %s", method);
  }
}

// =========================================================
// Exported entry point: base_measure_order
//   This is what you bind to `_kefV1_base_measure_order`.
// =========================================================

// [[Rcpp::export]]
List base_measure_order_cpp(IntegerMatrix samples,
                            std::string method = "auto",
                            int parallel_workers = 1,
                            int K = 100000,
                            int max_exact_factorial = 3628800,
                            int seed = -1) {

  // compute weights/probs using the C++ closest_perm_weights logic
  List res = closest_perm_weights_cpp(samples, method,
                                      K, parallel_workers,
                                      seed, max_exact_factorial);

  // Optional: also return an order of observations by decreasing probability
  NumericVector probs = res["probs"];
  int m = probs.size();
  IntegerVector ord = seq_len(m);
  // order by decreasing probs
  std::sort(ord.begin(), ord.end(),
            [&probs](int i, int j) {
              // i,j are 1-based indices
              return probs[i - 1] > probs[j - 1];
            });

  res["order"]   = ord;
  res["workers"] = parallel_workers;

  return res;
}

// ---------------------------------------------------------
// If you want *exactly* the symbol `_kefV1_base_measure_order`
// for .Call(), you can add this thin wrapper (and register it).
// ---------------------------------------------------------
extern "C" SEXP _kefV1_base_measure_order(SEXP samplesSEXP,
                                         SEXP methodSEXP,
                                         SEXP workersSEXP) {
  IntegerMatrix samples(samplesSEXP);
  std::string method = as<std::string>(methodSEXP);
  int workers = as<int>(workersSEXP);

  // defaults for K, max_exact_factorial, seed
  int K = 100000;
  int max_exact_factorial = 3628800;
  int seed = -1;

  return base_measure_order_cpp(samples, method,
                                workers,
                                K,
                                max_exact_factorial,
                                seed);
}
