#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <unordered_map>

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

// ================================
// Row grouping for "unique points"
// ================================

// Compute row groups: identical rows get same group ID
void compute_row_groups(const IntegerMatrix &obs,
                        std::vector<int> &group_id,
                        std::vector<int> &group_size) {
  int m = obs.nrow();
  int n = obs.ncol();

  group_id.assign(m, -1);
  group_size.clear();

  std::unordered_map<long long, int> key2id;
  key2id.reserve(m * 2);
  long long base = static_cast<long long>(n) + 1LL;

  for (int i = 0; i < m; ++i) {
    long long key = 0;
    for (int j = 0; j < n; ++j) {
      key = key * base + static_cast<long long>(obs(i, j));
    }
    auto it = key2id.find(key);
    int gid;
    if (it == key2id.end()) {
      gid = static_cast<int>(group_size.size());
      key2id[key] = gid;
      group_size.push_back(0);
    } else {
      gid = it->second;
    }
    group_id[i] = gid;
    group_size[gid]++;
  }
}

// =========================================================
// Tie-splitting accumulators
// =========================================================

// 1) Original method: tie-splitting over ALL winner rows
//    acc[winners] += 1/|winners|
void accumulate_winners_all(std::vector<double> &acc,
                            const std::vector<int> &d) {
  int m = static_cast<int>(d.size());
  if (m == 0) return;

  int md = d[0];
  for (int i = 1; i < m; ++i) {
    if (d[i] < md) md = d[i];
  }
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

// 2) New method: tie-splitting over UNIQUE POINTS (groups)
//    - winners = rows with minimal distance
//    - collapse by group_id -> winning groups
//    - each winning group g gets 1 / (#winning groups)
//    - every row in group g gets that same value (no splitting within group)
void accumulate_winners_unique_groups(std::vector<double> &acc_group,
                                      const std::vector<int> &d,
                                      const std::vector<int> &group_id) {
  int m = static_cast<int>(d.size());
  if (m == 0) return;

  int md = d[0];
  for (int i = 1; i < m; ++i) {
    if (d[i] < md) md = d[i];
  }

  int G = static_cast<int>(acc_group.size());
  std::vector<char> is_winning_group(G, 0);
  int uniq_cnt = 0;

  // mark winning groups
  for (int i = 0; i < m; ++i) {
    if (d[i] == md) {
      int gid = group_id[i];
      if (!is_winning_group[gid]) {
        is_winning_group[gid] = 1;
        ++uniq_cnt;
      }
    }
  }
  if (uniq_cnt == 0) return;

  double inc_group = 1.0 / static_cast<double>(uniq_cnt);

  // add same increment to each winning group
  for (int g = 0; g < G; ++g) {
    if (is_winning_group[g]) {
      acc_group[g] += inc_group;
    }
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
  double Wd = static_cast<double>(W);

  // handle m == 1, m == 2 — both measures coincide
  if (m == 1) {
    NumericVector w(1, Wd);
    NumericVector p(1, 1.0);
    return List::create(
      _["weights"]        = w,
      _["weights_unique"] = w,
      _["total"]          = Wd,
      _["probs"]          = p,
      _["probs_unique"]   = p,
      _["mode"]           = std::string("exact-serial")
    );
  }
  if (m == 2) {
    NumericVector w(2, Wd / 2.0);
    NumericVector p(2, 0.5);
    return List::create(
      _["weights"]        = w,
      _["weights_unique"] = w,
      _["total"]          = Wd,
      _["probs"]          = p,
      _["probs_unique"]   = p,
      _["mode"]           = std::string("exact-serial")
    );
  }

  // obs_inv: list of inverse permutations (1..n)
  std::vector< std::vector<int> > obs_inv(m);
  for (int i = 0; i < m; ++i) {
    std::vector<int> row(n);
    for (int j = 0; j < n; ++j) row[j] = obs(i, j);
    obs_inv[i] = invert_perm(row);
  }

  // row groups (unique points)
  std::vector<int> group_id, group_size;
  compute_row_groups(obs, group_id, group_size);
  int G = static_cast<int>(group_size.size());

  std::vector<double> acc_all(m, 0.0);      // row-level (old method)
  std::vector<double> acc_group(G, 0.0);    // group-level (new method)

  // enumerate all permutations of 1..n
  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = i + 1;

  do {
    std::vector<int> d = dist_to_obs(perm, obs_inv);
    accumulate_winners_all(acc_all, d);                         // method 1
    accumulate_winners_unique_groups(acc_group, d, group_id);   // method 2
  } while (std::next_permutation(perm.begin(), perm.end()));

  NumericVector weights(m), weights_unique(m);
  NumericVector probs(m), probs_unique(m);

  for (int i = 0; i < m; ++i) {
    int gid = group_id[i];

    // method 1: row-level
    weights[i] = acc_all[i];
    probs[i]   = acc_all[i] / Wd;

    // method 2: group-level, copied to all duplicates
    weights_unique[i] = acc_group[gid];
    probs_unique[i]   = acc_group[gid] / Wd;
  }

  return List::create(
    _["weights"]        = weights,
    _["weights_unique"] = weights_unique,
    _["total"]          = Wd,
    _["probs"]          = probs,
    _["probs_unique"]   = probs_unique,
    _["mode"]           = std::string("exact-serial")
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
  double total_perm_d = static_cast<double>(total_perm);
  double Kd = static_cast<double>(K);

  if (m == 1) {
    NumericVector w(1, Kd);
    NumericVector p(1, 1.0);
    return List::create(
      _["weights"]          = w,
      _["weights_unique"]   = w,
      _["total"]            = total_perm_d,
      _["total_K"]          = Kd,
      _["weights_K"]        = w,
      _["weights_unique_K"] = w,
      _["probs"]            = p,
      _["probs_unique"]     = p,
      _["mode"]             = std::string("mc-serial")
    );
  }
  if (m == 2) {
    NumericVector w(2, Kd / 2.0);
    NumericVector p(2, 0.5);
    return List::create(
      _["weights"]          = w,
      _["weights_unique"]   = w,
      _["total"]            = total_perm_d,
      _["total_K"]          = Kd,
      _["weights_K"]        = w,
      _["weights_unique_K"] = w,
      _["probs"]            = p,
      _["probs_unique"]     = p,
      _["mode"]             = std::string("mc-serial")
    );
  }

  // obs_inv
  std::vector< std::vector<int> > obs_inv(m);
  for (int i = 0; i < m; ++i) {
    std::vector<int> row(n);
    for (int j = 0; j < n; ++j) row[j] = obs(i, j);
    obs_inv[i] = invert_perm(row);
  }

  // row groups
  std::vector<int> group_id, group_size;
  compute_row_groups(obs, group_id, group_size);
  int G = static_cast<int>(group_size.size());

  std::vector<double> acc_all(m, 0.0);      // rows
  std::vector<double> acc_group(G, 0.0);    // groups

  // RNG
  std::mt19937 rng;
  if (seed > 0) {
    rng.seed(static_cast<unsigned int>(seed));
  } else {
    std::random_device rd;
    rng.seed(rd());
  }

  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 1); // 1..n

  for (int t = 0; t < K; ++t) {
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<int> d = dist_to_obs(perm, obs_inv);
    accumulate_winners_all(acc_all, d);
    accumulate_winners_unique_groups(acc_group, d, group_id);
  }

  NumericVector weights(m), weights_K(m), probs(m);
  NumericVector weights_unique(m), weights_unique_K(m), probs_unique(m);

  for (int i = 0; i < m; ++i) {
    int gid = group_id[i];

    // method 1: row-level splits
    weights_K[i] = acc_all[i];
    probs[i]     = acc_all[i] / Kd;
    weights[i]   = probs[i] * total_perm_d;

    // method 2: group-level, same for all duplicates
    weights_unique_K[i] = acc_group[gid];
    probs_unique[i]     = acc_group[gid] / Kd;
    weights_unique[i]   = probs_unique[i] * total_perm_d;
  }

  return List::create(
    _["weights"]          = weights,
    _["weights_unique"]   = weights_unique,
    _["weights_K"]        = weights_K,
    _["weights_unique_K"] = weights_unique_K,
    _["total"]            = total_perm_d,
    _["total_K"]          = Kd,
    _["probs"]            = probs,
    _["probs_unique"]     = probs_unique,
    _["mode"]             = std::string("mc-serial")
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

  // order is based on the "all-winners" probs (original method)
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
