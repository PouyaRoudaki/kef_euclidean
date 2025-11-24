// src/centered_kernel_matrix_order.cpp
#include <RcppArmadillo.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// Compute Kendall-tau "flips" between two permutations given as rows.
// Assumes permutations of 1..d (labels 1,2,...,d).
inline int kendall_flips(const arma::Row<int>& p,
                         const arma::Row<int>& q) {
  int d = p.n_elem;

  // pos_p[x] = position of label x in permutation p
  // pos_q[x] = position of label x in permutation q
  std::vector<int> pos_p(d + 1), pos_q(d + 1);

  for (int i = 0; i < d; ++i) {
    int lp = p(i);   // label at position i in p
    int lq = q(i);   // label at position i in q
    pos_p[lp] = i;
    pos_q[lq] = i;
  }

  int flips = 0;
  // Count discordant pairs (a,b) with 1 <= a < b <= d
  for (int a = 1; a <= d; ++a) {
    for (int b = a + 1; b <= d; ++b) {
      bool order_p = pos_p[a] < pos_p[b];
      bool order_q = pos_q[a] < pos_q[b];
      if (order_p != order_q) {
        ++flips;
      }
    }
  }

  return flips;
}

// [[Rcpp::export]]
arma::mat centered_kernel_matrix_order(const arma::imat& eval_points_1,
                                      const arma::imat& eval_points_2,
                                      const arma::imat& centering_grid,
                                      double hurst_coef) {

  int n1 = eval_points_1.n_rows;
  int n2 = eval_points_2.n_rows;
  int ng = centering_grid.n_rows;
  // int d  = eval_points_1.n_cols; // dimension of permutations (not used below directly)

  arma::mat h_xx_prime(n1, n2);
  arma::mat h_xz(n1, ng);
  arma::mat h_xprime_z(n2, ng);
  arma::mat h_zz(ng, ng);

  // Term 1: h(x, x')
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      int flips = kendall_flips(eval_points_1.row(i), eval_points_2.row(j));
      double d = std::sqrt(static_cast<double>(flips));      // sqrt of Kendall distance
      h_xx_prime(i, j) = std::pow(d, 2.0 * hurst_coef);      // (sqrt d)^(2H) = d^H
    }
  }

  // Term 2: h(x, z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < ng; j++) {
      int flips = kendall_flips(eval_points_1.row(i), centering_grid.row(j));
      double d = std::sqrt(static_cast<double>(flips));
      h_xz(i, j) = std::pow(d, 2.0 * hurst_coef);
    }
  }

  // Term 3: h(x', z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < ng; j++) {
      int flips = kendall_flips(eval_points_2.row(i), centering_grid.row(j));
      double d = std::sqrt(static_cast<double>(flips));
      h_xprime_z(i, j) = std::pow(d, 2.0 * hurst_coef);
    }
  }

  // Term 4: h(z, z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < ng; i++) {
    for (int j = 0; j < ng; j++) {
      int flips = kendall_flips(centering_grid.row(i), centering_grid.row(j));
      double d = std::sqrt(static_cast<double>(flips));
      h_zz(i, j) = std::pow(d, 2.0 * hurst_coef);
    }
  }

  // Centering as in your original code
  arma::vec mean_h_xz = arma::mean(h_xz, 1);              // n1 × 1
  arma::vec mean_h_xprime_z = arma::mean(h_xprime_z, 1);  // n2 × 1
  double mean_h_zz = arma::mean(arma::vectorise(h_zz));

  arma::mat centered_kernel = -0.5 * (
    h_xx_prime -
      arma::repmat(mean_h_xz, 1, n2) -
      arma::repmat(mean_h_xprime_z.t(), n1, 1) +
      mean_h_zz
  );

  return centered_kernel;
}
