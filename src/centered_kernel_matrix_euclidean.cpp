// src/centered_kernel_matrix.cpp
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat centered_kernel_matrix_euclidean(const arma::vec& eval_points_1,
                                 const arma::vec& eval_points_2,
                                 const arma::vec& centering_grid,
                                 double hurst_coef) {

  int n1 = eval_points_1.n_elem;
  int n2 = eval_points_2.n_elem;
  int ng = centering_grid.n_elem;

  arma::mat h_xx_prime(n1, n2);
  arma::mat h_xz(n1, ng);
  arma::mat h_xprime_z(n2, ng);
  arma::mat h_zz(ng, ng);

  // Term 1: h(x, x')
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      h_xx_prime(i, j) = std::pow(std::abs(eval_points_1(i) - eval_points_2(j)), 2 * hurst_coef);
    }
  }

  // Term 2: h(x, z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < ng; j++) {
      h_xz(i, j) = std::pow(std::abs(eval_points_1(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

  // Term 3: h(x', z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n2; i++) {
    for (int j = 0; j < ng; j++) {
      h_xprime_z(i, j) = std::pow(std::abs(eval_points_2(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

  // Term 4: h(z, z)
#pragma omp parallel for collapse(2)
  for (int i = 0; i < ng; i++) {
    for (int j = 0; j < ng; j++) {
      h_zz(i, j) = std::pow(std::abs(centering_grid(i) - centering_grid(j)), 2 * hurst_coef);
    }
  }

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
