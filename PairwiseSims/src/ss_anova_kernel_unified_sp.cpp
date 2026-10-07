#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// O(n^2) Double-Centering: Replaces C * M * C
inline arma::mat double_center(const arma::mat& M, int n) {
  arma::rowvec col_means = arma::mean(M, 0);
  arma::colvec row_means = arma::mean(M, 1);
  double grand_mean = arma::mean(col_means);
  
  arma::mat centered = M;
  centered.each_row() -= col_means;
  centered.each_col() -= row_means;
  centered += grand_mean;
  return centered;
}

// ---------------------------------------------------------
// 1. Unified Full Model (4 Parameters for Optimization)
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_full_model_sp(const arma::mat& Z_cov, const arma::mat& coords_sp, 
                                double log_rho_a, double log_lam_a, 
                                double log_rho_sp, double log_lam_sp, int n) {
  
  double rho_a = std::exp(log_rho_a);
  double lam_a = std::exp(log_lam_a);
  double rho_sp = std::exp(log_rho_sp);
  double lam_sp = std::exp(log_lam_sp);
  
  int P = Z_cov.n_cols;
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n), S_sq = arma::zeros<arma::mat>(n, n);
  arma::mat dS_1D = arma::zeros<arma::mat>(n, n), K_dK_sq = arma::zeros<arma::mat>(n, n);
  
  // SS-ANOVA Space: Compute distances on the fly
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    
    arma::mat K_raw = arma::exp(-D2 / rho_a);
    arma::mat dK_raw = K_raw % (D2 / rho_a);
    
    arma::mat K_j = double_center(K_raw, n);
    arma::mat dK_j = double_center(dK_raw, n);
    
    S_1D += K_j;
    S_sq += K_j % K_j;
    dS_1D += dK_j;
    K_dK_sq += K_j % dK_j;
  }
  
  arma::mat S_2D = 0.5 * (S_1D % S_1D - S_sq);
  arma::mat dS_2D = S_1D % dS_1D - K_dK_sq;
  
  arma::mat K_anova = (S_1D + S_2D) / L_total;
  arma::mat dK_anova = (dS_1D + dS_2D) / L_total;
  
  // Spatial Field: Compute distances on the fly
  arma::mat D2_sp(n, n, arma::fill::zeros);
  for(int k = 0; k < coords_sp.n_cols; k++) {
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = coords_sp(r, k) - coords_sp(c, k);
        D2_sp(r, c) += d * d;
      }
    }
  }
  
  arma::mat K_sp = arma::exp(-D2_sp / rho_sp);
  arma::mat dK_sp = K_sp % (D2_sp / rho_sp);
  
  arma::mat W = lam_a * K_anova + lam_sp * K_sp + arma::eye<arma::mat>(n, n);
  
  return List::create(
    Named("W") = W,
    Named("dW_dlogrho_a") = lam_a * dK_anova,
    Named("dW_dloglam_a") = lam_a * K_anova,
    Named("dW_dlogrho_sp") = lam_sp * dK_sp,
    Named("dW_dloglam_sp") = lam_sp * K_sp
  );
}

// ---------------------------------------------------------
// 2. Precomputations for Variable Importance Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_var_imp_fixed_rho_sp(const arma::mat& Z_cov, const arma::mat& coords_sp, 
                                       IntegerVector test_indices, 
                                       double log_rho_a, double log_rho_sp, int n) {
  double rho_a = std::exp(log_rho_a);
  double rho_sp = std::exp(log_rho_sp);
  
  int P = Z_cov.n_cols;
  std::vector<bool> is_test(P, false);
  for(int i = 0; i < test_indices.size(); i++) {
    if(test_indices[i] >= 0 && test_indices[i] < P) is_test[test_indices[i]] = true;
  }
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n), S_sq = arma::zeros<arma::mat>(n, n);
  arma::mat S_keep = arma::zeros<arma::mat>(n, n), S_keep_sq = arma::zeros<arma::mat>(n, n);
  
  double P_keep = 0.0;
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_raw = arma::exp(-D2 / rho_a);
    arma::mat K_j = double_center(K_raw, n);
    
    S_1D += K_j;
    S_sq += K_j % K_j;
    if(!is_test[j]) {
      S_keep += K_j;
      S_keep_sq += K_j % K_j;
      P_keep += 1.0;
    }
  }
  
  arma::mat K_null_a_sum = S_keep + 0.5 * (S_keep % S_keep - S_keep_sq);
  arma::mat S_all = S_1D + 0.5 * (S_1D % S_1D - S_sq);
  
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_null = P_keep + (P_keep * (P_keep - 1.0)) / 2.0;
  double L_test = L_total - L_null;
  
  arma::mat K_null_a = arma::zeros<arma::mat>(n, n);
  if (L_null > 0.5) K_null_a = K_null_a_sum / L_null;
  
  arma::mat K_test = arma::zeros<arma::mat>(n, n);
  if (L_test > 0.5) K_test = (S_all - K_null_a_sum) / L_test;
  
  arma::mat D2_sp(n, n, arma::fill::zeros);
  for(int k = 0; k < coords_sp.n_cols; k++) {
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = coords_sp(r, k) - coords_sp(c, k);
        D2_sp(r, c) += d * d;
      }
    }
  }
  arma::mat K_sp = arma::exp(-D2_sp / rho_sp);
  
  return List::create(
    Named("K_null_a") = K_null_a,
    Named("K_sp") = K_sp,
    Named("K_test") = K_test
  );
}

// ---------------------------------------------------------
// 3. Precomputations for Single Interaction Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_interaction_fixed_rho_sp(const arma::mat& Z_cov, const arma::mat& coords_sp, 
                                           int idx1, int idx2, 
                                           double log_rho_a, double log_rho_sp, int n) {
  double rho_a = std::exp(log_rho_a);
  double rho_sp = std::exp(log_rho_sp);
  
  int P = Z_cov.n_cols;
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_null = L_total - 1.0; 
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n), S_sq = arma::zeros<arma::mat>(n, n);
  arma::mat K_A = arma::zeros<arma::mat>(n, n), K_B = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_raw = arma::exp(-D2 / rho_a);
    arma::mat K_j = double_center(K_raw, n);
    
    if (j == idx1) K_A = K_j;
    if (j == idx2) K_B = K_j;
    
    S_1D += K_j;
    S_sq += K_j % K_j;
  }
  
  arma::mat S_all = S_1D + 0.5 * (S_1D % S_1D - S_sq);
  arma::mat K_test = K_A % K_B;
  
  arma::mat K_null_a = arma::zeros<arma::mat>(n, n);
  if (L_null > 0.5) K_null_a = (S_all - K_test) / L_null;
  
  arma::mat D2_sp(n, n, arma::fill::zeros);
  for(int k = 0; k < coords_sp.n_cols; k++) {
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = coords_sp(r, k) - coords_sp(c, k);
        D2_sp(r, c) += d * d;
      }
    }
  }
  arma::mat K_sp = arma::exp(-D2_sp / rho_sp);
  
  return List::create(
    Named("K_null_a") = K_null_a,
    Named("K_sp") = K_sp,
    Named("K_test") = K_test
  );
}

// ---------------------------------------------------------
// 4. Precomputations for All-Interactions Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_all_interactions_fixed_rho_sp(const arma::mat& Z_cov, const arma::mat& coords_sp, 
                                                double log_rho_a, double log_rho_sp, int n) {
  int P = Z_cov.n_cols;
  double rho_a = std::exp(log_rho_a);
  double rho_sp = std::exp(log_rho_sp);
  
  arma::mat S_1D = arma::zeros(n, n);
  arma::mat S_sq = arma::zeros(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_j = arma::exp(-D2 / rho_a);
    arma::mat K_jc = double_center(K_j, n);
    
    S_1D += K_jc;
    S_sq += K_jc % K_jc;
  }
  
  arma::mat K_null_a = S_1D; 
  arma::mat K_test = 0.5 * ((S_1D % S_1D) - S_sq);
  
  arma::mat D2_sp(n, n, arma::fill::zeros);
  for(int k = 0; k < coords_sp.n_cols; k++) {
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = coords_sp(r, k) - coords_sp(c, k);
        D2_sp(r, c) += d * d;
      }
    }
  }
  arma::mat K_sp = arma::exp(-D2_sp / rho_sp);
  
  return List::create(
    Named("K_null_a") = K_null_a,
    Named("K_sp") = K_sp,
    Named("K_test") = K_test
  );
}