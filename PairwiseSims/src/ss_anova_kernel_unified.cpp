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
// 1. Unified Full Model (For Global Lengthscale Estimation)
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_full_model(const arma::mat& Z_cov, double log_rho, double log_lam, int n) {
  double rho = std::exp(log_rho);
  double lam = std::exp(log_lam);
  
  int P = Z_cov.n_cols;
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n), S_sq = arma::zeros<arma::mat>(n, n);
  arma::mat dS_1D = arma::zeros<arma::mat>(n, n), K_dK_sq = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    
    arma::mat K_raw = arma::exp(-D2 / rho);
    arma::mat dK_raw = K_raw % (D2 / rho);
    
    arma::mat K_j = double_center(K_raw, n);
    arma::mat dK_j = double_center(dK_raw, n);
    
    S_1D += K_j;
    S_sq += K_j % K_j;
    dS_1D += dK_j;
    K_dK_sq += K_j % dK_j;
  }
  
  arma::mat S_2D = 0.5 * (S_1D % S_1D - S_sq);
  arma::mat dS_2D = S_1D % dS_1D - K_dK_sq;
  
  arma::mat K_full = (S_1D + S_2D) / L_total;
  arma::mat dK_full = (dS_1D + dS_2D) / L_total;
  
  arma::mat W = lam * K_full + arma::eye<arma::mat>(n, n);
  arma::mat dW_dlogrho = lam * dK_full;
  arma::mat dW_dloglam = lam * K_full;
  
  return List::create(
    Named("W") = W,
    Named("dW_dlogrho") = dW_dlogrho,
    Named("dW_dloglam") = dW_dloglam
  );
}

// ---------------------------------------------------------
// 2. Precomputations for Variable Importance Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_var_imp_fixed_rho(const arma::mat& Z_cov, IntegerVector test_indices, 
                                    double log_rho0, double log_rho1, 
                                    double log_lam, int n) {
  double rho0 = std::exp(log_rho0);
  double rho1 = std::exp(log_rho1);
  double lam = std::exp(log_lam);
  
  int P = Z_cov.n_cols;
  
  std::vector<bool> is_test(P, false);
  for(int i = 0; i < test_indices.size(); i++) {
    if(test_indices[i] >= 0 && test_indices[i] < P) {
      is_test[test_indices[i]] = true;
    }
  }
  
  arma::mat S_1D0 = arma::zeros<arma::mat>(n, n), S_sq0 = arma::zeros<arma::mat>(n, n);
  arma::mat S_1D1 = arma::zeros<arma::mat>(n, n), S_sq1 = arma::zeros<arma::mat>(n, n);
  arma::mat S_keep0 = arma::zeros<arma::mat>(n, n), S_keep_sq0 = arma::zeros<arma::mat>(n, n);
  
  double P_keep = 0.0;
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_j0 = double_center(arma::exp(-D2 / rho0), n);
    
    S_1D0 += K_j0;
    S_sq0 += K_j0 % K_j0;
    
    double rho_mixed = is_test[j] ? rho1 : rho0;
    arma::mat K_jm = double_center(arma::exp(-D2 / rho_mixed), n);
    S_1D1 += K_jm;
    S_sq1 += K_jm % K_jm;
    
    if(!is_test[j]) {
      S_keep0 += K_j0;
      S_keep_sq0 += K_j0 % K_j0;
      P_keep += 1.0;
    }
  }
  
  arma::mat K_null_sum0 = S_keep0 + 0.5 * (S_keep0 % S_keep0 - S_keep_sq0);
  arma::mat S_all1 = S_1D1 + 0.5 * (S_1D1 % S_1D1 - S_sq1);
  
  arma::mat K_test_raw = S_all1 - K_null_sum0; 
  
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_null = P_keep + (P_keep * (P_keep - 1.0)) / 2.0;
  double L_test = L_total - L_null;
  
  arma::mat K_null = arma::zeros<arma::mat>(n, n);
  if (L_null > 0.5) K_null = K_null_sum0 / L_null;
  
  arma::mat K_test = arma::zeros<arma::mat>(n, n);
  if (L_test > 0.5) K_test = K_test_raw / L_test;
  
  arma::mat W0 = lam * K_null + arma::eye<arma::mat>(n, n);
  arma::mat dW0_dloglam = lam * K_null;
  
  return List::create(
    Named("W0") = W0,
    Named("dW0_dloglam") = dW0_dloglam,
    Named("K_test") = K_test
  );
}

// ---------------------------------------------------------
// 3. Precomputations for Single Interaction Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_interaction_fixed_rho(const arma::mat& Z_cov, int idx1, int idx2, 
                                        double log_rho0, double log_rho1, 
                                        double log_lam, int n) {
  double rho0 = std::exp(log_rho0);
  double rho1 = std::exp(log_rho1);
  double lam = std::exp(log_lam);
  
  int P = Z_cov.n_cols;
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_null = L_total - 1.0; 
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n), S_sq = arma::zeros<arma::mat>(n, n);
  arma::mat K_A = arma::zeros<arma::mat>(n, n), K_B = arma::zeros<arma::mat>(n, n);
  arma::mat K_A0 = arma::zeros<arma::mat>(n, n), K_B0 = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_j0 = double_center(arma::exp(-D2 / rho0), n);
    
    if (j == idx1) {
      K_A0 = K_j0;
      K_A = double_center(arma::exp(-D2 / rho1), n);
    }
    if (j == idx2) {
      K_B0 = K_j0;
      K_B = double_center(arma::exp(-D2 / rho1), n);
    }
    
    S_1D += K_j0;
    S_sq += K_j0 % K_j0;
  }
  
  arma::mat S_all = S_1D + 0.5 * (S_1D % S_1D - S_sq);
  arma::mat K_test0 = K_A0 % K_B0;
  
  arma::mat K_test = K_A % K_B;
  arma::mat K_null = arma::zeros<arma::mat>(n, n);
  
  if (L_null > 0.5) {
    K_null = (S_all - K_test0) / L_null;
  }
  
  arma::mat W0 = lam * K_null + arma::eye<arma::mat>(n, n);
  arma::mat dW0_dloglam = lam * K_null;
  
  return List::create(
    Named("W0") = W0,
    Named("dW0_dloglam") = dW0_dloglam,
    Named("K_test") = K_test
  );
}

// ---------------------------------------------------------
// 4. Precomputations for All-Interactions Null Model
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_all_interactions_fixed_rho(const arma::mat& Z_cov, double fixed_log_rho, 
                                             double log_lam, int n) {
  int P = Z_cov.n_cols;
  double lam = std::exp(log_lam);
  double rho = std::exp(fixed_log_rho); 
  
  arma::mat S_1D = arma::zeros<arma::mat>(n, n);
  arma::mat S_sq = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_j = arma::exp(-D2 / rho);
    arma::mat K_jc = double_center(K_j, n);
    
    S_1D += K_jc;
    S_sq += K_jc % K_jc; 
  }
  
  arma::mat K_test = 0.5 * ((S_1D % S_1D) - S_sq);
  
  arma::mat dW0_dloglam = lam * S_1D;
  arma::mat W0 = arma::eye<arma::mat>(n, n) + lam * S_1D;
  
  return List::create(
    Named("W0") = W0,
    Named("dW0_dloglam") = dW0_dloglam,
    Named("K_test") = K_test
  );
}

// ---------------------------------------------------------
// 5. Out of Sample Prediction Logic
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_unified_K_cross_and_diag(const arma::mat& Z_new, const arma::mat& Z_train, double log_rho) {
  double rho = std::exp(log_rho);
  int n_train = Z_train.n_rows;
  int n_new = Z_new.n_rows;
  int P = Z_train.n_cols;
  
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  
  arma::mat S_1D_cross = arma::zeros<arma::mat>(n_new, n_train);
  arma::mat S_sq_cross = arma::zeros<arma::mat>(n_new, n_train);
  
  arma::vec S_1D_diag = arma::zeros<arma::vec>(n_new);
  arma::vec S_sq_diag = arma::zeros<arma::vec>(n_new);
  
  for(int j = 0; j < P; j++) {
    arma::mat D_train = arma::zeros<arma::mat>(n_train, n_train);
    for(int i=0; i<n_train; i++) {
      for(int k=0; k<n_train; k++) {
        double d = Z_train(i, j) - Z_train(k, j);
        D_train(i, k) = d * d;
      }
    }
    arma::mat R_train = arma::exp(-D_train / rho);
    arma::rowvec col_means_train = arma::mean(R_train, 0);
    double grand_mean = arma::mean(arma::mean(R_train));
    
    arma::mat D_cross = arma::zeros<arma::mat>(n_new, n_train);
    for(int i=0; i<n_new; i++) {
      for(int k=0; k<n_train; k++) {
        double d = Z_new(i, j) - Z_train(k, j);
        D_cross(i, k) = d * d;
      }
    }
    arma::mat R_cross = arma::exp(-D_cross / rho);
    arma::colvec row_means_cross = arma::mean(R_cross, 1);
    
    arma::mat K_cross_j = R_cross;
    for(int i=0; i<n_new; i++) {
      for(int k=0; k<n_train; k++) {
        K_cross_j(i, k) = R_cross(i, k) - row_means_cross(i) - col_means_train(k) + grand_mean;
      }
    }
    S_1D_cross += K_cross_j;
    S_sq_cross += K_cross_j % K_cross_j;
    
    arma::vec K_diag_j = arma::zeros<arma::vec>(n_new);
    for(int i=0; i<n_new; i++) {
      K_diag_j(i) = 1.0 - 2.0 * row_means_cross(i) + grand_mean;
    }
    S_1D_diag += K_diag_j;
    S_sq_diag += K_diag_j % K_diag_j;
  }
  
  arma::mat S_2D_cross = 0.5 * (S_1D_cross % S_1D_cross - S_sq_cross);
  arma::mat K_full_cross = (S_1D_cross + S_2D_cross) / L_total;
  
  arma::vec S_2D_diag = 0.5 * (S_1D_diag % S_1D_diag - S_sq_diag);
  arma::vec K_full_diag = (S_1D_diag + S_2D_diag) / L_total;
  
  return List::create(
    Named("K_cross") = K_full_cross,
    Named("K_diag") = K_full_diag
  );
}

// ---------------------------------------------------------
// 6. Separated Var Imp Model (Maintains Signature)
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_separated_var_imp(const arma::mat& Z_cov, IntegerVector test_indices, 
                            double log_rho0, double log_rho1, 
                            double log_lam0, double log_lam1, int n) {
  double rho0 = std::exp(log_rho0);
  double lam0 = std::exp(log_lam0);
  double rho1 = std::exp(log_rho1);
  double lam1 = std::exp(log_lam1);
  
  int P = Z_cov.n_cols;
  std::vector<bool> is_test(P, false);
  for(int i = 0; i < test_indices.size(); i++) {
    if(test_indices[i] >= 0 && test_indices[i] < P) is_test[test_indices[i]] = true;
  }
  
  arma::mat S_keep0 = arma::zeros<arma::mat>(n, n), S_keep_sq0 = arma::zeros<arma::mat>(n, n);
  arma::mat dS_keep0 = arma::zeros<arma::mat>(n, n), K_dK_keep_sq0 = arma::zeros<arma::mat>(n, n);
  double P_keep = 0.0;
  
  for(int j = 0; j < P; j++) {
    if(!is_test[j]) {
      arma::mat D2(n, n, arma::fill::none);
      for(int r = 0; r < n; r++) {
        for(int c = 0; c < n; c++) {
          double d = Z_cov(r, j) - Z_cov(c, j);
          D2(r, c) = d * d;
        }
      }
      arma::mat K_raw = arma::exp(-D2 / rho0);
      arma::mat dK_raw = K_raw % (D2 / rho0);
      
      arma::mat K_j = double_center(K_raw, n);
      arma::mat dK_j = double_center(dK_raw, n);
      
      S_keep0 += K_j;
      S_keep_sq0 += K_j % K_j;
      dS_keep0 += dK_j;
      K_dK_keep_sq0 += K_j % dK_j;
      P_keep += 1.0;
    }
  }
  
  double L_null = P_keep + (P_keep * (P_keep - 1.0)) / 2.0;
  arma::mat K_null = arma::zeros<arma::mat>(n, n);
  arma::mat dK_null_dlogrho0 = arma::zeros<arma::mat>(n, n);
  if (L_null > 0.5) {
    K_null = (S_keep0 + 0.5 * (S_keep0 % S_keep0 - S_keep_sq0)) / L_null;
    dK_null_dlogrho0 = (dS_keep0 + S_keep0 % dS_keep0 - K_dK_keep_sq0) / L_null;
  }
  
  arma::mat S_1D1 = arma::zeros<arma::mat>(n, n), S_sq1 = arma::zeros<arma::mat>(n, n);
  arma::mat dS_1D1 = arma::zeros<arma::mat>(n, n), K_dK_sq1 = arma::zeros<arma::mat>(n, n);
  arma::mat S_keep1 = arma::zeros<arma::mat>(n, n), S_keep_sq1 = arma::zeros<arma::mat>(n, n);
  arma::mat dS_keep1 = arma::zeros<arma::mat>(n, n), K_dK_keep_sq1 = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_raw = arma::exp(-D2 / rho1);
    arma::mat dK_raw = K_raw % (D2 / rho1);
    
    arma::mat K_j = double_center(K_raw, n);
    arma::mat dK_j = double_center(dK_raw, n);
    
    S_1D1 += K_j;
    S_sq1 += K_j % K_j;
    dS_1D1 += dK_j;
    K_dK_sq1 += K_j % dK_j;
    
    if(!is_test[j]) {
      S_keep1 += K_j;
      S_keep_sq1 += K_j % K_j;
      dS_keep1 += dK_j;
      K_dK_keep_sq1 += K_j % dK_j;
    }
  }
  
  arma::mat S_all1 = S_1D1 + 0.5 * (S_1D1 % S_1D1 - S_sq1);
  arma::mat dS_all1 = dS_1D1 + S_1D1 % dS_1D1 - K_dK_sq1;
  
  arma::mat K_null_sum1 = S_keep1 + 0.5 * (S_keep1 % S_keep1 - S_keep_sq1);
  arma::mat dK_null_sum1 = dS_keep1 + S_keep1 % dS_keep1 - K_dK_keep_sq1;
  
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_test = L_total - L_null;
  
  arma::mat K_test = arma::zeros<arma::mat>(n, n);
  arma::mat dK_test_dlogrho1 = arma::zeros<arma::mat>(n, n);
  
  if (L_test > 0.5) {
    K_test = (S_all1 - K_null_sum1) / L_test;
    dK_test_dlogrho1 = (dS_all1 - dK_null_sum1) / L_test;
  }
  
  arma::mat W = lam0 * K_null + lam1 * K_test + arma::eye<arma::mat>(n, n);
  arma::mat dW_dlogrho0 = lam0 * dK_null_dlogrho0;
  arma::mat dW_dlogrho1 = lam1 * dK_test_dlogrho1;
  arma::mat dW_dloglam0 = lam0 * K_null;
  arma::mat dW_dloglam1 = lam1 * K_test;
  
  return List::create(
    Named("W") = W,
    Named("dW_dlogrho0") = dW_dlogrho0,
    Named("dW_dlogrho1") = dW_dlogrho1,
    Named("dW_dloglam0") = dW_dloglam0,
    Named("dW_dloglam1") = dW_dloglam1
  );
}

// ---------------------------------------------------------
// 7. Separated Interaction Model (Maintains Signature)
// ---------------------------------------------------------
// [[Rcpp::export]]
List eval_separated_interaction(const arma::mat& Z_cov, int idx1, int idx2, 
                                double log_rho0, double log_rho1, 
                                double log_lam0, double log_lam1, int n) {
  double rho0 = std::exp(log_rho0);
  double lam0 = std::exp(log_lam0);
  double rho1 = std::exp(log_rho1);
  double lam1 = std::exp(log_lam1);
  
  int P = Z_cov.n_cols;
  double L_total = (double)P + (double)(P * (P - 1.0)) / 2.0;
  double L_null = L_total - 1.0; 
  
  arma::mat S_1D0 = arma::zeros<arma::mat>(n, n), S_sq0 = arma::zeros<arma::mat>(n, n);
  arma::mat dS_1D0 = arma::zeros<arma::mat>(n, n), K_dK_sq0 = arma::zeros<arma::mat>(n, n);
  arma::mat K_A0 = arma::zeros<arma::mat>(n, n), K_B0 = arma::zeros<arma::mat>(n, n);
  arma::mat dK_A0 = arma::zeros<arma::mat>(n, n), dK_B0 = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    arma::mat D2(n, n, arma::fill::none);
    for(int r = 0; r < n; r++) {
      for(int c = 0; c < n; c++) {
        double d = Z_cov(r, j) - Z_cov(c, j);
        D2(r, c) = d * d;
      }
    }
    arma::mat K_raw = arma::exp(-D2 / rho0);
    arma::mat dK_raw = K_raw % (D2 / rho0);
    
    arma::mat K_j = double_center(K_raw, n);
    arma::mat dK_j = double_center(dK_raw, n);
    
    if (j == idx1) { K_A0 = K_j; dK_A0 = dK_j; }
    if (j == idx2) { K_B0 = K_j; dK_B0 = dK_j; }
    
    S_1D0 += K_j;
    S_sq0 += K_j % K_j;
    dS_1D0 += dK_j;
    K_dK_sq0 += K_j % dK_j;
  }
  
  arma::mat S_all0 = S_1D0 + 0.5 * (S_1D0 % S_1D0 - S_sq0);
  arma::mat dS_all0 = dS_1D0 + S_1D0 % dS_1D0 - K_dK_sq0;
  
  arma::mat K_test0 = K_A0 % K_B0;
  arma::mat dK_test0 = dK_A0 % K_B0 + K_A0 % dK_B0;
  
  arma::mat K_null = arma::zeros<arma::mat>(n, n);
  arma::mat dK_null_dlogrho0 = arma::zeros<arma::mat>(n, n);
  if (L_null > 0.5) {
    K_null = (S_all0 - K_test0) / L_null;
    dK_null_dlogrho0 = (dS_all0 - dK_test0) / L_null;
  }
  
  arma::mat K_A1 = arma::zeros<arma::mat>(n, n), K_B1 = arma::zeros<arma::mat>(n, n);
  arma::mat dK_A1 = arma::zeros<arma::mat>(n, n), dK_B1 = arma::zeros<arma::mat>(n, n);
  
  for(int j = 0; j < P; j++) {
    if (j == idx1 || j == idx2) {
      arma::mat D2(n, n, arma::fill::none);
      for(int r = 0; r < n; r++) {
        for(int c = 0; c < n; c++) {
          double d = Z_cov(r, j) - Z_cov(c, j);
          D2(r, c) = d * d;
        }
      }
      arma::mat K_raw = arma::exp(-D2 / rho1);
      arma::mat dK_raw = K_raw % (D2 / rho1);
      
      arma::mat K_j = double_center(K_raw, n);
      arma::mat dK_j = double_center(dK_raw, n);
      
      if (j == idx1) { K_A1 = K_j; dK_A1 = dK_j; }
      if (j == idx2) { K_B1 = K_j; dK_B1 = dK_j; }
    }
  }
  
  arma::mat K_test = K_A1 % K_B1;
  arma::mat dK_test_dlogrho1 = dK_A1 % K_B1 + K_A1 % dK_B1;
  
  arma::mat W = lam0 * K_null + lam1 * K_test + arma::eye<arma::mat>(n, n);
  arma::mat dW_dlogrho0 = lam0 * dK_null_dlogrho0;
  arma::mat dW_dlogrho1 = lam1 * dK_test_dlogrho1;
  arma::mat dW_dloglam0 = lam0 * K_null;
  arma::mat dW_dloglam1 = lam1 * K_test;
  
  return List::create(
    Named("W") = W,
    Named("dW_dlogrho0") = dW_dlogrho0,
    Named("dW_dlogrho1") = dW_dlogrho1,
    Named("dW_dloglam0") = dW_dloglam0,
    Named("dW_dloglam1") = dW_dloglam1
  );
}