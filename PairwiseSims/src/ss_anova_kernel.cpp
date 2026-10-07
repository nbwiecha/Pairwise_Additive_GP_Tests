#include <RcppArmadillo.h>
#include <algorithm> // Required for std::find
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eval_var_effects_H0_multi_fast(List D_1D_list, IntegerVector test_indices, 
                                    double log_rho, double log_lam, int n) {
    double rho = exp(log_rho);
    double lam = exp(log_lam);
    
    int P = D_1D_list.size();
    int P_test = test_indices.size();
    double P_null = (double)(P - P_test);
    
    // Convert Rcpp IntegerVector to std::vector for fast element lookup
    std::vector<int> drop_vec = as<std::vector<int>>(test_indices);

    mat S = zeros<mat>(n, n);
    mat S2 = zeros<mat>(n, n);
    mat dS = zeros<mat>(n, n);
    mat K_dS = zeros<mat>(n, n);

    mat C = eye<mat>(n, n) - ones<mat>(n, n) / (double)n;

    // 1. Loop over ALL P variables
    for(int j = 0; j < P; j++) {
        // If j is in the drop list, skip it completely
        if(std::find(drop_vec.begin(), drop_vec.end(), j) != drop_vec.end()) {
            continue;
        }
        
        mat D = as<mat>(D_1D_list[j]);
        mat D2 = square(D);
        
        mat K_j_raw = exp(-D2 / rho);
        mat dK_j_raw = (D2 % K_j_raw) / (rho * rho);
        
        mat K_j = C * K_j_raw * C;
        mat dK_j = C * dK_j_raw * C;
        
        S += K_j;
        S2 += K_j % K_j;
        dS += dK_j;
        K_dS += K_j % dK_j;
    }

    // 2. Safely Construct the Null Space
    mat K_1D_null = zeros<mat>(n, n);
    mat dK_1D_null_drho = zeros<mat>(n, n);
    mat K_2D_null = zeros<mat>(n, n);
    mat dK_2D_null_drho = zeros<mat>(n, n);
    
    // Only build 1-D space if at least one background variable exists
    if (P_null > 0.5) {
        K_1D_null = S / P_null;
        dK_1D_null_drho = dS / P_null;
    }
    
    // Only build 2-D space if at least TWO background variables exist
    if (P_null > 1.5) {
        double M_2D_null = P_null * (P_null - 1.0) / 2.0;
        K_2D_null = 0.5 * (S % S - S2) / M_2D_null;
        dK_2D_null_drho = (S % dS - K_dS) / M_2D_null;
    }

    mat K_null = K_1D_null + K_2D_null;
    mat dK_null_drho = dK_1D_null_drho + dK_2D_null_drho;

    // 3. Construct Covariance and Gradients
    mat W0 = lam * K_null + eye<mat>(n, n);
    
    mat dW_dlogrho = lam * dK_null_drho * rho; 
    mat dW_dloglam = lam * K_null;

    return List::create(
        Named("W0") = W0,
        Named("dW_dlogrho") = dW_dlogrho,
        Named("dW_dloglam") = dW_dloglam
    );
}

#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eval_var_effects_H1_multi_split_fast(List D_1D_list, IntegerVector test_indices, 
                                          double log_rho, double log_lam_null, 
                                          double log_lam_test, int n) {
    double rho = exp(log_rho);
    double lam_null = exp(log_lam_null);
    double lam_test = exp(log_lam_test);
    
    int P = D_1D_list.size();
    int P_test_count = test_indices.size();
    
    double P_total = (double)P;
    double P_test = (double)P_test_count;
    double P_null = P_total - P_test;
    
    // Denominators for averaging
    double M_2D_total = (P_total * (P_total - 1.0)) / 2.0;
    
    double M_2D_null = 0.0;
    if (P_null > 1.5) {
        M_2D_null = (P_null * (P_null - 1.0)) / 2.0;
    }
    
    double M_2D_test = M_2D_total - M_2D_null;
    
    std::vector<int> test_vec = as<std::vector<int>>(test_indices);

    // Track Null and Test running sums separately
    mat S_null = zeros<mat>(n, n), S2_null = zeros<mat>(n, n);
    mat dS_null = zeros<mat>(n, n), K_dS_null = zeros<mat>(n, n);
    
    mat S_test = zeros<mat>(n, n), S2_test = zeros<mat>(n, n);
    mat dS_test = zeros<mat>(n, n), K_dS_test = zeros<mat>(n, n);

    mat C = eye<mat>(n, n) - ones<mat>(n, n) / (double)n;

    // 1. Loop over ALL P variables
    for(int j = 0; j < P; j++) {
        mat D = as<mat>(D_1D_list[j]);
        mat D2 = square(D);
        
        mat K_j_raw = exp(-D2 / rho);
        mat dK_j_raw = (D2 % K_j_raw) / (rho * rho);
        
        mat K_j = C * K_j_raw * C;
        mat dK_j = C * dK_j_raw * C;
        
        if(std::find(test_vec.begin(), test_vec.end(), j) != test_vec.end()) {
            S_test += K_j;
            S2_test += K_j % K_j;
            dS_test += dK_j;
            K_dS_test += K_j % dK_j;
        } else {
            S_null += K_j;
            S2_null += K_j % K_j;
            dS_null += dK_j;
            K_dS_null += K_j % dK_j;
        }
    }

    // 2. Construct Total Sums for the Fast Trick
    mat S_total = S_null + S_test;
    mat S2_total = S2_null + S2_test;
    mat dS_total = dS_null + dS_test;
    mat K_dS_total = K_dS_null + K_dS_test;

    // 3. Construct 1-D Spaces (Fixed Expression Templates)
    mat K_1D_null = zeros<mat>(n, n);
    mat dK_1D_null_drho = zeros<mat>(n, n);
    if (P_null > 0.5) {
        K_1D_null = S_null / P_null;
        dK_1D_null_drho = dS_null / P_null;
    }
    
    mat K_1D_test = zeros<mat>(n, n);
    mat dK_1D_test_drho = zeros<mat>(n, n);
    if (P_test > 0.5) {
        K_1D_test = S_test / P_test;
        dK_1D_test_drho = dS_test / P_test;
    }

    // 4. Construct 2-D Spaces via Geometric Partition (Fixed Expression Templates)
    mat K_2D_total_sum = 0.5 * (S_total % S_total - S2_total);
    mat dK_2D_total_sum_drho = S_total % dS_total - K_dS_total;
    
    mat K_2D_null_sum = zeros<mat>(n, n);
    mat dK_2D_null_sum_drho = zeros<mat>(n, n);
    if (P_null > 1.5) {
        K_2D_null_sum = 0.5 * (S_null % S_null - S2_null);
        dK_2D_null_sum_drho = S_null % dS_null - K_dS_null;
    }
    
    mat K_2D_test_sum = K_2D_total_sum - K_2D_null_sum;
    mat dK_2D_test_sum_drho = dK_2D_total_sum_drho - dK_2D_null_sum_drho;

    // Average the 2-D Spaces
    mat K_2D_null = zeros<mat>(n, n);
    mat dK_2D_null_drho = zeros<mat>(n, n);
    if (M_2D_null > 0.5) {
        K_2D_null = K_2D_null_sum / M_2D_null;
        dK_2D_null_drho = dK_2D_null_sum_drho / M_2D_null;
    }

    mat K_2D_test = zeros<mat>(n, n);
    mat dK_2D_test_drho = zeros<mat>(n, n);
    if (M_2D_test > 0.5) {
        K_2D_test = K_2D_test_sum / M_2D_test;
        dK_2D_test_drho = dK_2D_test_sum_drho / M_2D_test;
    }

    // 5. Aggregate Final Components
    mat K_null = K_1D_null + K_2D_null;
    mat dK_null_drho = dK_1D_null_drho + dK_2D_null_drho;

    mat K_test = K_1D_test + K_2D_test;
    mat dK_test_drho = dK_1D_test_drho + dK_2D_test_drho;

    // 6. Final Covariance Matrix and Gradients
    mat W = lam_null * K_null + lam_test * K_test + eye<mat>(n, n);
    
    mat dW_dlogrho = (lam_null * dK_null_drho + lam_test * dK_test_drho) * rho; 
    mat dW_dloglam_null = lam_null * K_null;
    mat dW_dloglam_test = lam_test * K_test;

    return List::create(
        Named("W") = W,
        Named("dW_dlogrho") = dW_dlogrho,
        Named("dW_dloglam_null") = dW_dloglam_null,
        Named("dW_dloglam_test") = dW_dloglam_test
    );
}

#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eval_ss_anova_score_components(List D_1D_list, IntegerVector test_indices, 
                                    double log_rho, double log_lam_null, int n) {
  double rho = exp(log_rho);
  double lam_null = exp(log_lam_null);
  
  int P = D_1D_list.size();
  int P_test_count = test_indices.size();
  
  double P_total = (double)P;
  double P_test = (double)P_test_count;
  double P_null = P_total - P_test;
  
  double M_2D_total = (P_total * (P_total - 1.0)) / 2.0;
  double M_2D_null = 0.0;
  if (P_null > 1.5) {
    M_2D_null = (P_null * (P_null - 1.0)) / 2.0;
  }
  double M_2D_test = M_2D_total - M_2D_null;
  
  std::vector<int> test_vec = as<std::vector<int>>(test_indices);
  
  mat S_null = zeros<mat>(n, n), S2_null = zeros<mat>(n, n);
  mat dS_null = zeros<mat>(n, n), K_dS_null = zeros<mat>(n, n);
  
  mat S_test = zeros<mat>(n, n), S2_test = zeros<mat>(n, n);
  mat dS_test = zeros<mat>(n, n), K_dS_test = zeros<mat>(n, n);
  
  mat C = eye<mat>(n, n) - ones<mat>(n, n) / (double)n;
  
  for(int j = 0; j < P; j++) {
    mat D = as<mat>(D_1D_list[j]);
    mat D2 = square(D);
    
    mat K_j_raw = exp(-D2 / rho);
    mat dK_j_raw = (D2 % K_j_raw) / (rho * rho);
    
    mat K_j = C * K_j_raw * C;
    mat dK_j = C * dK_j_raw * C;
    
    if(std::find(test_vec.begin(), test_vec.end(), j) != test_vec.end()) {
      S_test += K_j;
      S2_test += K_j % K_j;
    } else {
      S_null += K_j;
      S2_null += K_j % K_j;
      dS_null += dK_j;
      K_dS_null += K_j % dK_j;
    }
  }
  
  mat S_total = S_null + S_test;
  mat S2_total = S2_null + S2_test;
  
  mat K_1D_null = zeros<mat>(n, n);
  mat dK_1D_null_drho = zeros<mat>(n, n);
  if (P_null > 0.5) {
    K_1D_null = S_null / P_null;
    dK_1D_null_drho = dS_null / P_null;
  }
  
  mat K_1D_test = zeros<mat>(n, n);
  if (P_test > 0.5) {
    K_1D_test = S_test / P_test;
  }
  
  mat K_2D_total_sum = 0.5 * (S_total % S_total - S2_total);
  mat K_2D_null_sum = zeros<mat>(n, n);
  mat dK_2D_null_sum_drho = zeros<mat>(n, n);
  if (P_null > 1.5) {
    K_2D_null_sum = 0.5 * (S_null % S_null - S2_null);
    dK_2D_null_sum_drho = S_null % dS_null - K_dS_null;
  }
  
  mat K_2D_test_sum = K_2D_total_sum - K_2D_null_sum;
  
  mat K_2D_null = zeros<mat>(n, n);
  mat dK_2D_null_drho = zeros<mat>(n, n);
  if (M_2D_null > 0.5) {
    K_2D_null = K_2D_null_sum / M_2D_null;
    dK_2D_null_drho = dK_2D_null_sum_drho / M_2D_null;
  }
  
  mat K_2D_test = zeros<mat>(n, n);
  if (M_2D_test > 0.5) {
    K_2D_test = K_2D_test_sum / M_2D_test;
  }
  
  mat K_null = K_1D_null + K_2D_null;
  mat dK_null_drho = dK_1D_null_drho + dK_2D_null_drho;
  mat K_test = K_1D_test + K_2D_test;
  
  mat W0 = lam_null * K_null + eye<mat>(n, n);
  mat dW0_dlogrho = lam_null * dK_null_drho * rho; 
  mat dW0_dloglam = lam_null * K_null;
  
  return List::create(
    Named("W0") = W0,
    Named("dW0_dlogrho") = dW0_dlogrho,
    Named("dW0_dloglam") = dW0_dloglam,
    Named("K_test") = K_test
  );
}