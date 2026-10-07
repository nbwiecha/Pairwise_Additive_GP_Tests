#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List eval_ss_anova_interaction_fixed_rho(List D_1D_list, int idx1, int idx2, 
                                         double log_rho, double log_lam, int n) {
    double rho = exp(log_rho);
    double lam = exp(log_lam);
    
    int P = D_1D_list.size();
    double M_1D = (double)P;
    double M_2D_total = (P * (P - 1.0)) / 2.0;
    double M_2D_null = M_2D_total - 1.0;

    mat S = zeros<mat>(n, n);
    mat S2 = zeros<mat>(n, n);
    
    mat K_A = zeros<mat>(n, n);
    mat K_B = zeros<mat>(n, n);
    mat C = eye<mat>(n, n) - ones<mat>(n, n) / (double)n;

    // 1. Loop over variables to build centered 1-D kernels
    for(int j = 0; j < P; j++) {
        mat D = as<mat>(D_1D_list[j]);
        mat D2 = square(D);
        
        mat K_j_raw = exp(-D2 / rho);
        mat K_j = C * K_j_raw * C;
        
        if (j == idx1) K_A = K_j;
        if (j == idx2) K_B = K_j;
        
        S += K_j;
        S2 += K_j % K_j;
    }

    // 2. Construct Spaces
    mat K_1D = S / M_1D;
    mat K_2D_total_sum = 0.5 * (S % S - S2);
    
    // The exact interaction surface being tested
    mat K_test = K_A % K_B; 
    
    // The null background 2-D surface
    mat K_2D_null_sum = K_2D_total_sum - K_test;
    mat K_2D_null = zeros<mat>(n, n);
    if (M_2D_null > 0.5) {
        K_2D_null = K_2D_null_sum / M_2D_null;
    }

    // 3. Final Covariance and Derivatives
    mat K_null = K_1D + K_2D_null;
    mat W0 = lam * K_null + eye<mat>(n, n);
    mat dW0_dloglam = lam * K_null;

    return List::create(
        Named("W0") = W0,
        Named("dW0_dloglam") = dW0_dloglam,
        Named("K_test") = K_test
    );
}