#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec grad_loglik_cpp(arma::vec theta, 
                          arma::vec Y, 
                          arma::mat xx, 
                          arma::mat V, 
                          arma::mat K, 
                          arma::mat dK) {
  
  double lambda = theta[0];
  // double rho = theta[1]; // Implicitly used in K and dK
  // double sigma2 = theta[2]; // Implicitly used in V
  int n = Y.n_elem;
  
  // 1. Fast Symmetric Positive Definite Inverse
  // inv_sympd uses Cholesky decomposition, which is 2x faster than regular solve()
  arma::mat iV = arma::inv_sympd(V);
  
  // 2. Pre-compute common grouping terms
  // Let U = iV * xx^T
  arma::mat U = iV * xx.t(); 
  
  // Let W = xx * iV * xx^T
  arma::mat W = xx * U;      
  arma::mat invW = arma::inv_sympd(W); // tr1 in your R code
  
  // Let z = iV * Y
  arma::vec z = iV * Y; 
  
  // ==========================================
    // 1. Gradient w.r.t lambda (dV = K)
  // ==========================================
    // trace(iV * K) using the fast trace trick
  double term1_lam = arma::accu(iV % K); 
  
  // trace(invW * (U^T * K * U))
  arma::mat UKU = U.t() * K * U;
  double term2_lam = -arma::accu(invW % UKU);
  
  // Y^T * (-iV * K * iV) * Y  =>  -z^T * K * z
  double term3_lam = -arma::as_scalar(z.t() * K * z);
  
  double grad_lam = -0.5 * (term1_lam + term2_lam + term3_lam);
  
  // ==========================================
    // 2. Gradient w.r.t rho (dV = lambda * dK)
  // ==========================================
    arma::mat dV_rho = lambda * dK;
  
  double term1_rho = arma::accu(iV % dV_rho);
  
  arma::mat UdVU = U.t() * dV_rho * U;
  double term2_rho = -arma::accu(invW % UdVU);
  
  double term3_rho = -arma::as_scalar(z.t() * dV_rho * z);
  
  double grad_rho = -0.5 * (term1_rho + term2_rho + term3_rho);
  
  // ==========================================
    // 3. Gradient w.r.t sigma2 (dV = I)
  // ==========================================
    // trace(iV * I) = trace(iV)
  double term1_sig = arma::trace(iV);
  
  // trace(invW * (U^T * I * U)) = trace(invW * (U^T * U))
  arma::mat UtU = U.t() * U;
  double term2_sig = -arma::accu(invW % UtU);
  
  // -z^T * I * z = -z^T * z
  double term3_sig = -arma::dot(z, z);
  
  double grad_sigma2 = -0.5 * (term1_sig + term2_sig + term3_sig);
  
  // ==========================================
    // Return scaled gradient
  // ==========================================
    arma::vec gr = {grad_lam, grad_rho, grad_sigma2};
  
  return gr * (-1.0 / n);
}


#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double obj_loglik_cpp(arma::vec Y, 
                      arma::mat xx, 
                      arma::mat V) {
  
  int n = Y.n_elem;
  
  // 1. Fast Symmetric Positive Definite Inverse
  arma::mat iV = arma::inv_sympd(V);
  
  // 2. Safely compute log(det(V)) 
  // This prevents numerical overflow/Inf that occurs with standard det()
  double log_det_V;
  double sign_V;
  arma::log_det(log_det_V, sign_V, V);
  
  // 3. Compute W = xx * iV * xx^T  (Note: xx is P x N based on your R code)
  arma::mat W = xx * iV * xx.t();
  
  // 4. Safely compute log(det(W))
  double log_det_W;
  double sign_W;
  arma::log_det(log_det_W, sign_W, W);
  
  // 5. Compute the quadratic form: Y^T * iV * Y
  // as_scalar converts the 1x1 matrix result into a standard double
  double y_iV_y = arma::as_scalar(Y.t() * iV * Y);
  
  // 6. Combine into the final objective value
  double obj = log_det_V + log_det_W + y_iV_y;
  
  // Scale by 1 / (2n) to match the scaling in your gradient function
  return obj / (2.0 * n);
}

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List se_kernel_and_deriv_matrix_cpp(NumericMatrix D, double rho) {
  
  const int n_rows = D.nrow();
  const int n_cols = D.ncol();
  const int total_elements = n_rows * n_cols;
  
  // Get pointer to the input distance matrix
  const double* d_ptr = D.begin();
  
  // Initialize output matrices
  arma::mat K(n_rows, n_cols);
  arma::mat dK(n_rows, n_cols);
  double* K_ptr = K.memptr();
  double* dK_ptr = dK.memptr();
  
  // Pre-calculate constants for speed
  const double inv_rho = 1.0 / rho;
  const double scale_factor = 1.0 / (rho * rho);

  const int CHUNK_SIZE = 4096;

  #pragma omp parallel for schedule(static)
  for (int start_idx = 0; start_idx < total_elements; start_idx += CHUNK_SIZE) {
    int current_chunk_len = std::min(CHUNK_SIZE, total_elements - start_idx);
    
    for (int j = 0; j < current_chunk_len; j++) {
      int idx = start_idx + j;
      
      // 1. Get distance and square it
      double d = d_ptr[idx];
      double d2 = d * d; 
      
      // 2. Compute K = exp(-d^2 / rho)
      double k_val = std::exp(-d2 * inv_rho);
      K_ptr[idx] = k_val;
      
      // 3. Compute dK/drho = K * (d^2 / rho^2)
      dK_ptr[idx] = k_val * d2 * scale_factor;
    }
  }
  
  // Return both matrices in a named list
  return List::create(Named("K") = K, Named("dK") = dK);
}