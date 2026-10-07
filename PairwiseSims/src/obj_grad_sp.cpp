#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// ---------------------------------------------------------
// 1. Spatially Correlated Kernel Generation
// Computes both the covariate kernel (K_a) and spatial kernel (K_sp)
// ---------------------------------------------------------
// [[Rcpp::export]]
List se_kernel_and_deriv_matrix_sp_cpp(NumericMatrix D_cov, double rho_a, 
                                       NumericMatrix D_sp, double rho_sp) {
  
  const int n_rows = D_cov.nrow();
  const int n_cols = D_cov.ncol();
  const int total_elements = n_rows * n_cols;
  
  const double* d_cov_ptr = D_cov.begin();
  const double* d_sp_ptr = D_sp.begin();
  
  arma::mat K_a(n_rows, n_cols);
  arma::mat dK_a(n_rows, n_cols);
  arma::mat K_sp(n_rows, n_cols);
  arma::mat dK_sp(n_rows, n_cols);
  
  double* K_a_ptr = K_a.memptr();
  double* dK_a_ptr = dK_a.memptr();
  double* K_sp_ptr = K_sp.memptr();
  double* dK_sp_ptr = dK_sp.memptr();
  
  const double inv_rho_a = 1.0 / rho_a;
  const double scale_factor_a = 1.0 / (rho_a * rho_a);
  const double inv_rho_sp = 1.0 / rho_sp;
  const double scale_factor_sp = 1.0 / (rho_sp * rho_sp);
  
  const int CHUNK_SIZE = 4096;
  
#pragma omp parallel for schedule(static)
  for (int start_idx = 0; start_idx < total_elements; start_idx += CHUNK_SIZE) {
    int current_chunk_len = std::min(CHUNK_SIZE, total_elements - start_idx);
    
    for (int j = 0; j < current_chunk_len; j++) {
      int idx = start_idx + j;
      
      // Covariate Kernel
      double d_cov2 = d_cov_ptr[idx] * d_cov_ptr[idx];
      double k_a_val = std::exp(-d_cov2 * inv_rho_a);
      K_a_ptr[idx] = k_a_val;
      dK_a_ptr[idx] = k_a_val * d_cov2 * scale_factor_a;
      
      // Spatial Kernel
      double d_sp2 = d_sp_ptr[idx] * d_sp_ptr[idx];
      double k_sp_val = std::exp(-d_sp2 * inv_rho_sp);
      K_sp_ptr[idx] = k_sp_val;
      dK_sp_ptr[idx] = k_sp_val * d_sp2 * scale_factor_sp;
    }
  }
  
  return List::create(
    Named("K_a") = K_a, 
    Named("dK_a") = dK_a,
    Named("K_sp") = K_sp, 
    Named("dK_sp") = dK_sp
  );
}

// ---------------------------------------------------------
// 2. Spatial Objective Likelihood
// ---------------------------------------------------------
// [[Rcpp::export]]
double obj_loglik_sp_cpp(arma::vec Y, 
                         arma::mat xx, 
                         arma::mat V) { 
  
  int n = Y.n_elem;
  arma::mat iV = arma::inv_sympd(V);
  
  double log_det_V;
  double sign_V;
  arma::log_det(log_det_V, sign_V, V);
  
  arma::mat W = xx * iV * xx.t();
  
  double log_det_W;
  double sign_W;
  arma::log_det(log_det_W, sign_W, W);
  
  double y_iV_y = arma::as_scalar(Y.t() * iV * Y);
  double obj = log_det_V + log_det_W + y_iV_y;
  
  return obj / (2.0 * n);
}

// ---------------------------------------------------------
// 3. Spatial Gradient Function (5 Parameters)
// theta = [lam_a, rho_a, lam_sp, rho_sp, sigma2]
// ---------------------------------------------------------
// [[Rcpp::export]]
arma::vec grad_loglik_sp_cpp(arma::vec theta, 
                             arma::vec Y, 
                             arma::mat xx, 
                             arma::mat V, 
                             arma::mat K_a, 
                             arma::mat dK_a,
                             arma::mat K_sp, 
                             arma::mat dK_sp) {
  
  double lam_a = theta[0];
  double lam_sp = theta[2];
  int n = Y.n_elem;
  
  arma::mat iV = arma::inv_sympd(V);
  arma::mat U = iV * xx.t(); 
  arma::mat W = xx * U;      
  arma::mat invW = arma::inv_sympd(W); 
  arma::vec z = iV * Y; 
  
  auto compute_grad_component = [&](const arma::mat& dV_dParam) {
    double term1 = arma::accu(iV % dV_dParam); 
    arma::mat UdVU = U.t() * dV_dParam * U;
    double term2 = -arma::accu(invW % UdVU);
    double term3 = -arma::as_scalar(z.t() * dV_dParam * z);
    return -0.5 * (term1 + term2 + term3);
  };
  
  double grad_lam_a = compute_grad_component(K_a);
  
  arma::mat dV_rho_a = lam_a * dK_a;
  double grad_rho_a = compute_grad_component(dV_rho_a);
  
  double grad_lam_sp = compute_grad_component(K_sp);
  
  arma::mat dV_rho_sp = lam_sp * dK_sp;
  double grad_rho_sp = compute_grad_component(dV_rho_sp);
  
  double term1_sig = arma::trace(iV);
  arma::mat UtU = U.t() * U;
  double term2_sig = -arma::accu(invW % UtU);
  double term3_sig = -arma::dot(z, z);
  double grad_sigma2 = -0.5 * (term1_sig + term2_sig + term3_sig);
  
  arma::vec gr = {grad_lam_a, grad_rho_a, grad_lam_sp, grad_rho_sp, grad_sigma2};
  return gr * (-1.0 / n);
}