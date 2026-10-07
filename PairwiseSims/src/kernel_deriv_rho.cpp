#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat deriv_kernel_full_rho(List mat_list, double rho) {
  int n_vars = mat_list.size();
  if (n_vars < 1) return arma::mat();

  // SAFETY 1: Check the first matrix to get dimensions
  if (TYPEOF(mat_list[0]) != REALSXP) stop("Input matrices must be numeric (double).");
  NumericMatrix first_m = as<NumericMatrix>(mat_list[0]);
  const int n_rows = first_m.nrow();
  const int n_cols = first_m.ncol();
  const int total_elements = n_rows * n_cols;
  
  // OPTIMIZATION: Pre-calculate constants
  // Kernel is exp(-d^2 / rho)
  // Derivative factor w.r.t rho is (d^2 / rho^2) * Kernel
  const double inv_rho = 1.0 / rho;
  const double scale_factor = 1.0 / (rho * rho); // 1 / rho^2
  
  // SAFETY 2: Validate ALL matrices and get pointers
  std::vector<const double*> ptrs(n_vars);
  for(int k = 0; k < n_vars; k++) {
    SEXP s = mat_list[k];
    if (TYPEOF(s) != REALSXP) stop("Matrix %d is not numeric.", k + 1);
    ptrs[k] = as<NumericMatrix>(s).begin();
  }

  arma::mat result(n_rows, n_cols);
  double* out_ptr = result.memptr();
  
  // OPTIMIZATION: Cache Strip-Mining
  const int CHUNK_SIZE = 4096;

  #pragma omp parallel for schedule(static)
  for (int start_idx = 0; start_idx < total_elements; start_idx += CHUNK_SIZE) {
    int current_chunk_len = std::min(CHUNK_SIZE, total_elements - start_idx);
    
    // OPTIMIZATION: Vectorizable inner loop
    for (int j = 0; j < current_chunk_len; j++) {
      int idx = start_idx + j;
      
      double sum_E = 0.0;     // Sum of K
      double sum_ED = 0.0;    // Sum of K * d^2
      double sum_E2D = 0.0;   // Sum of K^2 * d^2
      
      for (int k = 0; k < n_vars; k++) {
        double d = ptrs[k][idx];
        double d2 = d * d;
        
        double val = std::exp(-(d2 * inv_rho));
        double val_d2 = val * d2;
        
        sum_E   += val;
        sum_ED  += val_d2;
        sum_E2D += val * val_d2; 
      }
      
      // Derivative Logic (Algebra Trick):
      // The derivative of the pairwise sum is: (Sum_E * Sum_ED) - Sum_E2D
      // The derivative of the singles sum is: Sum_ED
      // Total Raw Derivative = Pairwise + Singles
      
      double pairwise_part = (sum_E * sum_ED) - sum_E2D;
      double total_raw = pairwise_part + sum_ED;
      
      out_ptr[idx] = total_raw * scale_factor;
    }
  }

  // Normalize by the total count of kernels (Pairs + Singles)
  double total_count = (double)n_vars * (n_vars + 1) / 2.0;
  
  return result / total_count;
}

// [[Rcpp::export]]
arma::mat deriv_kernel_pairwise_rho(List mat_list, double rho) {
  int n_vars = mat_list.size();
  
  // SAFETY 1: Need at least 2 variables for pairwise interactions
  if (n_vars < 2) stop("At least 2 matrices are required for pairwise interactions.");
  
  // SAFETY 2: Check the first matrix to get dimensions
  if (TYPEOF(mat_list[0]) != REALSXP) stop("Input matrices must be numeric (double).");
  NumericMatrix first_m = as<NumericMatrix>(mat_list[0]);
  const int n_rows = first_m.nrow();
  const int n_cols = first_m.ncol();
  const int total_elements = n_rows * n_cols;
  
  // OPTIMIZATION: Pre-calculate constants
  const double inv_rho = 1.0 / rho;
  const double scale_factor = 1.0 / (rho * rho); // 1 / rho^2
  
  // SAFETY 3: Validate ALL matrices and get pointers
  std::vector<const double*> ptrs(n_vars);
  for(int k = 0; k < n_vars; k++) {
    SEXP s = mat_list[k];
    if (TYPEOF(s) != REALSXP) stop("Matrix %d is not numeric.", k + 1);
    ptrs[k] = as<NumericMatrix>(s).begin();
  }
  
  arma::mat result(n_rows, n_cols);
  double* out_ptr = result.memptr();
  
  // OPTIMIZATION: Cache Strip-Mining
  const int CHUNK_SIZE = 4096;
  
#pragma omp parallel for schedule(static)
  for (int start_idx = 0; start_idx < total_elements; start_idx += CHUNK_SIZE) {
    int current_chunk_len = std::min(CHUNK_SIZE, total_elements - start_idx);
    
    // OPTIMIZATION: Vectorizable inner loop
    for (int j = 0; j < current_chunk_len; j++) {
      int idx = start_idx + j;
      
      double sum_E = 0.0;     // Sum of K
      double sum_ED = 0.0;    // Sum of K * d^2
      double sum_E2D = 0.0;   // Sum of K^2 * d^2
      
      for (int k = 0; k < n_vars; k++) {
        double d = ptrs[k][idx];
        double d2 = d * d;
        
        double val = std::exp(-(d2 * inv_rho));
        double val_d2 = val * d2;
        
        sum_E   += val;
        sum_ED  += val_d2;
        sum_E2D += val * val_d2; 
      }
      
      // Derivative Logic (Algebra Trick for PAIRS ONLY):
      // The derivative of the pairwise sum is exactly: (Sum_E * Sum_ED) - Sum_E2D
      double pairwise_part = (sum_E * sum_ED) - sum_E2D;
      
      out_ptr[idx] = pairwise_part * scale_factor;
    }
  }
  
  // Normalize by the total count of PAIRWISE kernels
  // Count = P(P-1)/2
  double total_count = (double)n_vars * (n_vars - 1) / 2.0;
  
  return result / total_count;
}