#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mean_kernel_pairwise_rho(List mat_list, double rho) {
  int n_vars = mat_list.size();
  if (n_vars < 2) return arma::mat();

  // SAFETY 1: Check first matrix to get dimensions
  if (TYPEOF(mat_list[0]) != REALSXP) stop("Input matrices must be numeric (double).");
  NumericMatrix first_m = as<NumericMatrix>(mat_list[0]);
  const int n_rows = first_m.nrow();
  const int n_cols = first_m.ncol();
  const int total_elements = n_rows * n_cols;
  
  // OPTIMIZATION: Pre-calculate constants
  const double inv_rho = 1.0 / rho;
  
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
      double sum_K = 0.0;
      double sum_sq_K = 0.0;
      
      for (int k = 0; k < n_vars; k++) {
        double d = ptrs[k][idx];
        double val = std::exp(-(d * d * inv_rho)); // Kernel: exp(-d^2 / rho)
        
        sum_K += val;
        sum_sq_K += (val * val);
      }
      
      // Algebra Trick
      out_ptr[idx] = 0.5 * ((sum_K * sum_K) - sum_sq_K);
    }
  }

  double count = (double)n_vars * (n_vars - 1) / 2.0;
  return result / count;
}