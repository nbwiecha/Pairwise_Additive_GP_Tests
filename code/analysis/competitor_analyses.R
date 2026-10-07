run_competitor_analyses <- function(Y, X, S, 
                                    bkmr_iter = 10000, 
                                    bkmr_knots = NULL, 
                                    bkmr_groups = NULL,
                                    cl = NULL) {
  
  # Auto-detect cluster object in the global environment if not explicitly provided
  if (is.null(cl) && exists("cl", envir = .GlobalEnv)) {
    cl <- get("cl", envir = .GlobalEnv)
  }
  
  # Ensure standard formats
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  S <- as.matrix(S)
  p <- ncol(S)
  
  if (is.null(colnames(S))) colnames(S) <- paste0("S", 1:p)
  
  # Standardize S for BKMR (Standard practice for Gaussian Processes)
  S_scaled <- scale(S)
  
  # Initialize results list
  results <- list()
  
  # =========================================================
  # 1. Ordinary Least Squares (OLS)
  # =========================================================
  cat("[1/5] Running OLS...\n")
  results$OLS <- tryCatch({
    df_ols <- data.frame(Y = Y, X, S)
    fit_ols <- lm(Y ~ . - 1, data = df_ols)
    summary(fit_ols)$coefficients
  }, error = function(e) {
    message("OLS Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 2. Generalized Additive Model (GAM)
  # =========================================================
  cat("[2/5] Running GAM...\n")
  results$GAM <- tryCatch({
    df_gam <- data.frame(Y = Y, S)
    smooth_terms <- paste(paste0("s(", colnames(S), ")"), collapse = " + ")
    gam_form <- as.formula(paste("Y ~ X - 1 +", smooth_terms))
    fit_gam <- mgcv::gam(gam_form, data = df_gam, method = "REML")
    summary(fit_gam)$s.table
  }, error = function(e) {
    message("GAM Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 3. Previous KMR (Marginal Variable Importance) - PARALLELIZED
  # =========================================================
  cat("[3/5] Running Previous KMR...\n")
  results$Previous_KMR <- tryCatch({
    
    if (!is.null(cl)) {
      cat("      -> Utilizing parallel cluster for KMR...\n")
      # Ensure worker nodes have the required KMR functions
      # 1. Export the function and environment
      parallel::clusterExport(cl, varlist = c("Kernel.test.Gauss", "Kernel.reg.Gauss", 
                                              "kernel.Gauss", "obj.loglik.se", 
                                              "grad.loglik.se", "Trace", "se_kernel_and_deriv_matrix_cpp",
                                              "obj_loglik_cpp", "grad_loglik_cpp"), 
                              envir = .GlobalEnv)
      
      # 2. Use parSapply with a cleaner index
      # We define the indices as an explicit vector first
      idx_vec <- 1:p
      
      kmr_pvals <- parallel::parSapply(cl, X=idx_vec, FUN=function(j, Y, X1, S, p) {
        S_reord <- S[, c(j, setdiff(1:p, j))]
        test_kmr <- Kernel.test.Gauss(
          yy = Y, 
          zz = t(S_reord), 
          xx = t(X), 
          lam.min = .0001, rho.min = .0001, sig.min = .0001,
          maxiter = 500, rho.max = 500, lam.max = 500, 
          jitter = FALSE
        )
        return(test_kmr$pvalue[1])
      }, Y = Y, X1 = X, S = S, p = p)
      
    } else {
      cat("      -> No cluster detected. Running sequentially...\n")
      kmr_pvals <- numeric(p)
      for (j in 1:p) {
        S_reord <- S[, c(j, setdiff(1:p, j))]
        test_kmr <- Kernel.test.Gauss(
          yy = Y, 
          zz = t(S_reord), 
          xx = t(X), 
          lam.min = .0001, rho.min = .0001, sig.min = .0001,
          maxiter = 500, rho.max = 500, lam.max = 500, 
          jitter = FALSE
        )
        kmr_pvals[j] <- test_kmr$pvalue[1]
      }
    }
    
    names(kmr_pvals) <- colnames(S)
    kmr_pvals
  }, error = function(e) {
    message("Previous KMR Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 4. BKMR (Individual Selection)
  # =========================================================
  cat("[4/5] Running BKMR (Individual)...\n")
  results$BKMR_Individual <- tryCatch({
    
    knots_mat <- NULL
    if (!is.null(bkmr_knots)) {
      knots_mat <- fields::cover.design(S_scaled, nd = bkmr_knots)$design
    }
    
    fit_bkmr_indiv <- bkmr::kmbayes(
      y = Y, 
      Z = S_scaled, 
      X = X, 
      iter = bkmr_iter, 
      varsel = TRUE, 
      knots = knots_mat,
      verbose = FALSE
    )
    
    bkmr::ExtractPIPs(fit_bkmr_indiv)
  }, error = function(e) {
    message("BKMR (Individual) Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 5. BKMR (Grouped Selection)
  # =========================================================
  if (!is.null(bkmr_groups)) {
    cat("[5/5] Running BKMR (Grouped)...\n")
    results$BKMR_Grouped <- tryCatch({
      
      knots_mat <- NULL
      if (!is.null(bkmr_knots)) {
        knots_mat <- fields::cover.design(S_scaled, nd = bkmr_knots)$design
      }
      
      fit_bkmr_grouped <- bkmr::kmbayes(
        y = Y, 
        Z = S_scaled, 
        X = X, 
        iter = bkmr_iter, 
        groups = bkmr_groups,
        varsel = TRUE, 
        knots = knots_mat,
        verbose = FALSE
      )
      
      pips <- bkmr::ExtractPIPs(fit_bkmr_grouped)
      unique(pips[, c("group", "groupPIP")])
      
    }, error = function(e) {
      message("BKMR (Grouped) Failed: ", e$message)
      return(NULL)
    })
  } else {
    cat("[5/5] Skipping BKMR (Grouped) - No groups provided.\n")
  }
  
  cat("Done.\n")
  return(results)
}