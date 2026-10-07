# =========================================================
# FILE: competitor_analyses.R
# =========================================================
run_competitor_analyses <- function(Y, X, S, coords_sp,
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
  coords_sp <- as.matrix(coords_sp)
  p <- ncol(S)
  
  if (is.null(colnames(S))) colnames(S) <- paste0("S", 1:p)
  
  # Standardize S for BKMR (Standard practice for Gaussian Processes)
  S_scaled <- scale(S)
  
  # Initialize results list
  results <- list()
  
  # =========================================================
  # 1. Generalized Least Squares (GLS) via GpGp
  # =========================================================
  cat("[1/5] Running Spatial GLS (GpGp)...\n")
  results$GLS <- tryCatch({
    require(GpGp)
    
    # GpGp expects a design matrix. We combine X (fixed effects) and S (predictors of interest)
    X_gls <- cbind(X, S)
    
    # Fit the spatial model (Exponential Isotropic is standard for basic spatial autocorrelation)
    fit_gls <- fit_model(y = Y, locs = coords_sp, X = X_gls, covfun_name = "matern_isotropic", silent = TRUE)
    
    # Calculate Wald p-values for the fixed effects
    betas <- fit_gls$betahat
    se_betas <- fit_gls$sebeta
    p_vals <- 2 * pnorm(-abs(betas / se_betas))
    
    # Format and extract only the S variables
    res_gls <- data.frame(Estimate = betas, StdError = se_betas, P_value = p_vals)
    rownames(res_gls) <- colnames(X_gls)
    
    res_gls[colnames(S), ]
  }, error = function(e) {
    message("GLS Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 2. Generalized Additive Model (GAM) with Spatial Smooth
  # =========================================================
  cat("[2/5] Running GAM...\n")
  results$GAM <- tryCatch({
    require(mgcv)
    df_gam <- data.frame(Y = Y, S)
    
    # Safely extract latitude and longitude from the coords matrix
    df_gam$lat <- as.numeric(coords_sp[, 1])
    df_gam$lon <- as.numeric(coords_sp[, 2])
    
    if (!is.null(colnames(coords_sp))) {
      if ("lat" %in% tolower(colnames(coords_sp))) df_gam$lat <- as.numeric(coords_sp[, grep("lat", tolower(colnames(coords_sp)), ignore.case=TRUE)[1]])
      if ("lon" %in% tolower(colnames(coords_sp))) df_gam$lon <- as.numeric(coords_sp[, grep("lon", tolower(colnames(coords_sp)), ignore.case=TRUE)[1]])
    }
    
    smooth_terms <- paste(paste0("s(", colnames(S), ")"), collapse = " + ")
    spatial_terms <- "ti(lat, k=100) + ti(lon, k=100) + ti(lat, lon, k=20)"
    
    gam_form <- as.formula(paste("Y ~ X - 1 +", smooth_terms, "+", spatial_terms))
    fit_gam <- mgcv::gam(gam_form, data = df_gam, method = "REML")
    
    summary(fit_gam)$s.table
  }, error = function(e) {
    message("GAM Failed: ", e$message)
    return(NULL)
  })
  
  # =========================================================
  # 3. Spatial KMR (Marginal Variable Importance) - SEQUENTIAL
  # =========================================================
  cat("[3/5] Running Spatial KMR (Sequentially)...\n")
  results$Spatial_KMR <- tryCatch({
    kmr_pvals <- numeric(p)
    
    for (j in 1:p) {
      S_reord <- S[, c(j, setdiff(1:p, j)), drop = FALSE]
      test_kmr <- Kernel.test.Gauss.sp(
        yy = Y, 
        zz = t(S_reord), 
        xx = t(X),
        coords_sp = coords_sp,
        lam.min = .0001, rho.min = .0001, sig.min = .0001,
        maxiter = 500, rho.max = 500, lam.max = 500, 
        jitter = FALSE
      )
      kmr_pvals[j] <- test_kmr$pvalue[1]
    }
    
    names(kmr_pvals) <- colnames(S)
    kmr_pvals
  }, error = function(e) {
    message("Spatial KMR Failed: ", e$message)
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