#################
##### Kernel machine testing for Gaussian kernel
##### Using pairwise interactive gaussian kernel
###################

# Original source: https://www4.stat.ncsu.edu/~amaity//software/KMTest_garrote_header.r
# Archived on Wayback Machine
# Nate Wiecha changes: use pairwise Gaussian kernel
# For legacy compatibility using same input dimensions as original


create_d_mat_list_1d <- function(X1, X2){
  require(fields)
  p <- ncol(X1)
  
  d.list <- vector(mode="list", length=p)
  for(i in 1:length(d.list)){
    d.list[[i]] = rdist(X1[,i], X2[,i])
  }
  
  return(d.list)
}

kernel.Gauss.pairwise = mean_kernel_pairwise_rho


Kernel.reg.Gauss.pairwise = function(yy, zz, xx, zz.p=NULL, null=FALSE,
                                     starting_val=NULL, maxiter=200,
                                     lam.min=0.0001, lam.max=60,
                                     rho.min=0.0001, rho.max=40,
                                     sig.min=.0001, sig.max=10,
                                     criter.use="loglik", plot_conv=TRUE){
  
  
  
  kernel.fn <- kernel.Gauss.pairwise
  
  # Modifications by Nate Wiecha for simulation:
  # allow prediction at points zz.p
  # assuming only covariate is an intercept, ie. xx=(1, ..., 1)
  # allow 100 iterations
  
  ### This code does kernel machine regression for pairwise gaussian kernel
  ## yy: response nx1
  ## zz: matrix modeled nonparametrically (each col is a sample)
  ## xx: other variables (each col is a sample)
  ##
  ##Estimates: h(z), beta and rho (the kernel parameter)
  require(fields)
  if(is.null(zz.p)){
    zz.p<-zz
  }
  n = length(yy)
  
  ## iterative algo
  # dMat = sqare.diff.mat(zz)
  dMat.list = create_d_mat_list_1d(t(zz), t(zz))
  dMat.list.p = create_d_mat_list_1d(t(zz.p), t(zz))
  dMat.list.pp = create_d_mat_list_1d(t(zz.p), t(zz.p))
  
  if(null){
    
    dMat.list[[1]] <- matrix(0, 
                             ncol=ncol(dMat.list[[1]]), nrow=nrow(dMat.list[[1]]))
    
    dMat.list.p[[1]] <- matrix(0, 
                             ncol=ncol(dMat.list.p[[1]]), nrow=nrow(dMat.list.p[[1]]))
    
    dMat.list.pp[[1]] <- matrix(0, 
                             ncol=ncol(dMat.list.pp[[1]]), nrow=nrow(dMat.list.pp[[1]]))
  }
  # dMat.p = rdist(t(zz.p), t(zz))^2
  # dMat.pp = sqare.diff.mat(zz.p)
  
  bet.hat = matrix(.2, ncol=1, nrow=nrow(xx))
  ########### Step 1(a): fix sigma and rho and estimate lambda
  sig.hat = 1
  rho.hat = 1
  lam.hat = 0.01
  
  if(!is.null(starting_val)){
    bet.hat = starting_val$bet
    sig.hat = starting_val$sig
    rho.hat = starting_val$rho
    lam.hat = starting_val$lam
  }
  
  Kmat = kernel.fn(dMat.list, rho.hat)
  
  conv.crit = 0.00001
  criter = 999
  bet.old = bet.hat
  lam.old = lam.hat
  sig.old = sig.hat
  rho.old = rho.hat
  niter = 0
  
  criter.track <- rep(NA, maxiter)
  
  loglik.prev <- obj.loglik.fast(c(lam.hat, rho.hat, sig.hat), Y=(yy - t(xx)%*%bet.hat),
                                 xx=xx, dMat.list=dMat.list)
  
  pb <- txtProgressBar(max=maxiter, style=3)
  while( (criter > conv.crit) & (niter < maxiter) ){
    setTxtProgressBar(pb, value=niter)
    #print(c(niter, c(bet.old), lam.old, sig.old, rho.old))
    ###########################################
    ###### STEP 1 :fix  sigma, lambda and rho and estimate beta and h().
    ###########################################
    Kmat = kernel.fn(dMat.list, rho.hat)
    V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
    iV = solve(V)
    bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%yy
    h.hat = lam.hat*Kmat%*%iV%*% (yy - t(xx)%*%bet.hat)	
    
    ###########################################
    ###### STEP 2: fix beta and estimate sigma, lambda and rho
    ###########################################
    # Px = solve(xx%*%iV%*%t(xx))
    # H.mat = (t(xx)%*%Px%*%xx)%*%iV + 
    #   lam.hat*Kmat%*%iV%*%( diag(1, nrow=n, ncol=n) - t(xx)%*%Px%*%xx%*%iV )
    # 
    # sig.hat = sum( (yy - t(xx)%*%bet.hat - h.hat)^2 )/ (n - Trace(H.mat))
    # if((n - Trace(H.mat)) < 0)	{	
    #   sig.hat = (sum((yy - t(xx)%*%bet.hat - h.hat)^2)/ (n))	
    # }
    
    theta.hat <- optim(c(lam.hat, rho.hat, sig.hat), fn=obj.loglik.fast, gr=grad.loglik.fast,
                       method="L-BFGS-B",
                       xx=xx,  Y=(yy - t(xx)%*%bet.hat),
                       dMat.list=dMat.list,
                       lower=c(lam.min, rho.min, sig.min), upper=c(lam.max, rho.max, sig.max))
    lam.hat <- theta.hat$par[1]
    rho.hat <- theta.hat$par[2]
    sig.hat <- theta.hat$par[3]
    
    loglik.curr <- obj.loglik.fast(c(lam.hat, rho.hat, sig.hat), Y=(yy - t(xx)%*%bet.hat),
                                   xx=xx, dMat.list=dMat.list)
    # 
    # lam.hat = optimize(loglik.lam, interval=c(lam.min,lam.max), 
    #                    Y=(yy - t(xx)%*%bet.hat), 
    #                    Kmat=Kmat,
    #                    xx = xx,
    #                    sigmasq = sig.hat)$minimum 
    # 
    # 
    # rho.hat = optimize(loglik.rho.pairwise, interval=c(rho.min,rho.max), 
    #                    Y=(yy - t(xx)%*%bet.hat), 
    #                    dMat.list=dMat.list,
    #                    xx = xx,
    #                    lambda = lam.hat,
    #                    sigmasq = sig.hat,
    #                    null=null)$minimum
    ##
    
    
    criter1 = (sum(abs(bet.hat - bet.old)) + 
      abs(lam.hat - lam.old) +
      abs(sig.hat - sig.old) +
      abs(rho.hat - rho.old))
    
    criter2 <- abs(loglik.curr - loglik.prev)
    
    criter3 = sum(abs(lam.hat - lam.old) +
                 abs(sig.hat - sig.old) +
                 abs(rho.hat - rho.old))
    

    
    loglik.prev <- loglik.curr
    
    
    criter <- criter1
    
    if(criter.use=="loglik"){
      criter <- criter2
    }
    
    if(criter.use=="theta"){
      criter <- criter3
    }
    
    
    criter.track[niter] <- criter
    
    
    
    # lik.val = loglik.lam(lambda = lam.hat, Y=(yy - t(xx)%*%bet.hat), Kmat=Kmat, xx=xx, sigmasq=sig.hat)			
    
    bet.old = bet.hat
    lam.old = lam.hat
    sig.old = sig.hat
    rho.old = rho.hat
    niter = niter + 1
    if(niter %% 10 == 0 & plot_conv){
      plot(criter.track[max(niter - 10, 0):niter], type="l")
    }
    
  }
  # print(niter)
  close(pb)
  ## sanity check
  if(niter >= maxiter){
    print("Did not converge")
    # stop("Did not converge")
  }
  
  Kmat = kernel.fn(dMat.list, rho.hat)
  V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
  iV = solve(V)
  
  np <- ncol(zz.p)
  Kmat.p = kernel.fn(dMat.list.p, rho.hat)
  Kmat.pp = kernel.fn(dMat.list.pp, rho.hat)
  
  bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%yy
  h.hat = lam.hat * Kmat.p %*% iV %*% (yy - t(xx) %*% bet.hat)
  h.sd  = sqrt(diag(
    
    lam.hat*Kmat.pp - lam.hat^2*Kmat.p %*% iV %*% t(Kmat.p)
    
  ))
  
  return(list(h.hat = h.hat,sig=sig.hat, lam=lam.hat, rho=rho.hat, 
              bet=bet.hat, h.sd=h.sd,  criter=criter.track, null=null))				
}


V_fn <- function(d.list.1d, sigma2, lambda, rho){
  
  n <- nrow(d.list.1d[[1]])
  K <- kernel.Gauss.pairwise(d.list.1d, rho)
  V <- sigma2*diag(n) + lambda*K
  return(V)
  
}
obj.loglik.fast <- function(theta, Y, dMat.list, xx) {
  
  # Extract parameters
  lambda <- theta[1]
  rho    <- theta[2]
  sigma2 <- theta[3]
  n <- length(Y)
  
  # 1. Build the V matrix using your existing function
  # V <- V_fn(dMat.list, sigma2, lambda, rho)
  K  <- mean_kernel_pairwise_rho(dMat.list, rho) 
  V <- sigma2*diag(n) + lambda*K
  V <- (V + t(V))/2
  # 2. Call the fast C++ objective function
  # (Make sure Y is a standard vector, not a 1-column matrix)
  val <- obj_loglik_cpp(as.numeric(Y), xx, V)
  
  return(val)
}

Trace = function(M){sum(diag(M))}

d.rho.kernel.Gauss.pairwise = deriv_kernel_pairwise_rho 

grad.loglik.fast <- function(theta, Y, dMat.list, xx) {
  lambda <- theta[1]
  rho    <- theta[2]
  sigma2 <- theta[3]
  
  n <- length(Y)
  # 1. Generate matrices using your existing functions
  # V  <- V_fn(dMat.list, sigma2, lambda, rho)
  
  # Call the fast C++ functions we wrote earlier!
  K  <- mean_kernel_pairwise_rho(dMat.list, rho) 
  V <- sigma2*diag(n) + lambda*K
  V <- (V + t(V))/2
  dK <- deriv_kernel_pairwise_rho(dMat.list, rho)
  
  # 2. Pass them to the newly compiled C++ gradient function
  # Ensure Y is a numeric vector, not a matrix
  gr <- grad_loglik_cpp(theta, as.numeric(Y), xx, V, K, dK)
  
  # Rcpp returns a matrix column by default, cast it to a standard vector
  return(as.vector(gr)) 
}



###################### 
##### Testing using garrote pairwise kernel
#############################
Kernel.test.Gauss.pairwise = function(yy, zz, xx, zz.p=NULL, 
                                      maxiter=200, starting_val=NULL,
                                      lam.min=0.0001, lam.max=60,
                                      rho.min=0.0001, rho.max=40,
                                      sig.min=0.0001, sig.max=10,
                                      jitter=TRUE,
                                      pca=FALSE){
  # test using kernel machine for the effect of the first row of zz
  
  n = length(yy)
  # if(is.null(zz.p)){
  #   zz.p <- zz
  # }
  ### estimation under H0
  est.h0 = Kernel.reg.Gauss.pairwise(yy=yy, zz=zz, #zz=zz[-1,,drop=FALSE], 
                                     xx=xx, zz.p=zz.p,
                                     maxiter=maxiter, starting_val=starting_val,
                                     lam.min=lam.min, lam.max=lam.max,
                                     rho.min=rho.min, rho.max=rho.max,
                                     sig.min=sig.min, sig.max=sig.max,
                                     null=TRUE)
  h.hat = est.h0[[1]]
  sig.hat = est.h0[[2]]
  lam.hat = est.h0[[3]]
  rho.hat = est.h0[[4]]
  bet.hat = est.h0[[5]]
  #print(bet.hat)
  
  # dMat = sqare.diff.mat(zz[-1,,drop=FALSE])
  dMat.list = create_d_mat_list_1d(t(zz), t(zz))
  dMat.list[[1]] <- matrix(0, 
                           ncol=ncol(dMat.list[[1]]), nrow=nrow(dMat.list[[1]]))
  
  Kmat = kernel.Gauss.pairwise(dMat.list, rho.hat)
  
  mat = sig.hat*diag(1, nrow=n, ncol=n) + lam.hat*Kmat
  solve.mat = solve(mat)
  
  eps.hat = yy - t(xx)%*%bet.hat - h.hat
  
  ##derivative of the garrote pairwise kernel wrt garrote
  d2.z1 <- rdist(t(zz[1,,drop=FALSE]))^2
  K.uni <- exp(-rdist(t(zz[2,,drop=FALSE]))^2/rho.hat)
  if(nrow(zz) >= 3){
    for(j in 3:nrow(zz)){
      K.uni <- K.uni + exp(-rdist(t(zz[j,,drop=FALSE]))^2/rho.hat)
    }
  }
  p <- nrow(zz)
  Kmat.del <- - d2.z1 * K.uni / rho.hat / (p*(p-1)/2)
  
  # mul.mat = matrix(1, ncol=1, nrow=n)%*%zz[1,,drop=FALSE]
  # mul.mat = -(mul.mat - t(mul.mat))^2
  # Kmat.del = mul.mat * Kmat / rho.hat
  
  ## score stat
  score.chi = lam.hat * t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% Kmat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2
  
  ## derivatives of V ('mat' in this code)
  mat.lam = Kmat
  mat.sig = diag(1, nrow=n, ncol=n)
  mat.del = lam.hat*Kmat.del
  mat.rho = lam.hat * d.rho.kernel.Gauss.pairwise(dMat.list, rho.hat)
  
  #Info matrix
  P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
  
  #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
  I0 = matrix(NA, 4,4)
  
  I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del
  
  I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
  I0[2,1] = I0[1,2]
  
  I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
  I0[3,1] = I0[1,3]
  
  I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
  I0[4,1] = I0[1,4]
  
  #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
  #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]
  
  I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig
  
  I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
  I0[3,2] = I0[2,3]
  
  I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
  I0[4,2] = I0[2,4]
  
  #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
  #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]
  
  I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam
  
  I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
  I0[4,3] = I0[3,4]
  
  #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
  #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]
  
  I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2
  
  # Add jitter to ensure invertible
  if(jitter){
    I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
  }
  
  
  #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
  #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
  
  #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
  
  
  #Info for delta
  tot.dim = ncol(I0)
  I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
  #print(I.deldel)
  
  ####
  ## Chi-square approximation
  md = lam.hat*Trace(Kmat.del%*%P0.mat)/2
  #print(md)
  
  m.chi = I.deldel / (2*md)
  d.chi = md / m.chi
  
  pval.chi = 1- pchisq(score.chi/m.chi,d.chi)
  
  
  score.pca <- score.pca2 <- pval.pca <- 
    pval.pca2 <- m.pca <- m.pca2 <- d.pca <- d.pca2 <- NA
  
  if(pca){
    ###################################################
    ############################################
    ### PCA based reconstruction
    ############################################
    ###################################################
    svd.out = eigen(Kmat.del) 
    Kmat.del.plus = svd.out$vectors %*% diag(abs(svd.out$values), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    #print(svd.out$d)
    #print(eigen(Kmat.del)$values)
    
    ## derivatives of V ('mat' in this code)
    mat.lam = Kmat
    mat.sig = diag(1, nrow=n, ncol=n)
    mat.del = lam.hat*Kmat.del.plus
    mat.rho = d.rho.kernel.Gauss.pairwise(dMat.list, rho.hat)
    
    ## score stat
    score.pca = t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2
    
    #print(c(score.chi, score.pca))
    
    
    
    #Info matrix
    P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
    
    
    
    #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
    I0 = matrix(NA, 4,4)
    
    I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del
    
    I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
    I0[2,1] = I0[1,2]
    
    I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
    I0[3,1] = I0[1,3]
    
    I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
    I0[4,1] = I0[1,4]
    
    #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]
    
    I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig
    
    I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
    I0[3,2] = I0[2,3]
    
    I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
    I0[4,2] = I0[2,4]
    
    #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]
    
    I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam
    
    I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
    I0[4,3] = I0[3,4]
    
    #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]
    
    I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2
    
    # Add jitter to ensure invertible
    if(jitter){
      I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
    }
    
    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
    
    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
    
    
    #Info for delta
    tot.dim = ncol(I0)
    I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
    #print(I.deldel)
    
    ####
    ## Chi-square approximation
    md = Trace(mat.del%*%P0.mat)/2
    #print(md)
    
    m.pca = I.deldel / (2*md)
    d.pca = md / m.pca
    
    pval.pca = 1- pchisq(score.pca/m.pca,d.pca)
    
    
    
    ###################################################
    ############################################
    ### PCA based reconstruction (METHOD 2)
    ############################################
    ###################################################
    #svd.out = eigen(Kmat.del) 
    Kmat.del.zero = svd.out$vectors %*% diag(svd.out$values*(svd.out$values>0), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    #print(svd.out$d)
    #print(eigen(Kmat.del)$values)
    
    ## derivatives of V ('mat' in this code)
    mat.lam = Kmat
    mat.sig = diag(1, nrow=n, ncol=n)
    mat.del = lam.hat*Kmat.del.zero
    mat.rho = d.rho.kernel.Gauss.pairwise(dMat.list, rho.hat)
    ## score stat
    score.pca2 = t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2
    
    #print(c(score.chi, score.pca))
    
    
    
    #Info matrix
    P0.mat = solve.mat - solve.mat%*%t(xx)%*% solve(xx%*%solve.mat%*%t(xx)) %*% xx%*%solve.mat
    
    
    
    #I0 = matrix(NA, 4+nrow(xx),4+nrow(xx))
    I0 = matrix(NA, 4,4)
    
    I0[1,1] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.del )/2  ##del.del
    
    I0[1,2] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.sig )/2  ##del.sig
    I0[2,1] = I0[1,2]
    
    I0[1,3] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.lam )/2  ##del.lam
    I0[3,1] = I0[1,3]
    
    I0[1,4] = Trace( P0.mat%*%mat.del%*%P0.mat%*%mat.rho )/2  ##del.rho
    I0[4,1] = I0[1,4]
    
    #I0[1, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.del%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 1] = I0[1, 5:(4+nrow(xx))]
    
    I0[2,2] = Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.sig )/2  ##sig.sig
    
    I0[2,3] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.lam )/2  ##sig.lam
    I0[3,2] = I0[2,3]
    
    I0[2,4] =  Trace( P0.mat%*%mat.sig%*%P0.mat%*%mat.rho )/2  ##sig.rho
    I0[4,2] = I0[2,4]
    
    #I0[2, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.sig%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 2] = I0[2, 5:(4+nrow(xx))]
    
    I0[3,3] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.lam )/2  ##lam.lam
    
    I0[3,4] =  Trace( P0.mat%*%mat.lam%*%P0.mat%*%mat.rho )/2  ##lam.rho
    I0[4,3] = I0[3,4]
    
    #I0[3, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.lam%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 3] = I0[3, 5:(4+nrow(xx))]
    
    I0[4,4] =  Trace( P0.mat%*%mat.rho%*%P0.mat%*%mat.rho )/2
    
    # Add jitter to ensure invertible
    if(jitter){
      I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
    }
    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
    
    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
    
    
    #Info for delta
    tot.dim = ncol(I0)
    I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
    #print(I.deldel)
    
    ####
    ## Chi-square approximation
    md = Trace(mat.del%*%P0.mat)/2
    #print(md)
    
    m.pca2 = I.deldel / (2*md)
    d.pca2 = md / m.pca2
    
    pval.pca2 = 1- pchisq(score.pca2/m.pca2, d.pca2)
    
  }
  
  
  
  return(list(test.stat = c(score.chi, score.pca, score.pca2), 
              pvalue = c(pval.chi, pval.pca, pval.pca2), 
              m = c(m.chi, m.pca, m.pca2), 
              d = c(d.chi, d.pca, d.pca2), 
              h.hat=h.hat,
              fit.h0=est.h0))
  
  
}











