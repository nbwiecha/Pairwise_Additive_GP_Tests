#################
##### Kernel machine testing for Gaussian kernel
###################

# Original code by Arnab Maity
# Original source: https://www4.stat.ncsu.edu/~amaity//software/KMTest_garrote_header.r

# Modified by Nate Wiecha

library(PairwiseSims)
obj.loglik.se <- function(theta, Y, D, xx) {
  lambda <- theta[1]
  rho    <- theta[2]
  sigma2 <- theta[3]
  n <- length(Y)
  
  # Pass the single distance matrix directly
  kernel_out <- se_kernel_and_deriv_matrix_cpp(D, rho)
  K <- kernel_out$K
  
  # Build V
  V <- lambda * K + diag(sigma2, n)
  V <- (V + t(V))/2
  
  # Call the existing fast objective C++ function
  val <- obj_loglik_cpp(as.numeric(Y), xx, V)
  return(val)
}

grad.loglik.se <- function(theta, Y, D, xx) {
  lambda <- theta[1]
  rho    <- theta[2]
  sigma2 <- theta[3]
  n <- length(Y)
  
  # Get both K and dK in one fast pass
  kernel_out <- se_kernel_and_deriv_matrix_cpp(D, rho)
  K  <- kernel_out$K
  dK <- kernel_out$dK
  
  # Build V
  V <- lambda * K + diag(sigma2, n)
  V <- (V + t(V))/2

  # Call the existing fast gradient C++ function
  gr <- grad_loglik_cpp(theta, as.numeric(Y), xx, V, K, dK)
  return(as.vector(gr))
}

sqare.diff.mat = function(zz){
  # every column is a sample
  
  Kmat = matrix(NA, nrow=ncol(zz), ncol=ncol(zz))
  for (ii in 1:ncol(zz)){
    Kmat[ii,] = diag(t(zz-zz[,ii])%*%(zz-zz[,ii]))
    
  }
  
  return(Kmat)
}

kernel.Gauss = function(Kmat, scl){
  
  return(exp(-Kmat/scl))
  
}

Kernel.reg.Gauss = function(yy, zz, xx, maxiter=200,
                            starting_val=NULL,
                            sig.min=0.0001, sig.max=10,
                            lam.min=0.0001, lam.max=60,
                            rho.min=0.0001, rho.max=40,
                            criter.use="loglik",
                            zz.p=NULL){
  ### This code does kernel machine regression for gaussian kernel
  ## yy: response nx1
  ## zz: matrix modeled nonparametrically (each col is a sample)
  ## xx: other variables (each col is a sample)
  ##
  ##Estimates: h(z), beta and rho (the kernel parameter)
  
  # NW: change coord descent to bivariate optim wrt rho, lambda to try to
  # reduce number of NA test outputs
  require(fields)
  if(is.null(zz.p)){
    zz.p<-zz
  }
  
  n = length(yy)
  
  ## iterative algo
  dMat = rdist(t(zz))
  dMat.p = rdist(t(zz.p), t(zz))
  dMat.pp = rdist(t(zz.p))
  
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
  
  Kmat = kernel.Gauss(dMat^2, rho.hat)
  
  conv.crit = 0.00001
  criter = 999
  bet.old = bet.hat
  lam.old = lam.hat
  sig.old = sig.hat
  rho.old = rho.hat
  niter = 0
  
  criter.track <- rep(NA, maxiter)
  
  loglik.prev <- obj.loglik.se(c(lam.hat, rho.hat, sig.hat), Y=(yy - t(xx)%*%bet.hat),
                               xx=xx, D=dMat)
  
  pb <- txtProgressBar(max=maxiter, style=3)
  
  while( (criter > conv.crit) & (niter < maxiter) ){
    setTxtProgressBar(pb, niter)
    #print(c(niter, c(bet.old), lam.old, sig.old, rho.old))
    ###########################################
    ###### STEP 1 :fix  sigma, lambda and rho and estimate beta and h().
    ###########################################
    Kmat = kernel.Gauss(dMat^2, rho.hat)
    V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
    iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
    # bet.hat = solve(xx%*%iV%*%t(xx)) %*% xx%*%iV%*%yy
    bet_cov = tryCatch(solve(xx%*%iV%*%t(xx)), error = function(e) MASS::ginv(xx%*%iV%*%t(xx)))
    bet.hat = bet_cov %*% xx%*%iV%*%yy
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
    
    theta.hat <- optim(c(lam.hat, rho.hat, sig.hat), fn=obj.loglik.se, gr=grad.loglik.se,
                       method="L-BFGS-B",
                       xx=xx,  Y=(yy - t(xx)%*%bet.hat),
                       D=dMat,
                       lower=c(lam.min, rho.min, sig.min), upper=c(lam.max, rho.max, sig.max))
    lam.hat <- theta.hat$par[1]
    rho.hat <- theta.hat$par[2]
    sig.hat <- theta.hat$par[3]
    
    # lam.hat = optimize(loglik.lam, interval=c(0.0,60), 
    #                    Y=(yy - t(xx)%*%bet.hat), 
    #                    Kmat=Kmat,
    #                    xx = xx,
    #                    sigmasq = sig.hat)$minimum 
    # 
    # 
    # rho.hat = optimize(loglik.rho, interval=c(rho.min,rho.max), 
    #                    Y=(yy - t(xx)%*%bet.hat), 
    #                    dMat=dMat,
    #                    xx = xx,
    #                    lambda = lam.hat,
    #                    sigmasq = sig.hat)$minimum
    ##
    
    loglik.curr <- obj.loglik.se(c(lam.hat, rho.hat, sig.hat), Y=(yy - t(xx)%*%bet.hat),
                                 xx=xx, D=dMat)
    
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
    
    bet.old = bet.hat
    lam.old = lam.hat
    sig.old = sig.hat
    rho.old = rho.hat
    niter = niter + 1
    
    
  }
  print(niter)
  close(pb)
  ## sanity check
  if(niter >= maxiter){
    print("Did not converge")
    # stop("Did not converge")
  }
  
  Kmat = kernel.Gauss(dMat^2, rho.hat)
  V = sig.hat* diag(1, nrow=n, ncol=n) + lam.hat*Kmat
  iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
  
  np <- ncol(zz.p)
  Kmat.p = kernel.Gauss(dMat.p^2, rho.hat)
  Kmat.pp = kernel.Gauss(dMat.pp^2, rho.hat)
  
  h.hat = lam.hat * Kmat.p %*% iV %*% (yy - t(xx) %*% bet.hat)
  h.sd  = sqrt(diag(
    
    lam.hat*Kmat.pp - lam.hat^2*Kmat.p %*% iV %*% t(Kmat.p)
    
  ))
  
  return(list(h.hat = h.hat,sig=sig.hat, 
              lam=lam.hat, rho=rho.hat, bet=bet.hat, h.sd=h.sd))						
}
# 
# loglik.all <- function(theta, Y, dMat, xx, null){
#   
#   # theta = c(lam, rho)
#   lambda <- theta[1]
#   rho <- theta[2]
#   sigmasq <- theta[3]
#   
#   # if(null){
#   #   kernel.fn <- kernel.Gauss.pairwise.null
#   # }else{
#   #   kernel.fn <- kernel.Gauss.pairwise
#   # }
#   
#   n = length(Y)
#   Kmat = kernel.Gauss(dMat, rho)
#   V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
#   iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
#   
#   out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
#   return(-out/n)
#   
# }
# 
loglik.lam = function(lambda, Y, Kmat, xx, sigmasq){

  n = length(Y)
  V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
  iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))

  out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
  return(-out/n)
}
# 
# 
# loglik.sig = function(sigmasq, Y, Kmat, xx, lambda){
#   
#   n = length(Y)
#   V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
#   iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
#   
#   out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
#   return(-out/n)
# }
# 
# loglik.rho = function(rho, Y, dMat, xx, lambda, sigmasq){
#   
#   n = length(Y)
#   Kmat = kernel.Gauss(dMat, rho)
#   V = sigmasq* diag(1, nrow=n, ncol=n) + lambda*Kmat
#   iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
#   
#   out = -log(det(V))/2 - log( det(xx%*%iV%*%t(xx)) )/2 - t(Y)%*%iV%*%Y/2
#   return(-out/n)
# }

Trace = function(M){sum(diag(M))}





###################### 
##### Testing using garrote kernel
#############################
Kernel.test.Gauss = function(yy, zz, xx, maxiter=200,
                             starting_val=NULL,
                             sig.min=0.0001, sig.max=10,
                             lam.min=0.0001, lam.max=60,
                             rho.min=0.0001, rho.max=40,
                             pca=FALSE, jitter=FALSE){
  # test using kernel machine for the effect of the first row of zz
  
  n = length(yy)
  
  ### estimation under H0
  est.h0 = Kernel.reg.Gauss(yy, zz[-1,], xx, maxiter=maxiter,
                            starting_val=starting_val,
                            lam.min=lam.min, lam.max=lam.max,
                            rho.min=rho.min, rho.max=rho.max,
                            sig.min=sig.min, sig.max=sig.max)
  h.hat = est.h0[[1]]
  sig.hat = est.h0[[2]]
  lam.hat = est.h0[[3]]
  rho.hat = est.h0[[4]]
  bet.hat = est.h0[[5]]
  #print(bet.hat)
  
  dMat = rdist(t(zz[-1,]))
  Kmat = kernel.Gauss(dMat^2, rho.hat)
  
  mat = sig.hat*diag(1, nrow=n, ncol=n) + lam.hat*Kmat
  solve.mat = tryCatch(solve(mat), error = function(e) MASS::ginv(mat))  
  eps.hat = yy - t(xx)%*%bet.hat - h.hat
  
  ##derivative of the garrote kernel
  mul.mat = matrix(1, ncol=1, nrow=n)%*%zz[1,]
  mul.mat = -(mul.mat - t(mul.mat))^2
  Kmat.del = mul.mat * Kmat / rho.hat
  
  ## score stat
  score.chi = lam.hat * t(yy - t(xx)%*%bet.hat) %*% solve.mat %*% Kmat.del %*% solve.mat %*% (yy - t(xx)%*%bet.hat)/2
  
  
  
  ## derivatives of V ('mat' in this code)
  mat.lam = Kmat
  mat.sig = diag(1, nrow=n, ncol=n)
  mat.del = lam.hat*Kmat.del
  mat.rho = lam.hat * dMat^2 * Kmat / (rho.hat^2)
  
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
  
  if(jitter){
    I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
  }
  #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
  #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
  
  #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
  
  
  #Info for delta
  tot.dim = ncol(I0)
  # I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
  inv_I0 <- tryCatch(solve(I0[2:tot.dim,2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim,2:tot.dim]))
  I.deldel =  I0[1,1] - I0[1,2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim,1]
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
    mat.rho = lam.hat * dMat^2 * Kmat / (rho.hat^2)
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
    
    if(jitter){
      I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
    }
    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
    
    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
    
    
    #Info for delta
    tot.dim = ncol(I0)
    # I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
    inv_I0 <- tryCatch(solve(I0[2:tot.dim,2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim,2:tot.dim]))
    I.deldel =  I0[1,1] - I0[1,2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim,1]
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
    mat.rho = lam.hat * dMat^2 * Kmat / (rho.hat^2)
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
    
    if(jitter){
      I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
    }
    #I0[4, 5:(4+nrow(xx))] = c(xx%*%solve.mat%*%mat.rho%*%solve.mat%*%(yy - t(xx)%*%bet.hat))
    #I0[5:(4+nrow(xx)), 4] = I0[4, 5:(4+nrow(xx))]
    
    #I0[5:(4+nrow(xx)), 5:(4+nrow(xx))] = -xx%*%solve.mat%*%t(xx)
    
    
    #Info for delta
    tot.dim = ncol(I0)
    # I.deldel =  I0[1,1] - I0[1,2:tot.dim]%*%solve(I0[2:tot.dim,2:tot.dim])%*%I0[2:tot.dim,1] 
    inv_I0 <- tryCatch(solve(I0[2:tot.dim,2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim,2:tot.dim]))
    I.deldel =  I0[1,1] - I0[1,2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim,1]
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











