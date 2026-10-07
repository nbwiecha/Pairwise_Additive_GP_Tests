#################
##### Kernel machine testing for Gaussian kernel WITH SPATIAL RESIDUAL
###################

# Original code by Arnab Maity
# Original source: https://www4.stat.ncsu.edu/~amaity//software/KMTest_garrote_header.r

# Modified by Nate Wiecha; including addition of spatially-correlated residual

library(PairwiseSims)
library(fields)

obj.loglik.se.sp <- function(theta, Y, D_cov, D_sp, xx) {
  lam_a <- theta[1]
  rho_a <- theta[2]
  lam_sp <- theta[3]
  rho_sp <- theta[4]
  sigma2 <- theta[5]
  n <- length(Y)
  
  kernel_out <- se_kernel_and_deriv_matrix_sp_cpp(D_cov, rho_a, D_sp, rho_sp)
  
  V <- lam_a * kernel_out$K_a + lam_sp * kernel_out$K_sp + diag(sigma2, n)
  V <- (V + t(V))/2
  
  val <- obj_loglik_sp_cpp(as.numeric(Y), xx, V)
  return(val)
}

grad.loglik.se.sp <- function(theta, Y, D_cov, D_sp, xx) {
  lam_a <- theta[1]
  rho_a <- theta[2]
  lam_sp <- theta[3]
  rho_sp <- theta[4]
  sigma2 <- theta[5]
  n <- length(Y)
  
  kernel_out <- se_kernel_and_deriv_matrix_sp_cpp(D_cov, rho_a, D_sp, rho_sp)
  
  V <- lam_a * kernel_out$K_a + lam_sp * kernel_out$K_sp + diag(sigma2, n)
  V <- (V + t(V))/2
  
  gr <- grad_loglik_sp_cpp(theta, as.numeric(Y), xx, V, 
                           kernel_out$K_a, kernel_out$dK_a,
                           kernel_out$K_sp, kernel_out$dK_sp)
  return(as.vector(gr))
}

sqare.diff.mat = function(zz){
  Kmat = matrix(NA, nrow=ncol(zz), ncol=ncol(zz))
  for (ii in 1:ncol(zz)){
    Kmat[ii,] = diag(t(zz-zz[,ii]) %*% (zz-zz[,ii]))
  }
  return(Kmat)
}

kernel.Gauss = function(Kmat, scl){
  return(exp(-Kmat/scl))
}

Trace = function(M){sum(diag(M))}

Kernel.reg.Gauss.sp = function(yy, zz, xx, coords_sp, maxiter=200,
                               starting_val=NULL,
                               sig.min=0.0001, sig.max=10,
                               lam.min=0.0001, lam.max=60,
                               rho.min=0.0001, rho.max=40,
                               criter.use="loglik",
                               zz.p=NULL, coords_sp.p=NULL){
  
  if(is.null(zz.p)){
    zz.p <- zz
    coords_sp.p <- coords_sp
  }
  
  n = length(yy)
  
  dMat = rdist(t(zz))
  dMat.p = rdist(t(zz.p), t(zz))
  dMat.pp = rdist(t(zz.p))
  
  dMat_sp = rdist(coords_sp)
  dMat_sp.p = rdist(coords_sp.p, coords_sp)
  dMat_sp.pp = rdist(coords_sp.p)
  
  bet.hat = matrix(.2, ncol=1, nrow=nrow(xx))
  sig.hat = 1
  rho_a.hat = 1
  lam_a.hat = 0.01
  rho_sp.hat = 1
  lam_sp.hat = 0.01
  
  if(!is.null(starting_val)){
    bet.hat = starting_val$bet
    sig.hat = starting_val$sig
    rho_a.hat = starting_val$rho_a
    lam_a.hat = starting_val$lam_a
    rho_sp.hat = starting_val$rho_sp
    lam_sp.hat = starting_val$lam_sp
  }
  
  Kmat_a = kernel.Gauss(dMat^2, rho_a.hat)
  Kmat_sp = kernel.Gauss(dMat_sp^2, rho_sp.hat)
  
  conv.crit = 0.00001
  criter = 999
  
  bet.old = bet.hat
  lam_a.old = lam_a.hat
  rho_a.old = rho_a.hat
  lam_sp.old = lam_sp.hat
  rho_sp.old = rho_sp.hat
  sig.old = sig.hat
  niter = 0
  
  criter.track <- rep(NA, maxiter)
  
  loglik.prev <- obj.loglik.se.sp(c(lam_a.hat, rho_a.hat, lam_sp.hat, rho_sp.hat, sig.hat), 
                                  Y=(yy - t(xx) %*% bet.hat), xx=xx, D_cov=dMat, D_sp=dMat_sp)
  
  pb <- txtProgressBar(max=maxiter, style=3)
  
  while( (criter > conv.crit) & (niter < maxiter) ){
    setTxtProgressBar(pb, niter)
    
    Kmat_a = kernel.Gauss(dMat^2, rho_a.hat)
    Kmat_sp = kernel.Gauss(dMat_sp^2, rho_sp.hat)
    
    V = sig.hat* diag(1, nrow=n, ncol=n) + lam_a.hat*Kmat_a + lam_sp.hat*Kmat_sp
    iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
    
    bet_cov = tryCatch(solve(xx %*% iV %*% t(xx)), error = function(e) MASS::ginv(xx %*% iV %*% t(xx)))
    bet.hat = bet_cov %*% xx %*% iV %*% yy
    
    h.hat = (lam_a.hat*Kmat_a + lam_sp.hat*Kmat_sp) %*% iV %*% (yy - t(xx) %*% bet.hat)	
    
    theta.hat <- optim(c(lam_a.hat, rho_a.hat, lam_sp.hat, rho_sp.hat, sig.hat), 
                       fn=obj.loglik.se.sp, gr=grad.loglik.se.sp,
                       method="L-BFGS-B",
                       xx=xx, Y=(yy - t(xx) %*% bet.hat),
                       D_cov=dMat, D_sp=dMat_sp,
                       lower=c(lam.min, rho.min, lam.min, rho.min, sig.min), 
                       upper=c(lam.max, rho.max, lam.max, rho.max, sig.max))
    
    lam_a.hat <- theta.hat$par[1]
    rho_a.hat <- theta.hat$par[2]
    lam_sp.hat <- theta.hat$par[3]
    rho_sp.hat <- theta.hat$par[4]
    sig.hat <- theta.hat$par[5]
    
    loglik.curr <- obj.loglik.se.sp(c(lam_a.hat, rho_a.hat, lam_sp.hat, rho_sp.hat, sig.hat), 
                                    Y=(yy - t(xx) %*% bet.hat), xx=xx, D_cov=dMat, D_sp=dMat_sp)
    
    criter1 = (sum(abs(bet.hat - bet.old)) + abs(lam_a.hat - lam_a.old) + abs(rho_a.hat - rho_a.old) + abs(lam_sp.hat - lam_sp.old) + abs(rho_sp.hat - rho_sp.old) + abs(sig.hat - sig.old))
    criter2 <- abs(loglik.curr - loglik.prev)
    criter3 = sum(abs(lam_a.hat - lam_a.old) + abs(rho_a.hat - rho_a.old) + abs(lam_sp.hat - lam_sp.old) + abs(rho_sp.hat - rho_sp.old) + abs(sig.hat - sig.old))
    
    loglik.prev <- loglik.curr
    criter <- criter1
    if(criter.use=="loglik"){ criter <- criter2 }
    if(criter.use=="theta"){ criter <- criter3 }
    
    criter.track[niter] <- criter
    
    bet.old = bet.hat
    lam_a.old = lam_a.hat
    rho_a.old = rho_a.hat
    lam_sp.old = lam_sp.hat
    rho_sp.old = rho_sp.hat
    sig.old = sig.hat
    niter = niter + 1
  }
  
  print(niter)
  close(pb)
  if(niter >= maxiter){ print("Did not converge") }
  
  Kmat_a = kernel.Gauss(dMat^2, rho_a.hat)
  Kmat_sp = kernel.Gauss(dMat_sp^2, rho_sp.hat)
  
  V = sig.hat* diag(1, nrow=n, ncol=n) + lam_a.hat*Kmat_a + lam_sp.hat*Kmat_sp
  iV = tryCatch(solve(V), error = function(e) MASS::ginv(V))
  
  np <- ncol(zz.p)
  Kmat_a.p = kernel.Gauss(dMat.p^2, rho_a.hat)
  Kmat_a.pp = kernel.Gauss(dMat.pp^2, rho_a.hat)
  Kmat_sp.p = kernel.Gauss(dMat_sp.p^2, rho_sp.hat)
  Kmat_sp.pp = kernel.Gauss(dMat_sp.pp^2, rho_sp.hat)
  
  K_cross = lam_a.hat*Kmat_a.p + lam_sp.hat*Kmat_sp.p
  K_full_pp = lam_a.hat*Kmat_a.pp + lam_sp.hat*Kmat_sp.pp
  
  h.hat = K_cross %*% iV %*% (yy - t(xx) %*% bet.hat)
  h.sd  = sqrt(diag(K_full_pp - K_cross %*% iV %*% t(K_cross)))
  
  return(list(h.hat = h.hat, sig = sig.hat, 
              lam_a = lam_a.hat, rho_a = rho_a.hat, 
              lam_sp = lam_sp.hat, rho_sp = rho_sp.hat, 
              bet = bet.hat, h.sd = h.sd))						
}

Kernel.test.Gauss.sp = function(yy, zz, xx, coords_sp, maxiter=200,
                                starting_val=NULL,
                                sig.min=0.0001, sig.max=10,
                                lam.min=0.0001, lam.max=60,
                                rho.min=0.0001, rho.max=40,
                                pca=FALSE, jitter=FALSE){
  
  n = length(yy)
  
  est.h0 = Kernel.reg.Gauss.sp(yy, zz[-1, , drop=FALSE], xx, coords_sp=coords_sp, maxiter=maxiter,
                               starting_val=starting_val,
                               lam.min=lam.min, lam.max=lam.max,
                               rho.min=rho.min, rho.max=rho.max,
                               sig.min=sig.min, sig.max=sig.max)
  
  h.hat = est.h0$h.hat
  sig.hat = est.h0$sig
  lam_a.hat = est.h0$lam_a
  rho_a.hat = est.h0$rho_a
  lam_sp.hat = est.h0$lam_sp
  rho_sp.hat = est.h0$rho_sp
  bet.hat = est.h0$bet
  
  dMat = rdist(t(zz[-1, , drop=FALSE]))
  Kmat_a = kernel.Gauss(dMat^2, rho_a.hat)
  
  dMat_sp = rdist(coords_sp)
  Kmat_sp = kernel.Gauss(dMat_sp^2, rho_sp.hat)
  
  mat = sig.hat*diag(1, nrow=n, ncol=n) + lam_a.hat*Kmat_a + lam_sp.hat*Kmat_sp
  solve.mat = tryCatch(solve(mat), error = function(e) MASS::ginv(mat))  
  
  mul.mat = matrix(1, ncol=1, nrow=n) %*% zz[1,]
  mul.mat = -(mul.mat - t(mul.mat))^2
  Kmat.del = mul.mat * Kmat_a / rho_a.hat
  
  score.chi = lam_a.hat * t(yy - t(xx) %*% bet.hat) %*% solve.mat %*% Kmat.del %*% solve.mat %*% (yy - t(xx) %*% bet.hat) / 2
  
  mat.lam_a = Kmat_a
  mat.sig = diag(1, nrow=n, ncol=n)
  mat.del = lam_a.hat*Kmat.del
  mat.rho_a = lam_a.hat * dMat^2 * Kmat_a / (rho_a.hat^2)
  mat.lam_sp = Kmat_sp
  mat.rho_sp = lam_sp.hat * dMat_sp^2 * Kmat_sp / (rho_sp.hat^2)
  
  P0.mat = solve.mat - solve.mat %*% t(xx) %*% solve(xx %*% solve.mat %*% t(xx)) %*% xx %*% solve.mat
  
  I0 = matrix(NA, 6, 6)
  
  I0[1,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.del )/2 
  
  I0[1,2] = I0[2,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.sig )/2 
  I0[1,3] = I0[3,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_a )/2 
  I0[1,4] = I0[4,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_a )/2 
  I0[1,5] = I0[5,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_sp )/2 
  I0[1,6] = I0[6,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_sp )/2 
  
  I0[2,2] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.sig )/2 
  I0[2,3] = I0[3,2] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.lam_a )/2 
  I0[2,4] = I0[4,2] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.rho_a )/2 
  I0[2,5] = I0[5,2] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.lam_sp )/2 
  I0[2,6] = I0[6,2] = Trace( P0.mat %*% mat.sig %*% P0.mat %*% mat.rho_sp )/2 
  
  I0[3,3] = Trace( P0.mat %*% mat.lam_a %*% P0.mat %*% mat.lam_a )/2 
  I0[3,4] = I0[4,3] = Trace( P0.mat %*% mat.lam_a %*% P0.mat %*% mat.rho_a )/2 
  I0[3,5] = I0[5,3] = Trace( P0.mat %*% mat.lam_a %*% P0.mat %*% mat.lam_sp )/2 
  I0[3,6] = I0[6,3] = Trace( P0.mat %*% mat.lam_a %*% P0.mat %*% mat.rho_sp )/2 
  
  I0[4,4] = Trace( P0.mat %*% mat.rho_a %*% P0.mat %*% mat.rho_a )/2
  I0[4,5] = I0[5,4] = Trace( P0.mat %*% mat.rho_a %*% P0.mat %*% mat.lam_sp )/2 
  I0[4,6] = I0[6,4] = Trace( P0.mat %*% mat.rho_a %*% P0.mat %*% mat.rho_sp )/2 
  
  I0[5,5] = Trace( P0.mat %*% mat.lam_sp %*% P0.mat %*% mat.lam_sp )/2 
  I0[5,6] = I0[6,5] = Trace( P0.mat %*% mat.lam_sp %*% P0.mat %*% mat.rho_sp )/2 
  
  I0[6,6] = Trace( P0.mat %*% mat.rho_sp %*% P0.mat %*% mat.rho_sp )/2 
  
  if(jitter){
    I0 <- I0 + diag(10^(-8), nrow=nrow(I0), ncol=ncol(I0))
  }
  
  tot.dim = ncol(I0)
  inv_I0 <- tryCatch(solve(I0[2:tot.dim, 2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim, 2:tot.dim]))
  I.deldel =  I0[1,1] - I0[1, 2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim, 1]
  
  md = lam_a.hat*Trace(Kmat.del %*% P0.mat)/2
  
  m.chi = I.deldel / (2*md)
  d.chi = md / m.chi
  
  pval.chi = 1 - pchisq(score.chi/m.chi, d.chi)
  
  score.pca <- score.pca2 <- pval.pca <- pval.pca2 <- m.pca <- m.pca2 <- d.pca <- d.pca2 <- NA
  
  if(pca){
    svd.out = eigen(Kmat.del) 
    Kmat.del.plus = svd.out$vectors %*% diag(abs(svd.out$values), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    
    mat.del = lam_a.hat*Kmat.del.plus
    score.pca = t(yy - t(xx) %*% bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx) %*% bet.hat)/2
    
    I0[1,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.del )/2 
    I0[1,2] = I0[2,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.sig )/2 
    I0[1,3] = I0[3,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_a )/2 
    I0[1,4] = I0[4,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_a )/2 
    I0[1,5] = I0[5,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_sp )/2 
    I0[1,6] = I0[6,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_sp )/2 
    
    inv_I0 <- tryCatch(solve(I0[2:tot.dim, 2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim, 2:tot.dim]))
    I.deldel =  I0[1,1] - I0[1, 2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim, 1]
    
    md = Trace(mat.del %*% P0.mat)/2
    m.pca = I.deldel / (2*md)
    d.pca = md / m.pca
    pval.pca = 1- pchisq(score.pca/m.pca, d.pca)
    
    Kmat.del.zero = svd.out$vectors %*% diag(svd.out$values*(svd.out$values>0), ncol=length(svd.out$values), nrow=length(svd.out$values)) %*% t(svd.out$vectors)
    
    mat.del = lam_a.hat*Kmat.del.zero
    score.pca2 = t(yy - t(xx) %*% bet.hat) %*% solve.mat %*% mat.del %*% solve.mat %*% (yy - t(xx) %*% bet.hat)/2
    
    I0[1,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.del )/2 
    I0[1,2] = I0[2,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.sig )/2 
    I0[1,3] = I0[3,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_a )/2 
    I0[1,4] = I0[4,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_a )/2 
    I0[1,5] = I0[5,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.lam_sp )/2 
    I0[1,6] = I0[6,1] = Trace( P0.mat %*% mat.del %*% P0.mat %*% mat.rho_sp )/2 
    
    inv_I0 <- tryCatch(solve(I0[2:tot.dim, 2:tot.dim]), error = function(e) MASS::ginv(I0[2:tot.dim, 2:tot.dim]))
    I.deldel =  I0[1,1] - I0[1, 2:tot.dim] %*% inv_I0 %*% I0[2:tot.dim, 1]
    
    md = Trace(mat.del %*% P0.mat)/2
    m.pca2 = I.deldel / (2*md)
    d.pca2 = md / m.pca2
    pval.pca2 = 1 - pchisq(score.pca2/m.pca2, d.pca2)
  }
  
  return(list(test.stat = c(score.chi, score.pca, score.pca2), 
              pvalue = c(pval.chi, pval.pca, pval.pca2), 
              m = c(m.chi, m.pca, m.pca2), 
              d = c(d.chi, d.pca, d.pca2), 
              h.hat=h.hat,
              fit.h0=est.h0))
}