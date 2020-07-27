library(MASS)
library(Matrix)
## Generate simulated dataset
# For fixed X, use the same seed.x (default)
# For random X, specify different seed.x for different replications
gen.beta.Sigma <- function(p, rho, type){
  # Covariance matrix and true beta
  # beta0 = 0 # for the intercept term
  if(type == "Sparse-Ex1" | type == "Omit"){
    # Zhang, Li and Tsai (2010)
    Sigma = Matrix(rho^abs(outer(1:p, 1:p, "-")), sparse = TRUE)
    beta = Matrix(c(rev(seq(1, 3.5, by=0.5)), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "Sparse-Ex2"){
    Sigma = bdiag(list(rho^abs(outer(1:6, 1:6, "-")), rho^abs(outer(7:p, 7:p, "-"))))
    beta = Matrix(c(c(1,1,3,3,5,5), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "Sparse-Ex3"){
    Sigma = Diagonal(p)
    for(i in 1:6){
      Sigma[i,i+6] = Sigma[i+6,i] = rho
    }
    beta = Matrix(c(rep(1,6), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "Sparse-Ex4"){
    Sigma = Diagonal(p)
    for(i in 1:3){
      Sigma[2*i-1,2*i] = Sigma[2*i,2*i-1] = rho
    }
    beta = Matrix(c(c(1,-1,3,-3,5,-5), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "Dense-Ex1"){
    Sigma = Matrix(rho^abs(outer(1:p, 1:p, "-")), sparse=TRUE)
    kappa = 10
    beta = Matrix((-1)^seq(1,p) * exp(-seq(1,p)/kappa), ncol=1, sparse = TRUE)
  }else if(type == "Dense-Ex2"){
    Sigma = Matrix(rho^abs(outer(1:p, 1:p, "-")), sparse=TRUE)
    kappa = 50
    beta = Matrix((-1)^seq(1,p) * exp(-seq(1,p)/kappa), ncol=1, sparse = TRUE)
  }else if(type == "Exponential"){
    Sigma = Diagonal(p)
    #beta = beta0 = NULL
    beta = NULL
  }
  
  if(rho != 0){
    obj = svd(Sigma)
    Sigma.half = Matrix(obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v), sparse=TRUE)
  }else{
    Sigma.half = NULL
  }
  
  return(list(beta=beta, Sigma=Sigma, Sigma.half=Sigma.half))
}
gen.beta.Sigma.restriction <- function(p, rho, type){
  # Covariance matrix and true beta
  # Zhang, Li and Tsai (2010)
  Sigma = Matrix(rho^abs(outer(1:p, 1:p, "-")), sparse = TRUE)
  if(rho != 0){
    obj = svd(Sigma)
    Sigma.half = Matrix(obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v), sparse=TRUE)
  }else{
    Sigma.half = NULL
  }
  
  if(type == "Ex1"){
    beta = Matrix(c(2,2,2,1,1,1), ncol=1)
    # ## True restriction
    # R = Matrix(0, nrow=4, ncol=p)
    # R[1,1] = R[2,2] = R[3,4] = R[4,5] = 1
    # R[1,2] = R[2,3] = R[3,5] = R[4,6] = -1
    # r = rep(0, 4)
    # ## All candidate restrictions
    # m = nrow(R)
    
    # R_candidate = rbind(R, c(0,0,1,-1,0,0))
    # r_candidate = c(r, 0)
    # 
    # R_all = list(Matrix(, nrow=0, ncol=p)) # No restrictions
    # r_all = list(numeric(0)) # No restrictions
    # R_all = c(R_all, lapply(1:nrow(R_candidate), function(i){R_candidate[1:i, , drop=FALSE]}))
    # r_all = c(r_all, lapply(1:length(r_candidate), function(i){r_candidate[1:i]}))
    # R_all = c(R_all, list(Diagonal(p))) # null model
    # r_all = c(r_all, list(rep(0,p))) # null model
    # 
    # R_all = rev(R_all)
    # r_all = rev(r_all)
    
    # beta = Matrix(c(2,2,2,1,1,1,rep(0,p-6)), ncol=1)
    # R_candidate = Matrix(0, nrow=p, ncol=p)
    # R_candidate[1,1] = R_candidate[2,2] = R_candidate[3,4] = R_candidate[4,5] = 1
    # R_candidate[1,2] = R_candidate[2,3] = R_candidate[3,5] = R_candidate[4,6] = -1
    # for(i in 7:p){
    #   R_candidate[i-2,i] = 1
    # }
    # R_candidate[p-1, 1] = R_candidate[p, 4] = 1
    # r_candidate = rep(0, p)
    # 
    # 
    # R_all = list(Matrix(, nrow=0, ncol=p)) # No restrictions
    # r_all = list(numeric(0)) # No restrictions
    # R_all = c(R_all, lapply(1:nrow(R_candidate), function(i){R_candidate[1:i, , drop=FALSE]}))
    # r_all = c(r_all, lapply(1:length(r_candidate), function(i){r_candidate[1:i]}))

    R_correct = Matrix(rbind(c(1,-1,0,0,0,0), c(0,1,-1,0,0,0), c(0,0,0,1,-1,0), c(0,0,0,0,1,-1), c(1,0,0,-2,0,0)), sparse=TRUE)
    r_correct = rep(0, 5)
    R_wrong = Matrix(rbind(c(1,0,0,-1,0,0), c(0,1,0,0,-1,0), c(0,0,1,0,0,-1), c(1,-2,0,0,0,0), c(0,1,-2,0,0,0)), sparse=TRUE)
    r_wrong = rep(0, 5)
    
    # R_all = lapply(1:nrow(R), function(i){R[1:i, , drop=FALSE]})
    # r_all = lapply(1:length(r), function(i){r[1:i]})
    # R_all = c(R_all, list(rbind(R, c(1,0,0,-2,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # 
    # R_wrong = Matrix(0, nrow=5, ncol=p)
    # R_wrong[1,1] = R_wrong[2,2] = R_wrong[3,3] = R_wrong[4,1] = R_wrong[5,2] = 1
    # R_wrong[1,4] = R_wrong[2,5] = R_wrong[3,6] = -1
    # R_wrong[4,2] = R_wrong[5,3] = -2
    # r_wrong = rep(0, 5)
    # 
    # R_all = c(R_all, lapply(1:nrow(R_wrong), function(i){R_wrong[1:i, , drop=FALSE]}))
    # r_all = c(r_all, lapply(1:length(r_wrong), function(i){r_wrong[1:i]}))
    # 
    # R_all = c(R_all, list(Matrix(, nrow=0, ncol=p)))
    # r_all = c(r_all, list(numeric(0)))
    # R_all = c(R_all, list(Diagonal(p)))
    # r_all = c(r_all, list(rep(0,p)))

    # # Under-specified restrictions
    # R_all = list(matrix(, nrow=0, ncol=p)) # No restrictions
    # r_all = list(numeric(0)) # No restrictions
    # R_all = c(R_all, list(matrix(c(1,-1,0,0,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(rbind(c(1,-1,0,0,0,0), c(0,0,0,1,-1,0), c(0,0,0,0,1,-1))))
    # r_all = c(r_all, list(c(0,0,0)))
    # # for(i in 1:m){
    # #   tmp = combn(m,i)
    # #   R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R[ii, , drop=FALSE]}))
    # #   r_all = c(r_all, replicate(dim(tmp)[2], rep(0, i), simplify = FALSE))
    # # }
    # 
    # # Over-specified restrictions
    # R_all = c(R_all, list(R))
    # r_all = c(r_all, list(r))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,-2,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,1,0,0))))
    # r_all = c(r_all, list(c(r, 3)))
    # 
    # # Under-specified and wrong
    # R_all = c(R_all, list(matrix(c(1,0,0,-1,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(rbind(c(1,0,0,-1,0,0), c(0,1,0,0,-1,0))))
    # r_all = c(r_all, list(c(0,0)))
    # R_all = c(R_all, list(rbind(c(1,0,0,-1,0,0), c(0,1,0,0,-1,0), c(0,0,1,0,0,-1))))
    # r_all = c(r_all, list(c(0,0,0)))
    # 
    # # Over-specified and wrong
    # R_all = c(R_all, list(rbind(R, c(1,0,0,-1,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,-3,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,1,0,0))))
    # r_all = c(r_all, list(c(r, 1)))

  }else if(type == "Ex2"){
    beta = Matrix(c(-2,2,2,-1.5,-1.5,1), ncol=1)
    # ## True restriction
    # R = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1)), sparse=TRUE)
    # r = rep(0, 2)
    # ## All candidate restrictions
    # m = nrow(R)
    # 
    # R_candidate = rbind(R, c(0,0,1,1,0,0), c(0,1,1,0,0,0), c(0,0,0,1,1,0))
    # r_candidate = c(r, 0, 0, 0)
    # 
    # R_all = list(Matrix(, nrow=0, ncol=p)) # No restrictions
    # r_all = list(numeric(0)) # No restrictions
    # R_all = c(R_all, lapply(1:nrow(R_candidate), function(i){R_candidate[1:i, , drop=FALSE]}))
    # r_all = c(r_all, lapply(1:length(r_candidate), function(i){r_candidate[1:i]}))
    # R_all = c(R_all, list(Diagonal(p))) # null model
    # r_all = c(r_all, list(rep(0,p))) # null model
    # 
    # R_all = rev(R_all)
    # r_all = rev(r_all)
    
    R_correct = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), c(1,0,1,0,0,0), c(0,0,1,0,0,-2), c(0,0,0,1,-1,0)), sparse=TRUE)
    r_correct = rep(0, 5)
    R_wrong = Matrix(rbind(c(0,0,0,0,1,1), c(1,1,1,1,0,0), c(0,0,0,1,0,1), c(-2,0,0,1,0,0), c(1,-1,0,0,0,0)), sparse=TRUE)
    r_wrong = rep(0, 5)
       
    # # No restrictions
    # R_all = list(matrix(, nrow=0, ncol=p))
    # r_all = list(numeric(0))
    # # Under-specified restrictions
    # R_all = c(R_all, list(matrix(c(1,1,0,0,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(matrix(c(0,0,1,1,1,1), nrow=1)))
    # r_all = c(r_all, list(0))
    # # for(i in 1:m){
    # #   tmp = combn(m,i)
    # #   R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R[ii, , drop=FALSE]}))
    # #   r_all = c(r_all, replicate(dim(tmp)[2], rep(0, i), simplify = FALSE))
    # # }
    # 
    # # Over-specified restrictions
    # R_all = c(R_all, list(R))
    # r_all = c(r_all, list(r))
    # R_all = c(R_all, list(rbind(R, c(1,0,-1,0,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(0,0,0,1,-1,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # 
    # # Under-specified and wrong
    # R_all = c(R_all, list(matrix(c(1,-1,0,0,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(matrix(c(1,0,1,0,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(matrix(c(1,1,1,1,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # 
    # # Over-specified and wrong
    # R_all = c(R_all, list(rbind(R, c(1,0,-2,0,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(0,0,0,1,-2,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(0,0,1,1,0,0))))
    # r_all = c(r_all, list(c(r, 0)))

  }else if(type == "Ex3"){
    beta = Matrix(c(2,-2,1,-1,0.5,-0.5), ncol=1)
    
    # ## True restriction
    # R = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,0,0), c(0,0,0,0,1,1)), sparse=TRUE)
    # r = rep(0, 3)
    # ## All candidate restrictions
    # m = nrow(R)
    # 
    # R_candidate = rbind(R, c(0,1,1,0,0,0), c(0,0,0,1,1,0))
    # r_candidate = c(r, 0, 0)
    # 
    # R_all = list(Matrix(, nrow=0, ncol=p)) # No restrictions
    # r_all = list(numeric(0)) # No restrictions
    # R_all = c(R_all, lapply(1:nrow(R_candidate), function(i){R_candidate[1:i, , drop=FALSE]}))
    # r_all = c(r_all, lapply(1:length(r_candidate), function(i){r_candidate[1:i]}))
    # R_all = c(R_all, list(Diagonal(p))) # null model
    # r_all = c(r_all, list(rep(0,p))) # null model
    # 
    # R_all = rev(R_all)
    # r_all = rev(r_all)
    
    R_correct = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,0,0), c(0,0,0,0,1,1), c(1,0,-2,0,0,0), c(0,0,1,0,-2,0)), sparse=TRUE)
    r_correct = rep(0, 5)
    R_wrong = Matrix(rbind(c(1,0,0,1,0,0), c(0,1,0,0,1,0), c(0,0,1,0,0,1), c(1,-2,0,0,0,0), c(0,1,-2,0,0,0)), sparse=TRUE)
    r_wrong = rep(0, 5)
    
    # # No restrictions
    # R_all = list(matrix(, nrow=0, ncol=p))
    # r_all = list(numeric(0))
    # # Under-specified restrictions
    # # for(i in 1:m){
    # #   tmp = combn(m,i)
    # #   R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R[ii, , drop=FALSE]}))
    # #   r_all = c(r_all, replicate(dim(tmp)[2], rep(0, i), simplify = FALSE))
    # # }
    # R_all = c(R_all, list(matrix(c(1,1,0,0,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(rbind(c(1,1,0,0,0,0), c(0,0,1,1,0,0))))
    # r_all = c(r_all, list(c(0,0)))
    # 
    # # Over-specified restrictions
    # R_all = c(R_all, list(R))
    # r_all = c(r_all, list(r))
    # R_all = c(R_all, list(rbind(R, c(1,0,-2,0,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,-2,0,0,0), c(0,0,1,0,-2,0))))
    # r_all = c(r_all, list(c(r, 0, 0)))
    # 
    # # Under-specified and wrong
    # R_all = c(R_all, list(matrix(c(1,0,0,-1,0,0), nrow=1)))
    # r_all = c(r_all, list(0))
    # R_all = c(R_all, list(rbind(c(1,0,0,-1,0,0), c(0,1,0,0,-1,0))))
    # r_all = c(r_all, list(c(0,0)))
    # R_all = c(R_all, list(rbind(c(1,0,0,-1,0,0), c(0,0,1,0,0,-1))))
    # r_all = c(r_all, list(c(0,0)))
    # 
    # # Over-specified and wrong
    # R_all = c(R_all, list(rbind(R, c(1,0,-3,0,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,-2,0,0))))
    # r_all = c(r_all, list(c(r, 0)))
    # R_all = c(R_all, list(rbind(R, c(1,0,0,0,-2,0))))
    # r_all = c(r_all, list(c(r, 0)))
    
  }
  R_all = r_all = list()
  for(i in 1:5){
    tmp = combn(5,i)
    R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R_correct[ii, , drop=FALSE]}))
    r_all = c(r_all, lapply(split(tmp, col(tmp)), function(ii){r_correct[ii]}))
    R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R_wrong[ii, , drop=FALSE]}))
    r_all = c(r_all, lapply(split(tmp, col(tmp)), function(ii){r_wrong[ii]}))
  }
  R_all = c(R_all, list(Matrix(, nrow=0, ncol=p)))
  r_all = c(r_all, list(numeric(0)))
  R_all = c(R_all, list(Diagonal(p)))
  r_all = c(r_all, list(rep(0,p)))
  return(list(beta=beta, Sigma=Sigma, Sigma.half=Sigma.half, restriction_candidate=list(R=R_all, r=r_all)))
}
gen.data <- function(n, p, snr, type, beta, Sigma.half, seed.x, seed.y){
  # generate design matrix X
  # for all cases except Exponential model, the X is random
  # for the Exponential model, X is fixed
  set.seed(seed.x)
  # x = mvrnorm(n, mu=rep(0,p), Sigma=Sigma)
  x = matrix( rnorm(n*p), nrow=n, ncol=p )
  if(!is.null(Sigma.half)){
    x = x %*% Sigma.half
  }
  
  # if(type == "Exponential"){
  #   mu = Matrix(exp(4*seq(1,n) / n), ncol=1)
  # }else{
  mu = x %*% beta
  # }
  
  # # Exponential model
  # # Trigonometric design matrix, X has orthonormal columns
  # tmp = 2*pi*outer(1:n, 1:(p/2), "*") / n
  # x1 = sin(tmp)
  # x2 = cos(tmp)
  # x = matrix(NA, nrow=n, ncol=p)
  # x[, seq(1,p,by=2)] = x1
  # x[, seq(2,p,by=2)] = x2
  # x = replicate(nrep, x, simplify = FALSE)
  # beta  = beta0 = NULL
  # mu = exp(4*seq(1,n) / n)
  # mu = replicate(nrep, mu, simplify = FALSE)
  
  # sigma = mean( unlist( lapply(mu, function(xx){ sqrt(var(xx) / snr) }) ) )
  sigma = as.numeric(sqrt(var(as.vector(mu)) / snr))
  set.seed(seed.y)
  # generate the response
  y = mu + Matrix(rnorm(n, mean=0, sd=sigma), ncol=1)
  
  # # Oracle R squares (regression y upon true predictors)
  # if(print.r2){
  #   r.squares = 0
  #   for(rep in 1:nrep){
  #     if(type != "Exponential"){
  #       r.squares = r.squares + summary(lm(y[[rep]]~x[[rep]][,beta!=0]))$r.squared
  #     }else{
  #       r.squares = r.squares + summary(lm(y[[rep]]~x[[rep]]))$r.squared
  #     }
  #   }
  #   r.squares = r.squares/nrep
  #   print(paste("Average R squares of regressing upon true predictors is ", as.character(round(r.squares,2)),sep=""))
  # }
  
  # Treat the last true predictor as missing, for the omit predictor case
  # if(type == "Omit"){
  #   x = x[, -6]
  #   # beta = beta0 = NULL
  #   beta = NULL
  # }
  
  # return(list(x=x, y=y, beta=beta, beta0=beta0, mu=mu, sigma=sigma))
  return(list(x=x, y=y, mu=mu, sigma=sigma))
}

## Coefficient vectors of the nested models
coef.nest <- function(x, y, intercept=TRUE){
  n = dim(x)[1]
  p = dim(x)[2]
  maxstep = ifelse(intercept, min(n, p+1)-1, min(n, p))
  x = x[, 1:maxstep, drop=FALSE]
  
  # standardization
  if (intercept) {
    mean_x = colMeans(x)
    mean_y = mean(y)
    x = scale(x, center = mean_x, scale = FALSE)
    y = scale(y, center = mean_y, scale = FALSE)
  }else {
    mean_x = rep(0, maxstep)
    mean_y = 0
  }
  sd_demeanedx = sqrt(colSums(x^2))
  x = scale(x, center = FALSE, scale = sd_demeanedx)
  # fit LS on nested models via QR decomposition
  qr_decomp = qr(x)
  Q = qr.Q(qr_decomp)
  R = qr.R(qr_decomp)
  z = t(Q) %*% y
  
  # transform coef in Q space back to X space
  trans.q.to.x <- function(beta.q){
    beta.x = diag(1/sd_demeanedx) %*% backsolve(R, beta.q)
    beta.x = cbind(0, beta.x)
    if(intercept){
      beta.x = rbind((mean_y - mean_x %*% beta.x), beta.x)
    }
    return(beta.x)
  }
  
  beta_q = matrix(rep(z, maxstep), nrow=maxstep, byrow=FALSE)
  beta_q = beta_q * upper.tri(beta_q, diag=TRUE)
  beta_nest = trans.q.to.x(beta_q)
  
  if(maxstep < p){
    beta_nest = rbind(beta_nest, Matrix(0, nrow=p-maxstep, ncol=maxstep+1))
  }
  
  return(Matrix(beta_nest, sparse = TRUE))
}

## Restricted coefficient vectors
## Given a set of linear restrictions, returns a matrix of coefficient estimates
## (R, r) here are lists, they are also in pairs
coef.restrict <- function(x, y, R, r){
  # n = dim(x)[1]
  # p = dim(x)[2]
  
  ind_null = which( unlist(lapply(R, function(RR){dim(RR)[1]}))==0 )
  # if(length(ind_null) > 1){
  #   stop("only one NULL entry of R is allowed for now")
  # }
  # if(length(ind_null) > 0){
  #   R[ind_null] = NULL
  #   r[ind_null] = NULL
  # }
  xtx_inv = solve(t(x) %*% x)
  # OLS estimator on the full set of data
  betahat_f = xtx_inv %*% t(x) %*% y
  tmp_function <- function(RR, rr){
    betahat_f + xtx_inv %*% t(RR) %*% solve( RR %*% xtx_inv %*% t(RR) ) %*% (rr - RR %*% betahat_f)
  }
  betahat = do.call(cbind, Map(tmp_function, R[-ind_null], r[-ind_null]))
  # if(length(ind_null) > 0){
  #   betahat = append(betahat, list(betahat_f), after = ind_null-1)
  # }
  betahat_final = Matrix(0, nrow=p, ncol=length(R))
  betahat_final[,ind_null] = matrix(rep(betahat_f, length(ind_null)), nrow=p, byrow=FALSE)
  betahat_final[,-ind_null] = betahat
  return(Matrix(betahat_final, sparse = TRUE))
}

## Cross validation to choose the best candidate from the nested models
cv.nest <- function(x, y, n.folds=10, intercept=FALSE){
  n = dim(x)[1]
  p = dim(x)[2]

  if(n.folds == n){
    ## Use PRESS statistic for LOO CV
    betahat = coef.nest(x, y, intercept=intercept)
    if(!intercept){
      muhat = x %*% betahat
    }else{
      muhat = cbind(1, x) %*% betahat
    }
    qr_decomp = qr(x)
    Q = qr.Q(qr_decomp)
    H = apply(Q^2, 1, cumsum)
    H = rbind(rep(0, n), H)
    calc.press <- function(i){
      ( (y-muhat[,i]) / (1 - H[i,]) )^2
    }
    cv_tmp = lapply(1:(p+1), calc.press)
    cv = colMeans(do.call(cbind, cv_tmp))
  }else{
    fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
    cv_tmp = list()
    for(fold in 1:n.folds){
      # Split the training and testing sets
      test.index = which(fold.index==fold)
      x.test = x[test.index, , drop=FALSE]
      y.test = y[test.index]
      x.train = x[-test.index, , drop=FALSE]
      y.train = y[-test.index]
      # Fit the model on training set
      betahat = coef.nest(x.train, y.train, intercept=intercept)
      
      # Evaluate on the testing set
      if(!intercept){
        cv_tmp[[fold]] = colSums( sweep(x.test%*%betahat, 1, y.test, "-")^2 )
      }else{
        cv_tmp[[fold]] = colSums( sweep(cbind(rep(1,nrow(x.test)),x.test)%*%betahat, 1, y.test, "-")^2 )
      }
    }
    cv = colMeans(do.call(rbind, cv_tmp))
    # Fit on the full sample
    betahat = coef.nest(x, y, intercept=intercept)
  }
  return(list(betahat=betahat, i.min=which.min(cv)))
}

## Cross validation to choose the best candidate from a sequence of restricted least squares models
cv.restrict <- function(x, y, R, r, n.folds=10){
  n = dim(x)[1]
  p = dim(x)[2]
  
  if(n.folds == n){
    ind_null = which( unlist(lapply(R, function(RR){dim(RR)[1]}))==0 )
    ## Use PRESS statistic for LOO-CV
    ## See Tarpey (2000) for details
    xtx_inv = solve(t(x) %*% x)
    H_base = xtx_inv %*% t(x)
    betahat_f = H_base %*% y
    H = x %*% H_base
    calc.H <- function(RR){
      HRQ_base = xtx_inv %*% t(RR) %*% solve(RR %*% xtx_inv %*% t(RR)) 
      H_R = HRQ_base %*% RR
      H_Q = x %*% H_R %*% xtx_inv %*% t(x)
      return(list(HRQ_base=HRQ_base, H_R=H_R, H_Q=H_Q))
    }
    H_family = lapply(R[-ind_null], calc.H)
    calc.betahat <- function(HH_family, rr){
      (diag(p) - HH_family$H_R) %*% betahat_f + HH_family$HRQ_base %*% rr
    }
    betahat = Map(calc.betahat, H_family, r[-ind_null])
    muhat = lapply(betahat, function(bb){x %*% bb})
    calc.press <- function(HH_family, mm){
      ( (y - mm) / (1 - diag(H-HH_family$H_Q)) )^2
    }
    cv_tmp = do.call(cbind, Map(calc.press, H_family, muhat))
    
    cv_tmp_final = Matrix(0, nrow=n, ncol=length(R))
    cv_tmp_final[,ind_null] = Matrix(rep(( (y - x %*% betahat_f) / (1 - diag(H)) )^2, length(ind_null)), nrow=n, byrow=FALSE)
    cv_tmp_final[,-ind_null] = cv_tmp
    cv = colMeans(cv_tmp_final)
    
    betahat_final = Matrix(0, nrow=p, ncol=length(R))
    betahat_final[,ind_null] = matrix(rep(betahat_f, length(ind_null)), nrow=p, byrow=FALSE)
    betahat_final[,-ind_null] = do.call(cbind, betahat)
    betahat = Matrix(betahat_final, sparse = TRUE)

  }else{
    ## Numerically calculate the CV error for folds other than n
    cv_tmp = list()
    fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
    
    for(fold in 1:n.folds){
      # Split the training and testing sets
      test.index = which(fold.index==fold)
      x.test = x[test.index, , drop=FALSE]
      y.test = y[test.index]
      x.train = x[-test.index, , drop=FALSE]
      y.train = y[-test.index]
      # Fit the model on training set
      betahat = coef.restrict(x.train, y.train, R, r)
      
      # Evaluate on the testing set
      cv_tmp[[fold]] = colSums( sweep(x.test%*%betahat, 1, y.test, "-")^2 )
    }
    cv = colMeans(do.call(rbind, cv_tmp))
    # Fit on the full sample
    betahat = coef.restrict(x.train, y.train, R, r)
  }
  return(list(betahat=betahat, i.min=which.min(cv)))
}

## Information criteria (multiple) for all the candidates from the nested models
calc.ic.all <- function(fit, y, df, sigma.sq = NULL) {
  y = matrix(y, ncol = 1)
  df = matrix(df, nrow = 1)
  if (ncol(fit) != ncol(df)) {
    stop("the number of coef vectors does not match the number of df")
  }
  if (is.null(sigma.sq)) {
    stop("need to specify sigma for Mallow's Cp")
  }
  n = nrow(y)
  rss = colSums(sweep(fit, 1, y, "-")^2)
  
  ic = list()
  ic$aic = n*log(rss/n) + 2*(df + 1) + n
  ic$bic = n*log(rss/n) + log(n)*(df + 1) + n
  
  ic$cp = rss/n + 2*sigma.sq*df/n
  ic$rcp = rss/n + sigma.sq*df*(2 + (df + 1)/(n - df - 1))/n
  ic$cphat = rss/n + 2*rss*df/((n-df)*n)
  ic$sp = (n - 1)*rss/((n - df) * (n - df - 1))
  ic$gcv = n*rss/(n - df)^2
  
  # df[which(df >= n - 2)] = n - 3
  ic$aicc = n*log(rss/n) + n*(n + df)/(n - df - 2)
  ic$raicc = n*log(rss/n) + n^2*(n - 1)/((n - df - 1)*(n - df - 2))
  
  return(ic)
}

## Metrics to evaluate the fits
eval.metrics <- function(ii, muhat, betahat, x, y, sigma, mu, beta, Sigma, m){
  n = dim(x)[1]
  p = dim(x)[2]
  
  muhat = muhat[, ii, drop=FALSE]
  betahat = betahat[, ii, drop=FALSE]
  
  output = list()
  output$i_allmethods = ii
  ## Number of variables, Sparsistency and Number of extra variables
  output$nres = m[ii]
  output$nvar = colSums(betahat != 0)
  output$sparsistency = colSums(betahat[which(beta!=0), , drop=FALSE] != 0)
  output$extravariable = colSums(betahat[which(beta==0), , drop=FALSE] != 0)
  
  ## loss
  output$lossF = colSums(sweep(muhat, 1, mu,"-")^2)/n
  lossF_null = sum(mu^2)/n
  output$lossR = unlist( apply(betahat, 2, function(b){ as.numeric( t(b-beta) %*% Sigma %*% (b-beta) )  }) )
  lossR_null = as.numeric( t(beta) %*% Sigma %*% beta )
  output$rlossF = output$lossF/lossF_null
  output$rlossR = output$lossR/lossR_null
  ## KL
  rss = colSums(sweep(muhat, 1, y,"-")^2)
  rss_null = sum(y^2)
  sigmasqhat = rss/n
  sigmasqhat_null = rss_null/n
  Sigmahat = t(x)%*%x/n
  Sigmahat_inv = solve( Sigmahat)
  const = - n*log(sigma^2) - n
  output$KLF = n*log(sigmasqhat) + n*output$lossF/sigmasqhat + n*sigma^2/sigmasqhat + const
  KLF_null = n*log(sigmasqhat_null) + n*lossF_null/sigmasqhat_null + n*sigma^2/sigmasqhat_null + const
  const = n*as.numeric(determinant(Sigmahat, logarithm = TRUE)$modulus)  + n*sum(diag(Sigmahat_inv%*%Sigma)) + const - n*as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) - n*p
  output$KLR = n*log(sigmasqhat) + n*output$lossR/sigmasqhat + n*sigma^2/sigmasqhat + const
  KLR_null = n*log(sigmasqhat_null) + n*lossR_null/sigmasqhat_null + n*sigma^2/sigmasqhat_null + const
  output$rKLF = output$KLF/KLF_null
  output$rKLR = output$KLR/KLR_null
  
  return(output)
}

## Evaluate the selected subsets for all methods
run.all <- function(n, p, rho, snr, type, nrep, random.x=TRUE, restriction.general=FALSE, write.to.file=FALSE){
  
  ## Specify filenames, filepaths and etc
  if(random.x){
    seed_x = (nrep+1):(2*nrep)
    filepath = "randomx/"
  }else{
    seed_x = rep(nrep+1, nrep)
    filepath = "fixedx/"
  }
  snr_all = c(0.2, 1, 8.5)
  names(snr_all) = c("lsnr", "msnr", "hsnr")
  snr_name = names(snr_all)[which(snr_all == snr)]
  
  if(is.null(rho)){
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name)
  }else{
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
  }
  
  if(!restriction.general){
    filepath = paste0(filepath, "subset_selection")
    ## Generate beta and Sigma
    beta_Sigma = gen.beta.Sigma(p, rho, type)
    m = p:0
  }else{
    filepath = paste0(filepath, "general_restriction")
    ## Generate beta and Sigma
    beta_Sigma = gen.beta.Sigma.restriction(p, rho, type)
    m = unlist( lapply(beta_Sigma$restriction_candidate$R, nrow) )
  }
  
  ## Check if some intermediate results already exist
  result = betahat_all = list()
  # if(paste0(filename, ".rds") %in% list.files(paste0(base, "/tmp/", filepath, "/result_intermediate"))){
  #   result = readRDS(paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
  # }
  rep_start = length(result) + 1
  
  for(rep in rep_start:nrep){
    ## Generate data
    data = gen.data(n, p, snr, type, beta_Sigma$beta, beta_Sigma$Sigma.half, seed.y=rep, seed.x=seed_x[rep])
    
    i_allmethods = list()
    ## Cross-validation
    set.seed(2*nrep+rep)
    if(!restriction.general){
      # 10-fold CV
      restrict_cv = cv.nest(data$x, data$y, intercept = FALSE)
      i_allmethods$cv_tenfold = restrict_cv$i.min
      # LOO CV
      restrict_cv = cv.nest(data$x, data$y, n.folds=n, intercept = FALSE)
    }else{
      # 10-fold CV
      restrict_cv = cv.restrict(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r)
      i_allmethods$cv_tenfold = restrict_cv$i.min
      # LOO CV
      restrict_cv = cv.restrict(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r, n.folds=n)
    }
    i_allmethods$cv_loo = restrict_cv$i.min
    
    ## Coefficient vectors for all candidates and corresponding fitted values
    betahat = restrict_cv$betahat 
    muhat = data$x %*% betahat

    ## Information criteria
    sigma_sq_hat = sigma( lm(as.numeric(data$y)~as.matrix(data$x)-1) )^2
    ic_all = calc.ic.all(muhat, data$y, df = p-m , sigma.sq = sigma_sq_hat) 
    i_ic_all = lapply(ic_all, which.min)
    names(i_ic_all) = paste("ic", names(i_ic_all), sep="_")
    i_allmethods = c(i_allmethods, i_ic_all)
    
    ## Evaluate the selection rules
    result[[rep]] = eval.metrics(unlist(i_allmethods), muhat, betahat, data$x, data$y, data$sigma, data$mu, beta_Sigma$beta, beta_Sigma$Sigma, m)
    betahat_all[[rep]] = betahat
    
    if(rep %% 100 == 0){
      print(rep)
      ## Save temporary results
      if(write.to.file){
        saveRDS(result, paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
        saveRDS(betahat_all[(rep-100+1):rep], paste0(base, "/tmp/", filepath, "/betahat/", filename, "_rep", rep, ".rds"))
        betahat_all = list()
      }
    }
  }
  
  ## Adjust the layout of results for the final output
  allmethods = strsplit(names(i_allmethods), split="_")
  result_final = list()
  for(i in 1:length(allmethods)){
    tmp = do.call(cbind, lapply(1:nrep, function(rep){unlist(lapply(result[[rep]], "[[", i))  }))
    if(length(allmethods[[i]]) == 1){
      result_final[[ allmethods[[i]] ]] = split(tmp, row(tmp, as.factor=TRUE))
    }else{
      result_final[[allmethods[[i]][1]]][[allmethods[[i]][2]]] = split(tmp, row(tmp, as.factor=TRUE))
    }
  }
  if(write.to.file){
    saveRDS(result_final, paste0(base, "/results/", filepath, "/", filename, ".rds"))
  }
  
  return(result_final)
}

run.all.existingbetahat <- function(n, p, rho, snr, type, nrep, random.x=TRUE, restriction.general=FALSE, write.to.file=FALSE){

  ## Specify filenames, filepaths and etc
  if(random.x){
    seed_x = (nrep+1):(2*nrep)
    filepath = "randomx/"
  }else{
    seed_x = rep(nrep+1, nrep)
    filepath = "fixedx/"
  }
  snr_all = c(0.2, 1, 8.5)
  names(snr_all) = c("lsnr", "msnr", "hsnr")
  snr_name = names(snr_all)[which(snr_all == snr)]

  if(is.null(rho)){
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name)
  }else{
    filename = paste0(type, "_n", n, "_p", p, "_", snr_name, "_rho", gsub("[.]","",as.character(rho)))
  }

  if(!restriction.general){
    filepath = paste0(filepath, "subset_selection")
    ## Generate beta and Sigma
    beta_Sigma = gen.beta.Sigma(p, rho, type)
    m = p:0
  }else{
    filepath = paste0(filepath, "general_restriction")
    ## Generate beta and Sigma
    beta_Sigma = gen.beta.Sigma.restriction(p, rho, type)
    m = unlist( lapply(beta_Sigma$restriction_candidate$R, nrow) )
  }

  ## Check if some intermediate results already exist
  result = readRDS(paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
  # if(paste0(filename, ".rds") %in% list.files(paste0(base, "/tmp/", filepath, "/result_intermediate"))){
  #   result = readRDS(paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
  # }
  rep_start = 1

  for(rep in rep_start:nrep){
    if(rep %% 100 == 1){
      print(rep)
      betahat_all = list()
      betahat_all[rep:(rep+99)] = readRDS(paste0(base, "/tmp/", filepath, "/betahat/", filename, "_rep", rep+99, ".rds"))
    }

    ## Generate data
    data = gen.data(n, p, snr, type, beta_Sigma$beta, beta_Sigma$Sigma.half, seed.y=rep, seed.x=seed_x[rep])

    i_allmethods = result[[rep]]$i_allmethods

    ## Coefficient vectors for all candidates and corresponding fitted values
    betahat = betahat_all[[rep]]
    muhat = data$x %*% betahat

    ## Evaluate the selection rules
    result[[rep]] = eval.metrics(unlist(i_allmethods), muhat, betahat, data$x, data$y, data$sigma, data$mu, beta_Sigma$beta, beta_Sigma$Sigma, m)
    # betahat_all[[rep]] = betahat

  }

  ## Adjust the layout of results for the final output
  allmethods = strsplit(names(i_allmethods), split="_")
  result_final = list()
  for(i in 1:length(allmethods)){
    tmp = do.call(cbind, lapply(1:nrep, function(rep){unlist(lapply(result[[rep]], "[[", i))  }))
    # if(length(allmethods[[i]]) == 1){
    #   result_final[[ allmethods[[i]] ]] = split(tmp, row(tmp, as.factor=TRUE))
    # }else{
    result_final[[allmethods[[i]][1]]][[allmethods[[i]][2]]] = split(tmp, row(tmp, as.factor=TRUE))
    # }
  }
  if(write.to.file){
    saveRDS(result_final, paste0(base, "/results/", filepath, "/", filename, ".rds"))
  }

  return(result_final)
}
