library(MASS)
library(Matrix)

## Generate simulated dataset
gen.beta.Sigma <- function(p, rho, type){
  
  # Covariance matrix and true beta for each true model configuration
  Sigma = Matrix(rho^abs(outer(1:p, 1:p, "-")), sparse = TRUE)
  if(type == "VS-Ex1" ){
    beta = Matrix(c(rev(seq(1, 3.5, by=0.5)), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "VS-Ex2"){
    Sigma = bdiag(list(rho^abs(outer(1:6, 1:6, "-")), rho^abs(outer(7:p, 7:p, "-"))))
    beta = Matrix(c(c(1,1,3,3,5,5), rep(0,p-6)), ncol=1, sparse = TRUE)
  }else if(type == "VS-Ex3"){
    kappa = 10
    beta = Matrix((-1)^seq(1,p) * exp(-seq(1,p)/kappa), ncol=1, sparse = TRUE)
  }else if(type == "GR-Ex1"){
    beta = Matrix(c(2,2,2,1,1,1), ncol=1)
    R_candidate = Matrix(rbind(c(1,-1,0,0,0,0), c(0,1,-1,0,0,0), c(0,0,0,1,-1,0), 
                               c(0,0,0,0,1,-1), c(1,0,0,-1,0,0), c(1,-2,0,0,0,0)), sparse=TRUE)
    r_candidate = Matrix(rep(0, 6), ncol=1)
  }else if(type == "GR-Ex2"){
    beta = Matrix(c(-2,2,2,-1.5,-1.5,1), ncol=1)
    R_candidate = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,1,1), c(0,0,0,0,1,1), 
                               c(0,0,0,1,0,1), c(-2,0,0,1,0,0), c(1,-1,0,0,0,0)), sparse=TRUE)
    r_candidate = Matrix(rep(0, 6), ncol=1)
  }else if(type == "GR-Ex3"){
    beta = Matrix(c(2,-2,1,-1,0.5,-0.5), ncol=1)
    R_candidate = Matrix(rbind(c(1,1,0,0,0,0), c(0,0,1,1,0,0), c(0,0,0,0,1,1), 
                               c(1,0,0,1,0,0), c(0,1,0,0,1,0), c(1,-1,0,0,0,0)), sparse=TRUE)
    r_candidate = Matrix(rep(0, 6), ncol=1)
  }else if(type == "GR-Ex4"){
    beta = Matrix(c(2,2,2,1,1,1, rep(0,p-6)), ncol=1)
    R_candidate = cbind(Matrix(0, nrow=p-6, ncol=6), apply(Diagonal(p-6), 1, rev))
    R_candidate = rbind(R_candidate, 
                        cbind(Matrix(rbind(c(0,0,0,0,1,-1), c(0,0,0,1,-1,0), c(0,1,-1,0,0,0), 
                                           c(1,-1,0,0,0,0), c(1,-2,0,0,0,0), c(1,0,0,-1,0,0))), 
                              Matrix(0, nrow=6, ncol=p-6))
    )
    r_candidate = Matrix(rep(0, p), ncol=1)
  }else if(type == "GR-Ex5"){
    beta = Matrix(c(-2,2,2,-1.5,-1.5,1, rep(0,p-6)), ncol=1)
    R_candidate = cbind(Matrix(0, nrow=p-6, ncol=6), apply(Diagonal(p-6), 1, rev))
    R_candidate = rbind(R_candidate, 
                        cbind(Matrix(rbind(c(0,0,1,1,1,1), c(1,1,0,0,0,0), c(0,0,0,0,1,1), 
                                           c(0,0,0,1,0,1), c(-2,0,0,1,0,0), c(1,-1,0,0,0,0))), 
                              Matrix(0, nrow=6, ncol=p-6))
    )
    r_candidate = Matrix(rep(0, p), ncol=1)
  }else if(type == "GR-Ex6"){
    beta = Matrix(c(2,-2,1,-1,0.5,-0.5, rep(0,p-6)), ncol=1)
    R_candidate = cbind(Matrix(0, nrow=p-6, ncol=6), apply(Diagonal(p-6), 1, rev))
    R_candidate = rbind(R_candidate, 
                        cbind(Matrix(rbind(c(0,0,0,0,1,1), c(0,0,1,1,0,0), c(1,1,0,0,0,0), 
                                           c(0,1,0,0,1,0), c(1,0,0,1,0,0), c(1,-1,0,0,0,0))), 
                              Matrix(0, nrow=6, ncol=p-6))
    )
    r_candidate = Matrix(rep(0, p), ncol=1)
  }
  
  # SVD of the covariance matrix Sigma
  if(rho != 0){
    obj = svd(Sigma)
    Sigma.half = Matrix(obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v), sparse=TRUE)
  }else{
    Sigma.half = NULL
  }
  
  # Generate a list of restrictions (R, r) from the candidates
  R_all = r_all = list()
  if(type %in% paste0("GR-Ex", seq(1,3))){
    for(i in 0:nrow(R_candidate)){
      if(i == 0){
        R_all = c(R_all, list(Matrix(, nrow=0, ncol=p)))
        r_all = c(r_all, list(numeric(0)))
      }else{
        tmp = combn(nrow(R_candidate),i)
        R_all = c(R_all, lapply(split(tmp, col(tmp)), function(ii){R_candidate[ii, , drop=FALSE]}))
        r_all = c(r_all, lapply(split(tmp, col(tmp)), function(ii){r_candidate[ii, , drop=FALSE]}))
      }
    }
  }else if(type %in% paste0("GR-Ex", seq(4,6))){
    R_all = list(Matrix(, nrow=0, ncol=p))
    r_all = list(numeric(0))
    
    R_all = c(R_all, lapply(1:p, function(ii){R_candidate[1:ii, , drop=FALSE]}))
    r_all = c(r_all, lapply(1:p, function(ii){r_candidate[1:ii, , drop=FALSE]}))
    
    R_all = rev(R_all)
    r_all = rev(r_all)
  }
  
  return(list(beta=beta, Sigma=Sigma, Sigma.half=Sigma.half, restriction_candidate=list(R=R_all, r=r_all)))
}
gen.data <- function(n, p, snr, type, beta, Sigma.half, seed.x, seed.y){
  
  # generate design matrix X
  set.seed(seed.x)
  x = matrix( rnorm(n*p), nrow=n, ncol=p )
  if(!is.null(Sigma.half)){
    x = x %*% Sigma.half
  }
  mu = x %*% beta
  sigma = as.numeric(sqrt(var(as.vector(mu)) / snr))
  
  # generate the response
  set.seed(seed.y)
  y = mu + Matrix(rnorm(n, mean=0, sd=sigma), ncol=1)
  
  return(list(x=x, y=y, mu=mu, sigma=sigma))
}

## Fit LS on the nested subsets of X
coef.subset.nest <- function(x, y, intercept=TRUE){
  
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
## (R, r) here are in the format of list
coef.restrict <- function(x, y, R, r){          
  n = dim(x)[1]
  p = dim(x)[2]
  maxstep = min(n, p)

  # No-restriction case
  ind_null = which( unlist(lapply(R, function(RR){dim(RR)[1]}))==0 )
  # Calculate the coefficient vector
  xtx_inv = ginv(as.matrix(t(x) %*% x))
  # OLS estimator on the full set of data
  betahat_f = xtx_inv %*% t(x) %*% y
  tmp_function <- function(RR, rr){
    betahat_f + xtx_inv %*% t(RR) %*% ginv( as.matrix( RR %*% xtx_inv %*% t(RR) ) ) %*% (rr - RR %*% betahat_f)
  }
  betahat = do.call(cbind, Map(tmp_function, R[-ind_null], r[-ind_null]))

  # Unrestricted OLS for no-restriction cases
  betahat_final = Matrix(0, nrow=maxstep, ncol=length(R))
  betahat_final[,ind_null] = matrix(rep(betahat_f, length(ind_null)), nrow=maxstep, byrow=FALSE)
  betahat_final[,-ind_null] = betahat
  betahat_final = Matrix(betahat_final, sparse = TRUE)
  
  return(betahat_final)
}

## Combinations of general restriction and variable selection
## This function is by no means useful in practice.
## The intention is to only run fast.
## It is built on the specific designs of the simulations
## e.g., the general restrictions are imposed on first six predictors
## while subset restrictions (setting to be zero) are imposed to the rest of predictors
## also the restrictions are nested for the candidate models
coef.restrict.nest <- function(x, y, R, r){          
  n = dim(x)[1]
  p = dim(x)[2]
  
  # Take the general restrictions (remove all the restrictions w.r.t. variable selection)
  R_sub = R[[1]][(nrow(R[[1]])-6+1):p, 1:6]
  
  # Singular value decomposition
  xtx_sub = t(x[,1:6]) %*% x[,1:6]
  svd = eigen(xtx_sub)
  P = svd$vectors
  Linv = Diagonal(x=1/svd$values)
  xty_sub = t(x[,1:6]) %*% y
  
  # Compotents in the formula of betahat
  PLinvPt = P %*% Linv %*% t(P)
  RPLinvPt = R_sub %*% PLinvPt
  RPLinvPtRt = RPLinvPt %*% t(R_sub)
  betahat_f = PLinvPt %*% xty_sub
  
  # Iteratively calculate the inverse of R(X^T X)^{-1}R^T, using the property of the inverse of block matrices
  # and calculate betahat correspondingly
  Ainv = 1 / RPLinvPtRt[1,1]
  xtxinvRt_Ainv_Rxtxinv = t(RPLinvPt[1,,drop=FALSE]) %*% Ainv %*% RPLinvPt[1,,drop=FALSE]
  betahat = Matrix(0, nrow=6, ncol=6)
  betahat[,6] = betahat_f - xtxinvRt_Ainv_Rxtxinv %*% xty_sub
  for(i in 2:6){
    Ainvb = Ainv %*% RPLinvPtRt[1:(i-1),i]
    k = as.numeric( RPLinvPtRt[i,i] - t(RPLinvPtRt[1:(i-1),i]) %*% Ainvb )
    Ainv = rbind(cbind(Ainv  + Ainvb %*% t(Ainvb) / k, -Ainvb / k), cbind(-t(Ainvb) / k, 1 / k))
    xtxinvRt_Ainv_Rxtxinv = t(RPLinvPt[1:i,]) %*% Ainv %*% RPLinvPt[1:i,]
    betahat[,6-i+1] = betahat_f - xtxinvRt_Ainv_Rxtxinv %*% xty_sub
  }
  
  # Work out the coefficients for the remaining restrictions (variable selection)
  betahat_final = coef.subset.nest(x, y, intercept = FALSE)
  betahat_final[1:6,1:6] = betahat
 
  return(betahat_final)
}

## Cross-validation to choose the best candidate from the nested models
cv.subset.nest <- function(x, y, n.folds=10, intercept=FALSE){
  n = dim(x)[1]
  p = dim(x)[2]

  if(n.folds == n){
    ## Use PRESS statistic for LOO CV
    betahat = coef.subset.nest(x, y, intercept=intercept)
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
      betahat = coef.subset.nest(x.train, y.train, intercept=intercept)
      
      # Evaluate on the testing set
      if(!intercept){
        cv_tmp[[fold]] = colSums( sweep(x.test%*%betahat, 1, y.test, "-")^2 )
      }else{
        cv_tmp[[fold]] = colSums( sweep(cbind(rep(1,nrow(x.test)),x.test)%*%betahat, 1, y.test, "-")^2 )
      }
    }
    cv = colMeans(do.call(rbind, cv_tmp))
    # Fit on the full sample
    betahat = coef.subset.nest(x, y, intercept=intercept)
  }
  return(list(betahat=betahat, i.min=which.min(cv), cvm=cv))
}

## Cross-validation to choose the best model, each of which is under a set of restrictions
cv.restrict <- function(x, y, R, r, n.folds=10){
  n = dim(x)[1]
  p = dim(x)[2]
  
  if(n.folds == n){
    # Use PRESS statistic for LOO-CV
    # See Tarpey (2000) for details
    xtx_inv = solve(t(x) %*% x)
    H_base = xtx_inv %*% t(x)
    betahat_f = H_base %*% y
    H = x %*% H_base
    ind_null = which( unlist(lapply(R, function(RR){dim(RR)[1]}))==0 )
    
    # Expressions in Tarpey (2000)
    calc.H <- function(RR){
      tmp1 = xtx_inv %*% t(RR)
      tmp2 = x %*% tmp1
      tmp3 = solve(RR %*% tmp1)
      
      HRQ_base = tmp1 %*% tmp3
      H_R = HRQ_base %*% RR
      H_Q = tmp2 %*% tmp3 %*% t(tmp2)
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

## Cross-validation to choose the best model, each of which is under a set of restrict.
## Again, similar to 'coef.restrict.nest', this function is intended to run fast, 
## and does not generally have practical value
## It is built on knowing the first 6 restrictions are general, and the rest are variable selection
cv.restrict.nest <- function(x, y, R, r, n.folds=10){
  n = dim(x)[1]
  p = dim(x)[2]
  
  if(n.folds == n){
    # LOO CV
    
    R_sub = R[[1]][(nrow(R[[1]])-6+1):p, 1:6]
    # Singular value decomposition
    xtx_sub = t(x[,1:6]) %*% x[,1:6]
    svd = eigen(xtx_sub)
    P = svd$vectors
    Linv = Diagonal(x=1/svd$values)
    xty_sub = t(x[,1:6]) %*% y
    
    # compotents in the formula of betahat
    PLinvPt = P %*% Linv %*% t(P)
    RPLinvPt = R_sub %*% PLinvPt
    RPLinvPtRt = RPLinvPt %*% t(R_sub)
    H_sub = x[,1:6] %*% PLinvPt %*% t(x[,1:6])
    betahat_f = PLinvPt %*% xty_sub
    
    # iteratively calculate the inverse of R(X^T X)^{-1}R^T, using the property of the inverse of block matrices
    # and betahat
    Ainv = 1 / RPLinvPtRt[1,1]
    xtxinvRt_Ainv_Rxtxinv = t(RPLinvPt[1,,drop=FALSE]) %*% Ainv %*% RPLinvPt[1,,drop=FALSE]
    H_Q = x[,1:6] %*% xtxinvRt_Ainv_Rxtxinv %*% t(x[,1:6])
    betahat = Matrix(0, nrow=6, ncol=6)
    cv = rep(0, p+1)
    betahat[,6] = betahat_f - xtxinvRt_Ainv_Rxtxinv %*% xty_sub
    cv[6] = mean(((y - x[,1:6] %*% betahat[,6]) / (1 - diag(H_sub) + diag(H_Q)))^2)
    for(i in 2:6){
      Ainvb = Ainv %*% RPLinvPtRt[1:(i-1),i]
      k = as.numeric( RPLinvPtRt[i,i] - t(RPLinvPtRt[1:(i-1),i]) %*% Ainvb )
      Ainv = rbind(cbind(Ainv  + Ainvb %*% t(Ainvb) / k, -Ainvb / k), cbind(-t(Ainvb) / k, 1 / k))
      xtxinvRt_Ainv_Rxtxinv = t(RPLinvPt[1:i,]) %*% Ainv %*% RPLinvPt[1:i,]
      H_Q = x[,1:6] %*% xtxinvRt_Ainv_Rxtxinv %*% t(x[,1:6])
      betahat[,6-i+1] = betahat_f - xtxinvRt_Ainv_Rxtxinv %*% xty_sub
      cv[6-i+1] = mean( ((y - x[,1:6] %*% betahat[,6-i+1]) / (1 - diag(H_sub) + diag(H_Q)))^2 )
    }
    
    # calculate the CV error for the remaining restrictions
    tmp = cv.subset.nest(x, y, n.folds=n)
    tmp$betahat[1:6,1:6] = betahat
    betahat = tmp$betahat
    cv[7:(p+1)] = tmp$cvm[7:(p+1)]
    # print(cv)
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
      betahat = coef.restrict.nest(x.train, y.train, R, r)
      # Evaluate on the testing set
      cv_tmp[[fold]] = colSums( sweep(x.test%*%betahat, 1, y.test, "-")^2 )
    }
    cv = colMeans(do.call(rbind, cv_tmp))
    # Fit on the full sample
    betahat = coef.restrict.nest(x.train, y.train, R, r)
  }
  return(list(betahat=betahat, i.min=which.min(cv)))
}

## Calculate a list of information criteria for the candidate models
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
  
  df[which(df >= n - 2)] = n - 3
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

## This is the main function to fit the models and evaluate the performance
run.all <- function(n, p, rho, snr, type, nrep, random.x=TRUE, write.to.file=FALSE){
  
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
  
  restriction.general = grepl("GR-", type)
  beta_Sigma = gen.beta.Sigma(p, rho, type)
  if(!restriction.general){
    m = p:0
  }else{
    m = unlist( lapply(beta_Sigma$restriction_candidate$R, nrow) )
  }
  is_nest = ifelse(type %in% paste0("GR-Ex", seq(4,6)), TRUE, FALSE)

  ## Check if some intermediate results already exist
  result = betahat_all = list()
  if(paste0(filename, ".rds") %in% list.files(paste0(base, "/tmp/", filepath, "/result_intermediate"))){
    result = readRDS(paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
    i_allmethods = result[[1]]$i_allmethods
  }
  rep = length(result) + 1
  
  while(rep <= nrep){
    ## Generate data
    data = gen.data(n, p, snr, type, beta_Sigma$beta, beta_Sigma$Sigma.half, seed.y=rep, seed.x=seed_x[rep])
    
    i_allmethods = list()
    ## Cross-validation
    set.seed(2*nrep+rep)
    if(!restriction.general){
      # 10-fold CV
      restrict_cv = cv.subset.nest(data$x, data$y, intercept = FALSE)
      i_allmethods$cv_tenfold = restrict_cv$i.min
      # LOO CV
      restrict_cv = cv.subset.nest(data$x, data$y, n.folds=n, intercept = FALSE)
    }else{
      if(is_nest){
        # 10-fold CV
        restrict_cv = cv.restrict.nest(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r)
        i_allmethods$cv_tenfold = restrict_cv$i.min
        # LOO CV
        restrict_cv = cv.restrict.nest(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r, n.folds=n)
      }else{
        # 10-fold CV
        restrict_cv = cv.restrict(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r)
        i_allmethods$cv_tenfold = restrict_cv$i.min
        # LOO CV
        restrict_cv = cv.restrict(data$x, data$y, beta_Sigma$restriction_candidate$R, beta_Sigma$restriction_candidate$r, n.folds=n)
      }

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
      ## Save intermediate results
      if(write.to.file){
        saveRDS(result, paste0(base, "/tmp/", filepath, "/result_intermediate/", filename, ".rds"))
        saveRDS(betahat_all[(rep-100+1):rep], paste0(base, "/tmp/", filepath, "/betahat/", filename, "_rep", rep, ".rds"))
        betahat_all = list()
      }
    }
    
    rep = rep + 1
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
