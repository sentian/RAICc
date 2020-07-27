library(MASS)
library(limSolve)
n = 20
p = 3
# x = matrix(rnorm(n*p), nrow = n)

# mu = x %*% beta
# y = mu + rnorm(n)

R = matrix(rnorm(2*p), nrow=2)
r = matrix(rnorm(2), ncol=1)
beta = Solve(R, r)

## beta has to satisfy R*beta = r
R = rbind(c(0,1,0), c(0,0,1))
r = matrix(c(1,-1), ncol=1)
beta = matrix(c(1, 1, -1), ncol=1)

R = matrix(c(1,-1,0),nrow=1)
r = matrix(0, ncol=1)
beta = matrix(c(1, 1, 0), ncol=1)

R = matrix(c(1,1,1),nrow=1)
r = matrix(3, ncol=1)
beta = matrix(c(1, 1, 1), ncol=1)

# covmatrix = diag(p)
covmatrix = matrix(0.2, nrow=p, ncol=p)
diag(covmatrix) = 1
nrep = 100000
#tmp1 = tmp2 = tmp3 = list()

tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = c()
tmp = c()
for(rep in 1:nrep){
  x = mvrnorm(n, mu=rep(0,p), Sigma=covmatrix)
  mu = x %*% beta
  xtx = t(x) %*% x
  xtx_inv = ginv(xtx)
  y = mu + rnorm(n)
  
  betahat_f = xtx_inv %*% t(x) %*% y
  betahat = betahat_f - xtx_inv %*% t(R) %*% ginv(R %*% xtx_inv %*% t(R)) %*% (R %*% betahat_f - r )
  rss = sum((y - x%*%betahat)^2)
  # tmp1[rep] = as.numeric( t(y-mu) %*% (y-mu) ) # chi-square(n)
  # tmp2[rep] = as.numeric( t(y-x%*%betahat_f) %*% (y-x%*%betahat_f) ) # chi-square(n-p)
  # tmp3[rep] = as.numeric( t(betahat - beta) %*% xtx %*% (betahat - beta) ) # chi-square(p-m)
  # tmp4[rep] = as.numeric( t(betahat - betahat_f) %*% xtx %*% (betahat - betahat_f) ) # chi-square(m)
  # tmp5[rep] = as.numeric( t(y - x %*% betahat) %*% (y - x %*% betahat) ) # chi-square(n-p+m)
  # tmp6[rep] = as.numeric( t(betahat - beta) %*% covmatrix %*% (betahat - beta) / rss ) 
  tmp7[rep] = as.numeric( t(betahat - beta) %*% covmatrix %*% (betahat - beta) ) 
  tmp8[rep] = sum(diag( covmatrix %*% (xtx_inv - xtx_inv %*% t(R) %*% ginv(R %*% xtx_inv %*% t(R)) %*% R %*% xtx_inv) ))
  tmp[rep] = as.numeric( n^2*t(betahat - beta) %*% covmatrix %*% (betahat - beta) / rss + n^2 / rss) 
  
  # tmp[[rep]] = covmatrix %*% (xtx_inv - xtx_inv %*% t(R) %*% ginv(R %*% xtx_inv %*% t(R)) %*% R %*% xtx_inv)
}
mean(tmp)
n^2*(n-1) / ((n - p + nrow(R) - 1)*(n - p + nrow(R) - 2)) 

mean(tmp1)
n
mean(tmp2)
n - p
mean(tmp3)
p - nrow(R)
mean(tmp4)
nrow(R)
mean(tmp5)
n - p + nrow(R)
mean(tmp6)
(p - nrow(R)) / ((n - p + nrow(R) - 1)*(n - p + nrow(R) - 2))
mean(tmp7)
mean(tmp8)
(p - nrow(R)) / (n - p + nrow(R) - 1)

(p - nrow(R)) / (n - p - 1)

proj = cbind(matrix(0, nrow = p-nrow(R), ncol = nrow(R)), diag(p-nrow(R)))
covmatrix = matrix(0, nrow=p, ncol=p)
diag(covmatrix) = 1
nrep = 100000
target = list()
for(rep in 1:nrep){
  x = mvrnorm(n, mu=rep(0,p), Sigma=covmatrix)
  xtx = t(x) %*% x
  tmp = eigen(xtx)
  P = tmp$vectors
  V_inv = P %*% diag(1/sqrt(tmp$values)) %*% t(P)
  H_1 = R %*% V_inv
  tmp = eigen( t(H_1) %*% solve(H_1 %*% t(H_1)) %*% H_1 )
  Q = tmp$vectors
  
  M = t(Q) %*% V_inv %*% covmatrix %*% V_inv %*% Q
  target[[rep]] = M
}
sum(diag(Reduce("+", target) / nrep)[(nrow(R)+1):p])
(p - nrow(R)) / (n - p + nrow(R) - 1)

sum( diag( Reduce("+", lapply(target, function(xx){ proj %*% xx %*% t(proj) })) / nrep ) )
(p - nrow(R)) / (n - p - 1)

diag( Reduce("+", target) / nrep )


x1 = x[,1:2]
x2 = x[,3:10]
h1 = x1 %*% solve(t(x1) %*% x1) %*% t(x1)
h2 = x2 %*% solve(t(x2) %*% x2) %*% t(x2)
solve(t(x2) %*% (diag(n) - h1) %*% x2) %*% t(x2) %*% h1 %*% (diag(n) - h2) %*% x1

xtx = t(x) %*% x
R = rbind(c(1,2,0), c(0,1,0))
Rc = matrix(c(0,0,1), nrow=1)
r = matrix(c(1,1), ncol=1)

Rc %*% xtx %*% t(R) %*% solve(R %*% xtx %*% t(R)) %*% R %*% xtx %*% t(Rc) %*% solve(Rc %*% xtx %*% t(Rc))

xtilde1 = x %*% t(R)
xtilde2 = x %*% t(Rc)
xtilde = cbind(xtilde1, xtilde2)
Rtilde = rbind(R, Rc)
mtilde = R %*% t(Rtilde)
htilde1 = xtilde1 %*% solve(t(xtilde1) %*% xtilde1) %*% t(xtilde1)
htilde2 = xtilde2 %*% solve(t(xtilde2) %*% xtilde2) %*% t(xtilde2)

solve(t(xtilde1) %*% (diag(n) - htilde2) %*% xtilde1) %*% (t(xtilde1) %*% y - t(xtilde1) %*% htilde2 %*% y)
solve(t(xtilde2) %*% (diag(n) - htilde1) %*% xtilde2) %*% (t(xtilde2) %*% y - t(xtilde2) %*% htilde1 %*% y)

solve(t(xtilde2) %*% (diag(n) - htilde1) %*% xtilde2) %*%  t(xtilde2) %*% (htilde1 - htilde1 %*% htilde2) %*% xtilde1 %*% solve(R %*% t(R))

solve(t(xtilde) %*% xtilde) %*% t(mtilde) %*% solve(mtilde %*% solve(t(xtilde) %*% xtilde) %*% t(mtilde)) %*% r



xtildetxtilde_inv = solve(t(xtilde) %*% xtilde)
xtildetxtilde_inv %*% t(xtilde) %*% y
xtildetxtilde_inv %*% t(xtilde) %*% y - xtildetxtilde_inv %*% t(mtilde) %*% solve(mtilde %*% xtildetxtilde_inv %*% t(mtilde)) %*% (mtilde %*% xtildetxtilde_inv %*% t(xtilde) %*% y - r)

lm(y ~ xtilde1 - 1) 
  
solve(t(xtilde2) %*% xtilde2) %*% t(xtilde2) %*% xtilde1 %*% solve(R %*% t(R)) %*% r

solve(t(xtilde2) %*% (diag(n) - htilde1) %*% xtilde2) %*% t(xtilde2) %*% htilde1 %*% (diag(n) - htilde2) %*% xtilde1 %*% solve(R %*% t(R)) %*% r

t(xtilde2) %*% htilde1 %*% (diag(n) - htilde2) %*% xtilde1
t(xtilde2) %*% (diag(n) - htilde1) %*% xtilde2 %*% solve(t(xtilde2) %*% xtilde2) %*% t(xtilde2) %*% xtilde1 

R %*% covmatrix %*% t(R) %*% solve(covmatrix[2:3,2:3])

R
Rc
