#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the MSE of transformed EBLUP based on a univariate Fay-Herriot model. For sampled domains, MSE is estimated using modified parametric bootstrap approach proposed by Butar & Lahiri. For non-sampled domains, MSE is estimated using modified approach proposed by Haris & Ubaidillah.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param sigma2u.awal vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param df vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param maxiter vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param tolerance vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param data vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @export mse.tsae
#' @import utils
#' @import stats



##MSE EBLUP t-sae
mse.tsae<-function (formula,sigma2u.awal,vardir,df,maxiter,tolerance, data)
{
  result <- list(est = NA, mse = NA)
  ##namevar <- deparse(substitute(vardir))

  if(!missing(data)){
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    ##vardir <- data[,namevar]
  } else{
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula, data)
  }
  y <- formuladata[, 1]
  result$est <- ecm.tsae(y ~ 0 + X, vardir, sigma2u.awal, df, maxiter,tolerance,data)
  sigma2u <- result$est$fit$refvar
  beta <- result$est$fit$estcoef$beta
  m <- dim(X)[1]
  p <- dim(X)[2]
  v<-rep(df,m)
  satu<-rep(1,m)
  var.u<-rep(sigma2u,m)
  g1i <- rep(0, m)
  g2i <- rep(0, m)
  g3i <- rep(0, m)
  g4i <- rep(0, m)
  mse2i <- rep(0, m)
  Vi<-1/(var.u+vardir)
  Bi<- vardir/(var.u + vardir)
  XtVi<-t(Vi*X)
  Q <-XtVi %*% X
  Q1 <- solve(Q)
  w1 <-(df+3)/(df+1)

  Xbeta <- X %*% beta
  resid.EM <- y - Xbeta
  z2<-resid.EM^2*Vi
  omega<-1/((1/var.u)+(1/vardir))
  tau<-(v+satu)/(v+z2)
  Wi <- as.vector(tau * Vi)
  XtWi <- t(Wi * X)
  Q3 <- solve(XtWi %*% X)
  Wi2 <- as.vector(tau * Vi^2)
  Vinv<-diag(Vi)
  w2<- df/(df-2)
  w3<- 1/(df*(df+3))
  w4<- (df+2)/(df*(df+3))
  w5<- (df+1)*(df+2)/(df*(df+3))
  Winv<-diag(Wi)
  W<-solve(Winv)
  Vinv<-diag(Vi)
  Vinv2<-diag(Vi^2)
  P <- Winv - Winv %*% X %*% Q3 %*% XtWi
  D <- sum(diag(Vinv %*% Vinv))
  E <- sum(diag( W %*% Vinv %*% Vinv %*% Vinv %*% W %*% P))
  F <- sum(diag( W %*% Vinv %*% Vinv %*% W %*% P))
  J <- sum(diag( W %*% Vinv %*% Vinv %*% W %*% P %*% W %*% Vinv %*%
                   Vinv %*% W %*% P))
  A <- (- 0.5)* D + E - (0.5)* w4 *(2*J + F^2)
  c <- (- 0.5)*sum(diag(Vinv)) + (0.5)* F
  cc <- (c^2) + (0.5)*(w3*F^2) +(.5)* w5 * J
  b <- c/A #asymptotic bias of sigm2u
  VarSigma2u<-cc/(A^2) #asymptotic variance of sigma2u
  for (i in 1:m) {
    g1i[i] <- w2 * omega[i]
    xi <- matrix(X[i, ], nrow = 1, ncol = p)
    g2i[i] <- w1*(Bi[i]^2) * xi %*% Q1 %*% t(xi)
    g3i[i]<-(Bi[i]^2)*VarSigma2u*sum(diag(W%*%Vinv%*%Vinv%*%W%*%P))
    g4i[i] <- w4*(Bi[i]^2)*VarSigma2u /(sigma2u + vardir[i])
    mse2i[i] <- g1i[i] + g2i[i]+ g3i[i] + g4i[i] - b*(Bi[i]^2)*w4
  }
  result$mse <- mse2i
  return(result)
}
