#' @title Parametric Bootstrap Mean Squared Error of EBLUPs based on a Univariate Fay Herriot model with Additive Logistic Transformation for Non-Sampled Data
#' @description This function gives the MSE of transformed EBLUP based on a univariate Fay-Herriot model. For sampled domains, MSE is estimated using modified parametric bootstrap approach proposed by Butar & Lahiri. For non-sampled domains, MSE is estimated using modified approach proposed by Haris & Ubaidillah.
#' @param formula an object of class \code{\link[stats]{formula}} that describe the fitted model.
#' @param vardir vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param sigma2u vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param df vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param maxiter vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param tolerance vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
#' @param data vector containing the sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.
ecm.tsae<- function(formula,vardir,sigma2u,df=3,maxiter=1000,tolerance=10^(-4), data){


  result <- list(eblup.t = NA, fit = list(iterations=0, estcoef = NA, refvar = NA, goodness = NA))
  if(!missing(data)){
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
  }

  else{
    formuladata <- model.frame(formula, na.action = na.omit)
    X <- model.matrix(formula, data)
  }

  y <- formuladata[, 1]

  if (attr(attributes(formuladata)$terms, "response") == 1){
    textformula <- paste(formula[2], formula[1], formula[3])
  } else{ textformula <- paste(formula[1], formula[2])
  }

  if (length(na.action(formuladata)) > 0){
    stop("Argument formula=", textformula, " contains NA values.")
  }

  if (any(is.na(vardir))){
    stop("Argument vardir=", namevar, " contains NA values.")
  }

  m<-length(y)
  p <- dim(X)[2]
  v<-rep(df,m)
  satu<-rep(1,m)
  conv <- 1

  ##BETA akan dihitung disini
  rlb<-lm(formula)
  beta<-matrix(ncol=1,nrow=length(coef(rlb)))
  for (i in 1:nrow(beta)){
    beta[i]=coef(rlb)[i]
  }

  L0 <- loglike(y,X, beta, vardir, sigma2u, df)
  i<-0
  repeat {
    if(i>maxiter) {conv<-0
    break
    }else {
      Xb<- X %*% (beta)
      resid<- y - Xb
      var.u <- rep(sigma2u,m)
      Vi<-1/(var.u + vardir)
      Bi<- vardir/(var.u + vardir)
      delta<-resid^2*Vi
      #E-step
      u<-resid*(1-Bi)
      omega<-1/((1/var.u)+(1/vardir))
      tau<-(v+satu)/(v+delta)
      #CM-step
      w<-as.vector(tau/vardir)
      Xtw<-t(w*X)
      Q<-solve(Xtw %*% X)
      beta<-Q %*% Xtw %*% (y-u)
      sigma2.u<-sum(tau*(u^2) + omega)/m
      sigma2u<-max(sigma2.u,0)
      L1 <- loglike(y, X, beta, vardir, sigma2u, df)
      if(L1 < L0) { print("log-likelihood must increase, llikel < llikeO, break.")
        conv <- 0
        break
      }
    }
    i <- i + 1
    result$fit$iterations <- i
    if(abs(L1 - L0) < tolerance) {break} #check for convergence
    L0 <- L1
  }
  var.u<-rep(sigma2u,m)
  Vi<-1/(var.u+vardir)
  XtVi<-t(Vi*X)
  w3<-(df+1)/(df+3)
  Q1<-XtVi %*% X
  w3Q1<-w3*Q1
  Q2<-solve(w3Q1)
  Iinv<-Q2
  std.error.beta<-sqrt(diag(Iinv))
  tvalue <- beta/std.error.beta
  pvalue <- 2 * pnorm(abs(tvalue), lower.tail = FALSE)
  Xbeta <- X %*% beta
  resid.EM <- y - Xbeta
  AIC <- (-2) * L0 + 2 * (p + 1)
  BIC <- (-2) * L0 + (p + 1) * log(m)
  goodness <- c(Loglikelihood = L0, AIC = AIC, BIC = BIC)
  coef <- data.frame(beta = beta, std.error = std.error.beta,tvalue,
                     pvalue)
  variance <- sigma2u
  EBLUPt <- Xbeta + sigma2u * Vi * resid.EM
  result$fit$estcoef <- coef
  result$fit$refvar <- variance
  result$fit$goodness <- goodness
  result$eblup.t <- EBLUPt
  return(result)
}
#########################################################
# loglike calculates the LogLikelihood #
#########################################################

loglike<- function(y,X, beta, vardir, sigma2u, df)
{
  m<- length(y)
  v<-rep(df,m)
  satu<-rep(1,m)
  Xb<- X %*% beta
  resid<- y - Xb
  var.u <- rep(sigma2u,m)
  Vi<-1/(var.u + vardir)
  Bi<- vardir/(var.u + vardir)
  delta<-resid^2*Vi
  l<-(-.5)*log(1/Vi)-(.5)*(v+satu)*log(1+(delta)/v)
  sum(l)
}
