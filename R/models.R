pgVBbase <- function(X, Y, sig2b=10, eps=0.01, wgt=NULL){
  n <- length(Y)
  if(is.null(wgt)) wgt <- rep(1, n)
  p <- ncol(X)
  C <- X
  MUbu <- rep(1,p) 
  SIGbu <- Diagonal(p)
  checkOld <- Inf
  checkVec <- c()
  iter <- 1
  repeat{
    ## Latent Variables
    xi <- as.numeric(sqrt(colSums(t(C)*(SIGbu%*%t(C))) + (C%*%MUbu)^2))
    
    ## Regression coefficients
    Zbar <- Diagonal(n,(wgt*0.5/xi)*tanh(0.5*xi))
    SIGbu <- solve((1/sig2b)*Diagonal(p) + t(C)%*%Zbar%*%C)
    MUbu <-  SIGbu %*% t(C) %*% as.numeric(wgt*(Y - 0.5))
    
    
    
    PPrec <- (1/sig2b)*Diagonal(p)
    
    ## Check for convergence
    checkNew  <- 0.5*(p) + 0.5*determinant(SIGbu, logarithm = T)$modulus  - 
      0.5*t(MUbu)%*%PPrec%*%(MUbu) + sum(wgt*(Y-0.5)*as.numeric(C%*%MUbu) +wgt*log(sigmoid(xi)) - 0.5*wgt*xi) - 
      0.5*sum(diag(PPrec %*% SIGbu)) 
    checkVec <- c(checkVec, as.numeric(checkNew))
    
    
    if (as.logical(abs(checkOld - checkNew) < eps)) break
    iter <- iter + 1
    checkOld <- checkNew
  }
  return(list(iter=iter, BU=list(type="Normal", mean=MUbu, Cov=SIGbu, p=p), Elbo=checkVec))
}


pgVB <- function(X, Z, Y, sig2b=10, Au=0.5, Bu=0.5, eps=0.01, wgt=NULL){
  n <- length(Y)
  if(is.null(wgt)) wgt <- rep(1, n)
  p <- ncol(X)
  r <- ncol(Z)
  C <- cbind(X,Z)
  Bsig2u <- 1
  MUbu <- rep(1,p+r) 
  SIGbu <- Diagonal(p+r)
  checkOld <- Inf
  checkVec <- c()
  iter <- 1
  repeat{
    ## Latent Variables
    xi <- as.numeric(sqrt(colSums(t(C)*(SIGbu%*%t(C))) + (C%*%MUbu)^2))
    
    ## Regression coefficients
    Zbar <- Diagonal(n,(wgt*0.5/xi)*tanh(0.5*xi))
    SIGbu <- solve(bdiag((1/sig2b)*Diagonal(p), 
                         as.numeric((Au+r/2)/(Bsig2u))*Diagonal(r)) +
                     t(C)%*%Zbar%*%C)
    MUbu <-  SIGbu %*% t(C) %*% as.numeric(wgt*(Y - 0.5))
    
    
    
    ## RE Variance
    Bsig2u <- Bu + 0.5*(t(MUbu[-c(1:p)])%*%(MUbu[-c(1:p)]) + sum(diag(SIGbu[-c(1:p),-c(1:p)])))
    
    PPrec <- bdiag((1/sig2b)*Diagonal(p), 
                   as.numeric((Au+r/2)/(Bsig2u))*Diagonal(r))
    
    ## Check for convergence
    checkNew  <- 0.5*(p+r) + 0.5*determinant(SIGbu, logarithm = T)$modulus  - 
      0.5*t(MUbu)%*%PPrec%*%(MUbu) + sum(wgt*(Y-0.5)*as.numeric(C%*%MUbu) +wgt*log(sigmoid(xi)) - 0.5*wgt*xi) - 
      0.5*sum(diag(PPrec %*% SIGbu)) - log(Bsig2u)
    checkVec <- c(checkVec, as.numeric(checkNew))
    
    
    if (as.logical(abs(checkOld - checkNew) < eps)) break
    iter <- iter + 1
    checkOld <- checkNew
  }
  return(list(iter=iter, BU=list(type="Normal", mean=MUbu, Cov=SIGbu, p=p, r=r),
              Bsig2u=list(type="Inv. Gamma", A=Au + r/2, B=Bsig2u), Elbo=checkVec))
}