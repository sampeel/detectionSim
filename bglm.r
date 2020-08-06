gearglm <- function(a.formula, b.formula, data, g.formula=NULL, area=1, control=glm.control(),
                    method=c("model","bglmQR","bglmQRC","bglmC","bglm")) {
  
  method <- match.arg(method)
  
  ## Construct A matrix
  mf <- model.frame(a.formula,data,na.action=na.pass)
  A <- model.matrix(a.formula,mf)
  Y <- model.response(mf)
  if ( is.null(ncol(Y)) ) Y <- as.matrix(Y)
  
  ## Construct B matrix
  mf <- model.frame(b.formula,data,na.action=na.pass)
  B <- model.matrix(b.formula,mf)
  isP <- model.response(mf)
  
  ## Construct the gear matrix G, if necessary.
  if ( is.null(g.formula) ) {
    G <- NULL
  } else {
    mf <- model.frame(g.formula,data,na.action=na.pass)
    G <- model.matrix(g.formula,mf)
    for(j in seq_len(ncol(G))) {
      G[,j] <- ifelse(isP,0,G[,j])
    }
  }
  
  ## Adjust intercepts
  as <- attr(A,"assign")
  A <- cbind(A,`(Bias)`=ifelse(isP,1,0))
  attr(A,"assign") <- c(as,max(as)+1)
  if(attr(B,"assign")[1]==0) {
    as <- attr(B,"assign")
    B <- B[,-1,drop=FALSE]
    attr(B,"assign") <- as[-1]
  }
  if ( ! is.null(g.formula) ) {
    if(attr(G,"assign")[1]==0) {
      as <- attr(G,"assign")
      G <- G[,-1,drop=FALSE]
      attr(G,"assign") <- as[-1]
    }
  
    ## Augment A with gear effects
    A <- cbind(A,G)
  }
  
  family <- fithian(isP)
  offset <- log(area)
  
  
  # Block IRLS
  fit <- switch(method, 
          model=list(A=A,B=B,Y=Y,offset=offset,family=family,control=control),
          bglmQR=bglmQR.fit(A,B,Y,offset=offset,family=family,control=control,method="bglm"),
          bglmQRC=bglmQR.fit(A,B,Y,offset=offset,family=family,control=control,method="bglmC"),
          bglmC=bglmC.fit(A,B,Y,offset=offset,family=family,control=control),
          bglm=bglm.fit(A,B,Y,offset=offset,family=family,control=control))
  
  # Return value.
  return(fit)
  
}

#-----------------------------------------------------------------------------------------

bglm.fit <- function(A, B, Y, weights=1, offset=0, family, control=glm.control()) {
  
  invperm <- function(p) {
    p[p] <- seq_along(p)
    p
  }
  
  nobs <- nrow(Y)
  mY <- ncol(Y)
  mA <- ncol(A)
  mB <- ncol(B)
  mbeta <- mY*mA+mB
  
  yNA <- is.na(Y)
  
  
  ## Initialize
  beta <- double(mbeta)
  js <- (mY*mA)+seq_len(mB)
  deviance.old <- Inf
  converged <- FALSE
  
  ## IRLS iterations
  for(iter in seq_len(control$maxit)) {
    
    ## Construct the (weighted) normal equation XTWX = XTWZ
    XTWX <- matrix(0,mbeta,mbeta)
    XTWZ <- double(mbeta)
    for(k in seq_len(mY)) {
      is <- ((k-1)*mA)+seq_len(mA)
      
      ## IRLS weights and adjusted response
      y <- Y[,k]
      if(iter>1) {
        eta <-  A%*%beta[is] + B%*%beta[js] + offset
      } else {
        eta <- local({
          weights <- if(is.matrix(weights)) weights[,k] else weights
          weights <- rep(weights,length.out=nobs)
          weights[yNA[,k]] <- 0
          eval(family$initialize)
          family$linkfun(mustart)
        })
      }
      mu <-  family$linkinv(eta)
      mu.eta <- family$mu.eta(eta)
      varmu <- family$variance(mu)
      z <- (eta-offset) + (y-mu)/mu.eta
      w <- if(is.matrix(weights)) weights[,k] else weights
      w <- w*(mu.eta^2)/varmu
      w[yNA[,k]] <- 0
      
      ## Contribution to XTWX 
      XTWX[is,is] <- crossprod(A,w*A)
      XTWX[js,is] <- t(XTWX[is,js] <- crossprod(A,w*B))
      XTWX[js,js] <- XTWX[js,js]+crossprod(B,w*B)
      
      ## Contribution to XTWZ
      wz <- w*z
      wz[yNA[,k]] <- 0
      XTWZ[is] <- crossprod(A,wz)
      XTWZ[js] <- XTWZ[js]+crossprod(B,wz)
    }
    beta.old <- beta
    
    ## Solve by Cholesky decomposition
    U <- suppressWarnings(chol(XTWX,pivot=TRUE))
    rank <- attr(U,"rank")
    p <- attr(U,"pivot")
    if(rank < mbeta) {
      U <- U[seq_len(rank),seq_len(rank)]  
      p <- p[seq_len(rank)]
    }
    beta0 <- backsolve(U,backsolve(U,XTWZ[p],transpose=TRUE))
    beta <- double(mbeta)
    beta[p] <- beta0
    #converged <- sqrt(crossprod(beta-beta.old)) < control$epsilon
    #if(converged) break
    
    ## Compute the deviance for this iteration.
    deviance <- 0
    Mu <- Y
    for(k in seq_len(mY)) {
      is <- ((k-1)*mA)+seq_len(mA)
      eta <-  A%*%beta[is] + B%*%beta[js] + offset
      Mu[,k] <-  family$linkinv(eta)
      wts <- if(is.matrix(weights)) weights[,k] else rep(weights,length.out=length(y))
      devres <- family$dev.resid(Y[,k],Mu[,k],wts)
      deviance <- deviance + sum(devres[!yNA[,k]])
    }
    
    ## Has method converged?
    if ( (abs(deviance - deviance.old)/(0.1 + abs(deviance))) < control$epsilon ) {
      converged <- TRUE
      break
    } else {
      deviance.old <- deviance
    }
  }
  
  beta <- rep(NA,mbeta)
  beta[p] <- beta0
  V <- matrix(NA,mbeta,mbeta)
  V[p,p] <- chol2inv(U)
  
  list(coefficients=beta,
       cov=V,
       deviance=deviance,
       fitted.values=Mu,
       prior.weights=weights,
       converged=converged,
       iter=iter)
}

#-----------------------------------------------------------------------------------------
bglmQR.fit <- function(A, B, Y, weights=1, offset=0, family, control=glm.control(),method=c("bglm","bglmC")) {
  method <- match.arg(method)
  
  nms <- if(!is.null(colnames(A)) && !is.null(colnames(B))) {
    sp <- if(!is.null(colnames(Y))) colnames(Y) else paste0("sp",seq_len(ncol(Y)))
    c(outer(colnames(A),sp,function(x,s) paste(s,x,sep=":")),colnames(B))
  }
  
  ## Compute the QR decompositions of A, B.
  A.qr <- qr(A)
  A.rank <- A.qr$rank
  A.Q <- qr.Q(A.qr)
  if(A.rank < ncol(A.Q))
    A.Q <- A.Q[,seq_len(A.rank)]
  
  AB.R <- crossprod(A.Q,B)
  
  B.qr <- qr(B-A.Q%*%AB.R)
  B.rank <- B.qr$rank
  B.Q <- qr.Q(B.qr)
  if(B.rank < ncol(B.Q)) {
    B.Q <- B.Q[,seq_len(B.rank)]
    AB.R <- AB.R[,seq_len(B.rank)]
  }
  
  Y.n <- ncol(Y)
  
  ## Fit with orthogonalized A, B
  fit <- switch(method,
                bglm=bglm.fit(A.Q,B.Q,Y,weights=weights,offset=offset,family=family,control=control),
                bglmC=bglmC.fit(A.Q,B.Q,Y,weights=weights,offset=offset,family=family,control=control))
  
  ## Rescale by R
  A.Rinv <- solve(qr.R(A.qr)[seq_len(A.rank),seq_len(A.rank)])
  B.Rinv <- solve(qr.R(B.qr)[seq_len(B.rank),seq_len(B.rank)])
  AB.Rinv <- -A.Rinv %*% AB.R %*% B.Rinv
  Rinv <- matrix(0,Y.n*A.rank+B.rank,Y.n*A.rank+B.rank)
  js <- seq.int(Y.n*A.rank+1,length.out=B.rank)
  Rinv[js,js] <- B.Rinv
  for(k in seq_len(Y.n)) {
    is <- seq.int((k-1)*A.rank+1,length.out=A.rank)
    Rinv[is,is] <- A.Rinv
    Rinv[is,js] <- AB.Rinv
    
  }
  fit$coefficients <- drop(Rinv%*%fit$coefficients)
  fit$cov <- Rinv%*%tcrossprod(fit$cov,Rinv)
  
  ## Insert NA for rank deficient A,B 
  A.is <- matrix(NA,ncol(A),Y.n)
  A.is[A.qr$pivot[seq_len(A.rank)],] <- seq_len(Y.n*A.rank)
  B.is <- rep(NA,ncol(B))
  B.is[B.qr$pivot[seq_len(B.rank)]] <- seq.int(Y.n*A.rank+1,length.out=B.rank)
  is <- c(A.is,B.is)
  fit$coefficients <- fit$coefficients[is]
  fit$cov <- fit$cov[is,is]
  
  names(fit$coefficients) <- nms
  dimnames(fit$cov) <- list(nms,nms)
  
  #if ( any(is.na(fit$coefficients)) ) fit$converged <- FALSE   ### Temporary line to suss out rank deficient issue.
  fit
}


#-----------------------------------------------------------------------------------------

bglmC.fit <- function(A, B, Y, weights=1, offset=0, family, control=glm.control()) {
  
  invperm <- function(p) {
    p[p] <- seq_along(p)
    p
  }
  
  nobs <- nrow(Y)
  mY <- ncol(Y)
  mA <- ncol(A)
  mB <- ncol(B)
  mbeta <- mY*mA+mB
  
  yNA <- is.na(Y)
  
  
  ## Initialize
  beta <- double(mbeta)
  js <- (mY*mA)+seq_len(mB)
  deviance.old <- Inf
  
  ## IRLS iterations
  for(iter in seq_len(control$maxit)) {
    
    ## Construct the (weighted) normal equation XTWX = XTWZ
    U <- matrix(0,mbeta,mbeta)
    XTWZ <- double(mbeta)
    for(k in seq_len(mY)) {
      is <- ((k-1)*mA)+seq_len(mA)
      
      ## IRLS weights and adjusted response
      y <- Y[,k]
      if(iter>1) {
        eta <-  A%*%beta[is] + B%*%beta[js] + offset
      } else {
        eta <- local({
          weights <- if(is.matrix(weights)) weights[,k] else weights
          weights <- rep(weights,length.out=nobs)
          weights[yNA[,k]] <- 0
          eval(family$initialize)
          family$linkfun(mustart)
        })
      }
      mu <-  family$linkinv(eta)
      mu.eta <- family$mu.eta(eta)
      varmu <- family$variance(mu)
      z <- (eta-offset) + (y-mu)/mu.eta
      w <- if(is.matrix(weights)) weights[,k] else weights
      w <- w*(mu.eta^2)/varmu
      w[yNA[,k]] <- 0
      
      ## Contribution to U = chol(XTWX) 
      U[is,js] <- backsolve((U[is,is] <- chol(crossprod(A,w*A))),crossprod(A,w*B),transpose=T)
      U[js,js] <- U[js,js] + crossprod(B,w*B) - crossprod(U[is,js])
      
      ## Contribution to XTWZ
      wz <- w*z
      wz[yNA[,k]] <- 0
      XTWZ[is] <- crossprod(A,wz)
      XTWZ[js] <- XTWZ[js]+crossprod(B,wz)
    }
    U[js,js] <- chol(U[js,js])
    beta.old <- beta
    
    ## Solve by Cholesky decomposition
    beta <- backsolve(U,backsolve(U,XTWZ,transpose=TRUE))
  #  converged <- sqrt(crossprod(beta-beta.old)) < control$epsilon
  #  if(converged) break
    
    
    ## Compute the deviance for this iteration
    deviance <- 0
    Mu <- Y
    for(k in seq_len(mY)) {
      is <- ((k-1)*mA)+seq_len(mA)
      eta <-  A%*%beta[is] + B%*%beta[js] + offset
      Mu[,k] <-  family$linkinv(eta)
      wts <- if(is.matrix(weights)) weights[,k] else rep(weights,length.out=length(y))
      devres <- family$dev.resid(Y[,k],Mu[,k],wts)
      deviance <- deviance + sum(devres[!yNA[,k]])
    }
    
    ## Has method converged?
    if ( (abs(deviance - deviance.old)/(0.1 + abs(deviance))) < control$epsilon ) {
      converged <- TRUE
      break
    } else {
      deviance.old <- deviance
    }
  }
  
  list(coefficients=beta,
       cov=chol2inv(U),
       deviance=deviance,
       fitted.values=Mu,
       prior.weights=weights,
       converged=converged,
       iter=iter)
}

#-----------------------------------------------------------------------------------------

fithian <- function(isP) {
  N <- length(isP)
  isB <- !isP
  P <- poisson(link=log)
  B <- binomial(link=cloglog)
  
  linkfun <- function(mu) {
    R <- double(N)
    R[isP] <- P$linkfun(mu[isP])
    R[isB] <- B$linkfun(mu[isB])
    R
  }
  linkinv <- function(eta) {
    R <- double(N)
    R[isP] <- P$linkinv(eta[isP])
    R[isB] <- B$linkinv(eta[isB])
    R
  }
  mu.eta <- function(eta) {
    R <- double(N)
    R[isP] <- P$mu.eta(eta[isP])
    R[isB] <- B$mu.eta(eta[isB])
    R
  }
  valideta <- function(eta) TRUE
  variance <- function(mu) {
    R <- double(N)
    R[isP] <- P$variance(mu[isP])
    R[isB] <- B$variance(mu[isB])
    R
  }
  validmu <- function(mu) all(is.finite(mu)) && all(mu>0 & (isP | mu<1))
  dev.resids <- function(y,mu,wt) {
    R <- double(N)
    R[isP] <- P$dev.resids(y[isP],mu[isP],wt[isP])
    R[isB] <- B$dev.resids(y[isB],mu[isB],wt[isB])
    R
  }
  
  aic <- function(y,n,mu,wt,dev) {
    P$aic(y[isP],n[isB],mu[isP],wt[isP],dev)+B$aic(y[isB],n[isB],mu[isB],wt[isB],dev)
  }
  
  ## Hard to get this right
  initialize <- substitute({
    n <- rep.int(1, nobs)
    mustart <- ifelse(isP,y+0.1,(y+0.5)/2)
  },list(isP=isP))
  
  structure(list(family = "fithian",
                 link = "fithian",
                 linkfun = linkfun,
                 linkinv = linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 aic = aic,
                 mu.eta = mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = valideta),
            class = "family")  
}

#-----------------------------------------------------------------------------------------