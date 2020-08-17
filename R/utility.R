linkVStructures <- function(graph){
  result = graph
  falseind = c(NA, NA)
  
  for(i in 1:(ncol(graph)-1)){
    atrib.to = which(graph[,i] == 1)
    for(j in (i+1):ncol(graph)){
      if(!(j %in% atrib.to)){
        atrib.to.j = which(graph[,j] == 1)
        if(any(atrib.to %in% atrib.to.j)){
          result[i,j] = 1
        }
      } else {
        atrib.to.j = which(graph[,j] == 1)
        if(any(atrib.to %in% atrib.to.j)){
          falseind = rbind(falseind, c(i, j))
        }
      }
    }
  }
  
  ind = which((graph - result) != 0, arr.ind=T)
  
  output=list()
  output[["linkedGraph"]] = result
  output[["vstruc"]] = ind
  output[["falsevstruc"]] = falseind
  
  return(output)
}


convert <- function(res, G){
  for(i in 1:nrow(res)){
    res.row = unlist(res[i,])
    
    if(!is.na(res.row[1]) & !is.na(res.row[2])){
      G[res.row[1], res.row[2]] = res.row[3]
      G[res.row[2], res.row[1]] = res.row[3]
    }
  }
  
  return(G)
}

rama <- function(n){
  return(n*log(n) - n + (log(n*(1+4*n*(1+2*n))))/6 + log(pi)/2)
}

choose.approx <- function(n,k){
  return(rama(n) - rama(k) - rama(n-k))
}

hdist <- function(mat){
  if(class(mat)!='matrix'){
    mat = as.matrix(mat)
  }
  return(sum(abs((mat!=0 | t(mat)!=0) - ((graph!=0 | t(graph)!=0)))))
}

# For Step 1: Neighborhood Selection

#-----------------------------------
# Functions for numerical optimization. 
#   - Input details can be found under penZINB
#-----------------------------------

### Function to update likelihood
updateObsLik <- function(X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau){
  # ############# Update mean estimates
  linkobj.z = binomial()
  linkobj.count = negative.binomial(theta = theta)

  etaz = X%*%gammas + offsetz
  pij = linkobj.z$linkinv(etaz)

  eta = X%*%betas + offsetx
  mu = linkobj.count$linkinv(eta)

  ############# Calculate obs log likelihood
  loglik0 <- log( pij + exp( log(1-pij) + 
        suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1-pij) + 
          suppressWarnings(dnbinom(y, size = theta, mu = mu, log = TRUE))

  loglik.obs = sum(loglik0[Y0]) + sum(loglik1[-Y0]) 
  
  return(loglik.obs)
}

  # Wrapper functions around updateObsLik to make it suitable for deriving numerical hessian/score
updateObsLik.beta <- function(betas, X, gammas, offsetz, offsetx, y, Y0, theta, lambda, tau){
  return(updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau))
}

updateObsLik.gamma <- function(gammas, X, offsetz, betas, offsetx, y, Y0, theta, lambda, tau){
  return(updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau))
}
  
  # Wrapper function around updateLik to make it suitable for updating 1 set of beta/gamma at a time
optimLikFunc <- function(input, k, X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau, betaWeight, gammaWeight){
  betas[k] = input[1]
  gammas[k] = input[2]
  pen = sum(lambda*log(abs(betaWeight[-1]*betas[-1])+abs(gammaWeight[-1]*gammas[-1])+tau))

  loglik = updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, theta, 
    lambda, tau) - pen
  return(loglik)
}

  # Wrapper function around updateLik to make it suitable for updating theta
optimThetaFunc <- function(input, X, gammas, offsetz, betas, offsetx, y, Y0, lambda, tau, betaWeight, gammaWeight){
  theta = input
  pen = sum(lambda*log(abs(betaWeight[-1]*betas[-1])+abs(gammaWeight[-1]*gammas[-1])+tau))
  return(updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau) - pen)
}

  # Gradient function for optim
gradNegBin <- function(k, X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau, betaWeight, gammaWeight, type){
  Y1 <- y != 0
  Y0 <- y == 0

  ## count mean
  eta <- as.vector(X %*% betas + offsetx)
  mu <- negative.binomial(theta)$linkinv(eta)
  ## binary mean
  etaz <- as.vector(X %*% gammas + offsetz)
  muz <- binomial()$linkinv(etaz)

  ## densities at 0
  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- muz * (y==0) + exp(log(1 - muz) + clogdens0)

  if(type == 'coef'){
    wres_count <- ifelse(Y1, y - mu * (y + theta)/(mu + theta), -exp(-log(dens0) +
      log(1 - muz) + clogdens0 + log(theta) - log(mu + theta) + log(mu)))
    wres_count_pen <- lambda*betaWeight*sign(betas)/(betaWeight*abs(betas) + gammaWeight*abs(gammas))
    wres_count_pen[1] = 0

    wres_zero <- ifelse(Y1, -1/(1-muz) * binomial()$mu.eta(etaz),
      (binomial()$mu.eta(etaz) - exp(clogdens0) * binomial()$mu.eta(etaz))/dens0)
    wres_zero_pen <- lambda*gammaWeight*sign(gammas)/(betaWeight*abs(betas) + gammaWeight*abs(gammas))
    wres_zero_pen[1] = 0

    return(c(sum(wres_count*X[,k]) - wres_count_pen[k], sum(wres_zero*X[,k]) - wres_zero_pen[k]))
  } else {
    wres_theta <- theta * ifelse(Y1, digamma(y + theta) - digamma(theta) +
      log(theta) - log(mu + theta) + 1 - (y + theta)/(mu + theta),
      exp(-log(dens0) + log(1 - muz) + clogdens0) *
      (log(theta) - log(mu + theta) + 1 - theta/(mu + theta)))
    return(sum(wres_theta))
  }
}

  # wrapper for coefficient gradient:
gradCoef <- function(input, k, X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau, betaWeight, gammaWeight){
    betas[k] = input[1]
    gammas[k] = input[2]
    return(gradNegBin(k, X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau, betaWeight, gammaWeight, type = "coef"))
}
  # wrapper for theta gradient
gradTheta <- function(input, k, X, gammas, offsetz, betas, offsetx, y, Y0, lambda, tau, betaWeight, gammaWeight){
    theta = input
    return(gradNegBin(k, X, gammas, offsetz, betas, offsetx, y, Y0, theta, lambda, tau, betaWeight, gammaWeight, type = "theta"))
}

### Function to update coefficient for 1 gene (k) using optim
updateOptim.k <- function(y, betas, gammas, theta, k, offsetx, offsetz, lambda, tau, X, theta.st, betaWeight, gammaWeight, useBFGS){
  
  ############################ Intermediate variables
  M = ncol(X)
  N = length(y)
  Y0 <- which(y <= 0)
  Y1 <- which(y > 0)

  ############################ Initialize 
  b = betas[k]
  g = gammas[k]

  ############################ Optimize
  if(useBFGS == 1){
    result <- try(optim(c(b,g), optimLikFunc, gr = gradCoef, k = k, X = X, gammas = gammas, 
                offsetz = offsetz, betas = betas, offsetx = offsetx, y = y, 
                Y0 = Y0, theta = theta, lambda = lambda, tau = tau, 
                betaWeight = betaWeight, gammaWeight = gammaWeight,
                control = list(fnscale = -1), method = 'BFGS'), silent = TRUE)
    if(inherits(result, 'try-error')){
      result <- try(optim(c(b,g), optimLikFunc, k = k, X = X, gammas = gammas, 
                 offsetz = offsetz, betas = betas, offsetx = offsetx, y = y, 
                 Y0 = Y0, theta = theta, lambda = lambda, tau = tau, 
                 betaWeight = betaWeight, gammaWeight = gammaWeight,
                 control = list(fnscale = -1)), silent = TRUE)
      useBFGS = 0
    }
  } else {
   result <- try(optim(c(b,g), optimLikFunc, k = k, X = X, gammas = gammas, 
                 offsetz = offsetz, betas = betas, offsetx = offsetx, y = y, 
                 Y0 = Y0, theta = theta, lambda = lambda, tau = tau, 
                 betaWeight = betaWeight, gammaWeight = gammaWeight,
                 control = list(fnscale = -1)), silent = TRUE) 
  }

  if(inherits(result, 'try-error')){
    result = list(par = c(b,g), convergence = 0)
  }

  return(c(result$par, result$convergence, useBFGS))
}

updateTheta <- function(betas, gammas, offsetx, offsetz, lambda, tau, theta.st, theta, Y0, X, zGy, y, model = 1){

  ############################ Calculate current mean estimates
  linkobj.z = binomial()
  linkobj.count = negative.binomial(theta = theta)

  if(model == 1){
    etaz = X%*%gammas + offsetz
    pij = linkobj.z$linkinv(etaz)   
  }

  if(model == 2){
    zGy = rep(0,length(y))
  }
  
  eta = X%*%betas + offsetx
  mu = linkobj.count$linkinv(eta)

  ############################ Initialize
  thetamethod = 0

  theta.s = theta
  # if(is.null(theta.st)){
    theta <- try(theta_ml(theta.s, y, mu, weights = (1-zGy)), silent = TRUE)
    if(inherits(theta, "try-error")){
      theta = theta.s
    }
  # } else {
  #  theta = theta.st
  #}

  # if(is.null(theta.st)){
  #   theta.temp = theta.ml(y, mu, weights = (1-zGy))
  #   if(is.null(attr(theta.temp, "warn"))){
  #     theta = theta.temp
  #     thetamethod = 1
  #   } else {
  #     theta.temp = theta.mm(y, mu, dfr = sum(1-zGy) - sum((1-zGy)*(betas != 0)) - sum((1-zGy)*(gammas != 0))-2, weights = 1-zGy)
  #     if(is.null(attr(theta.temp, "warn"))){
  #       theta = theta.temp
  #       thetamethod = 2
  #     }
  #   }        
  # } else {
  #   theta = theta.st
  # }

  return(c(max(0.1, min(100,theta)), thetamethod))
}

#-----------------------------------
# Hessian functions of penalized ZINB (used for penalty weighting)
#-----------------------------------

hessBeta <- function(X, gammas, offsetz, betas, offsetx, y, theta){
  phi = 1/theta
  eta = as.vector(X%*%betas + offsetx)
  mu = negative.binomial(theta)$linkinv(eta)

  etaz = as.vector(X%*%gammas + offsetz)
  pij = binomial()$linkinv(etaz)

  t1 = (1+phi*mu)
  tt = (y==0)*(((phi+1)*(t1^(-theta-2))*pij + phi*(t1^(-2*theta - 2))*(1-pij))*(mu^2)*(1-pij))/
        (pij+(1-pij)*(t1^(-theta)))^2 - 
       (y!=0)*(mu*(1+phi*y))/((1+phi*mu)^2)
 
  return(-colSums(sweep(X^2,MARGIN=1,tt,`*`)))
}

hessGamma <- function(X, gammas, offsetz, betas, offsetx, y, theta){
  phi = 1/theta
  eta = as.vector(X%*%betas + offsetx)
  mu = negative.binomial(theta)$linkinv(eta)

  etaz = as.vector(X%*%gammas + offsetz)
  pij = binomial()$linkinv(etaz)

  tt = (y==0)*(exp(etaz)*(1+phi*mu)^(theta))/(exp(etaz)*(1+phi*mu)^(theta)+1)^2 - (exp(etaz)/(1+exp(etaz))^2)
  
  return(-colSums(sweep(X^2,MARGIN=1,tt,`*`)))
}

# For step 2: Testing Conditional Independence ------
# Purpose : search for connected components of the vertices set x
# Input ------------------------------------------------------------------------
#     - G : undirected adjacency matrix (edges =TRUE, gaps = FALSE)
#     - x : set of vertices
# Output------------------------------------------------------------------------ 
#     - the union of connected components of all elements of x (x included)
# ------------------------------------------------------------------------------

connectedComp <- function(G,x) {
  if (!identical(dim(G)[1],dim(G)[2])) stop('G must be square matrix')
  if (!all(G==t(G))) stop('G must be symmetric')
  if (sum(diag(G))!=0) stop('diag(G) must be all zero')
  if (!is.vector(x)) stop('x must be a vector')
  p = nrow(G)
  id = seq_len(p)
  
  num_addcomp = 1L
  Ccomp = upx = x # connected components to be updated
  while (num_addcomp != 0) {
    adjmat = as.matrix(G[,upx])
    candx = id[rowSums(adjmat)!=0]
    upx = candx[!candx%in%Ccomp]
    #print(upx)
    num_addcomp = length(upx)
    Ccomp = unique(c(Ccomp,upx))
  }
  return((id %in% Ccomp))
}

