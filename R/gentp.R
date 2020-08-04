#' Generate tuning parameters for variable selection in the zero-inflated 
#' negative binomial regression
#'
#' Use a dynamic tuning parameter algorithm to generate lambads and taus 
#' in penalized likelihood estimation of the ZINB model. 
#' 
#' @inheritParams penZINB
#' @param offsetx,offsetz Two vector of observations that are included in the 
#' linear predictors of negative binomial regression and logistic regression 
#' respectively. Default are \code{NULL}. 
#' @param betaweight,gammaweight Weights of the coefficients for the 
#' penalization in negative binomial regression and logistic regression 
#' respectively. Default are 1.
#' @export
#' @return A matrix with two columns. Each row represents a pair of tuning
#' parameters: lambda and tau. 
gentp <- function(y, X, nlambda = 30, ntau = 5, unpenalizedx = NULL, 
                  unpenalizedz = NULL, offsetx = NULL, offsetz = NULL,
                  betaweight = 1, gammaweight = 1,
                  pfactor = 1e-2, penType = 1) {
  if(.Platform$OS.type == "unix"){
    sinkNul = "/dev/null"
  } else {
    sinkNul = 'NUL'
  }
  
  XR = scale(X, center = FALSE) 
  N = length(y)
  # Create tuning parameter grid if not given
  
  if(is.null(offsetx)){
    offsetx = rep(0, length(y))
  } 
  
  if(is.null(offsetz)){
    offsetz = rep(0, length(y))
  }
  
  
  fit <- try(zeroinfl(y~1 + offset(offsetx) | 1 + offset(offsetz), 
                      dist="negbin", link = "logit"), silent = TRUE)
  
  if(inherits(fit, 'try-error')){
    capture.output(fit <- penZINB(y, 1, unpenalizedx = unpenalizedx, 
                                  unpenalizedz = unpenalizedz,
                                  lambdas=0, taus=.1, maxIT=100, 
                                  maxOptimIT = 0, eps = 2e-3, theta = NULL, 
                                  warmStart = 'cond', pfactor = 1e-7, 
                                  oneTheta = TRUE, convType = 1, 
                                  start = 'jumpstart', order = FALSE), 
                   file = sinkNul)
    
    pi. <- binomial()$linkinv(fit[[1]]$gammas)
    mu. <- negative.binomial(1)$linkinv(fit[[1]]$betas)
  } else {
    pi. <- binomial()$linkinv(fit$coefficient$zero)
    mu. <- negative.binomial(1)$linkinv(fit$coefficient$count)      
  }  
  
  w <- pi.*(1-pi.)
  r = (y - pi.)/w
  
  l1.g <- 2 * abs( crossprod(XR,w*r) / gammaweight )
  
  w <- (mu./(1 + mu./1))
  r = (y-mu.)/mu.
  
  l1.b <- 2 * abs(crossprod(XR,w*r) / betaweight )
  l1.max = max(l1.g, l1.b)
  
  
  Thresholds <- c(exp(seq(log(l1.max),log(pfactor*l1.max),len=nlambda)))
  # Thresholds <- c(seq(l1.max,pfactor*l1.max,len=nlambda))   
  
  if(penType == 1){
    maxTau = min(l1.max/N, 0.1)
    tauset = c(exp(seq(log(1e-6),log(maxTau),len=ntau)))
    
    lambda_max <- maxTau * l1.max
    tmp <- penZINB(y, XR, lambdas = lambda_max, taus = max(tauset) )
    while(length(tmp[[1]]$betas.w) > 1 | length(tmp[[1]]$gammas.w) > 1){
      lambda_max <- lambda_max + 2
      tmp <- penZINB(y, XR, lambdas = lambda_max, taus = max(tauset) )
    }
    
    l1.max <- lambda_max / maxTau
    Thresholds <- c(exp(seq(log(l1.max),log(pfactor*l1.max),len=nlambda)))
    
    taus = rep(tauset, each = nlambda)
    lambdas = taus*Thresholds  
  } else if(penType == 2){
    taus = rep(1, nlambda)
    lambdas = Thresholds
  }
  
  # Always include case of no penalty 
  if(is.na(any(lambdas/taus == 0) == FALSE)){
    lambdas = lambdas[!is.na(lambdas) & !is.na(taus)]
    taus = taus[!is.na(lambdas) & !is.na(taus)]
  }
  if(any(lambdas/taus == 0) == FALSE){
    if(!is.null(lambdas)){
      if(length(lambdas > 1)){
        lambdas = c(lambdas, 0)        
      }
    }
    
    if(!is.null(taus)){
      if(length(taus) > 1){
        taus = c(taus, 0.1)        
      }
    }    
  }
  
  # sort tuning parameters from smallest to largest for warm start 
  
  tp.ratio = lambdas/taus
  tp.idx = sort.int(tp.ratio, index.return = TRUE)$ix
  
  lambdas = lambdas[tp.idx]
  taus = taus[tp.idx]
  
  return(cbind(lambdas,taus))
}





