#' Perform the conditional independence test based on the ZINB model
#'
#' Check whether x and y are independent conditional on z in the ZINB model.
#' 
#' @param x,y variables for testing independence
#' @param z conditioning set of variables
#' @param alpha the level of the significance for test of independence
#' @export
#' @return returns the p-value of the test
cond_ZINB <- function(x, y, z, alpha){
  
  x = as.matrix(x)
  y = as.matrix(y)
  
  if(length(z) == 0){
    zwo = 1
    zw = numeric(0)
  } else {
    zwo = zw = as.matrix(log(z+1/6))
  }
  
  prec = penZINB(y, (cbind(zw, log(x+1/6))), unpenalizedx = NULL, 
                unpenalizedz = NULL, 
                maxOptimIT = 0, 
                eps = 1e-5, theta = NULL, 
                warmStart = FALSE, lambdas = 0, 
                taus = 1e-6, pfactor = 1e-4, 
                oneTheta = FALSE, convType = 1, 
                start = 'jumpstart')
  
  precwo = penZINB(y, (zwo), unpenalizedx = NULL, 
                  unpenalizedz = NULL, 
                  maxOptimIT = 0, 
                  eps = 1e-5, theta = NULL, warmStart = FALSE, 
                  lambdas = 0, taus = 1e-6, pfactor = 1e-7, 
                  oneTheta = FALSE, convType = 1, 
                  start = 'jumpstart')
  
  loglikw = prec[[1]]$loglik.obs
  loglikwo = precwo[[1]]$loglik.obs
  
  stat = 2*(loglikw - loglikwo)
  pval = 1-pchisq(stat, df = 2)
  
  if(pval > alpha) {
    prec = penZINB(x, (cbind(zw, log(y+1/6))), unpenalizedx = NULL, 
                  unpenalizedz = NULL, maxOptimIT = 0, 
                  eps = 1e-5, theta = NULL, 
                  warmStart = FALSE, lambdas = 0, 
                  taus = 1e-6, pfactor = 1e-7, 
                  oneTheta = FALSE, convType = 1, 
                  start = 'jumpstart')
    
    precwo = penZINB(x, (zwo), unpenalizedx = NULL, 
                    unpenalizedz = NULL,  
                    maxOptimIT = 0, 
                    eps = 1e-5, theta = NULL, 
                    warmStart = FALSE, lambdas = 0, 
                    taus = 1e-6, pfactor = 1e-7, 
                    oneTheta = FALSE, convType = 1, 
                    start = 'jumpstart')
    
    loglikw = prec[[1]]$loglik.obs
    loglikwo = precwo[[1]]$loglik.obs
    
    stat = 2*(loglikw - loglikwo)
    pval = min(1-pchisq(stat, df = 2), pval)
  }
  
  return(pval)
}