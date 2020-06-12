#' Perform the conditional independence test based on the linear model
#'
#' Check whether x and y are independent conditional on z in the linear model.
#' 
#' @param x,y variables for testing independence
#' @param z conditioning set of variables
#' @param alpha the level of the significance for test of independence
#' @export
#' @return returns the p-value of the test
cond_LM <- function(x, y, z, alpha){
  
  x = as.matrix(log(x + 1 / 6))
  y = as.matrix(log(y + 1 / 6))
  
  if(length(z) == 0){
    zwo = rep(1, nrow(x))
    zw = numeric(0)
  } else {
    zwo = zw = as.matrix(log(z+1/6))
  }
  
  prec <- glm(y ~ cbind(zw, x), family = "gaussian")
  
  precwo <- glm(y ~ zwo, family = "gaussian")
  
  res <- anova(prec, precwo, test = "LRT")
  
  pval = res$`Pr(>Chi)`[2]
  
  if(pval > alpha) {
    prec <- glm(x ~ cbind(zw, y), family = "gaussian")
    
    precwo <- glm(x ~ zwo, family = "gaussian")
    
    res <- anova(prec, precwo, test = "LRT")
    
    pval = res$`Pr(>Chi)`[2]
  }
  
  return(pval)
}





