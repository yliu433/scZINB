#' Create a graph based on the output from the nsZINB function
#' 
#' @param models a list of p models. The jth model is the output from the nsZINB 
#' function using the jth variable as the response and rest (p-1) 
#' variables as the covariates. 
#' @param filter the filter matrix prior to performing the variable selection. 
#' @param bic the criterion used to create the graph, could be "BIC", "extBIC"
#' or "extBICGG". Default is "extBIC". 
#' @export
#' @return returns the estimated moral graph. 
createGraph <- function(models, filter, bic = "extBIC"){
  
  p <- length(models)
  
  sapply(1:p, function(i){
    
    model <- models[[i]]
    
    model1 <- model[[which.min(sapply(model, "[[", paste(bic)))]]
    
    idx1 <- unique(model1$betas.w, model1$gammas.w)
    
    idx <- which(filter[, i] != 0)
    
    res <- rep(0, p)
    
    res[idx[idx1 - 1]] <- 1
    
    res
    
  })
  
}

