#' Compare the estimated graph with the true graph
#' 
#' Both graphs should be undirected graphs. 
#'
#' @param gl the estimated graph
#' @param gt the true graph
#' @export
#' @return A vector of three components: true positive rate, false positive 
#' rate, true discovery rate. 
compareGraphs <- function(gl, gt){
  
  p <- ncol(gl)
  ac <- p^2 - p # the number of all connections
  
  # tpr
  tpr <- sum(gl[which(gt == 1)]) / sum(gt)
  
  # fpr
  fpr <- sum(gl[which(gt == 0)]) / (ac - sum(gt))
  
  # tdr
  tdr <- sum(gl[which(gt == 1)]) / sum(gl)
  
  c(tpr = tpr, fpr = fpr, tdr = tdr)
  
}