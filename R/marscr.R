#' Perform the marginal screening prior to penalized regression 
#'
#' This function removes the covariates that are not even weakly associated with 
#' the response variable prior to penalized regression to improve the performance
#' of the variable selection and reduce the runtime.
#' 
#' @param data the sample data 
#' @param thres either a number between 0 and 1 or an positive integer. If it is 
#' a number between 0 and 1, the function selects the variables whose 
#' correlation with the response are greater than the number; 
#' if it is an integer, the function selects the top thres variables. 
#' Default is 0.1.  
#' @export
#' @return A square matrix. The ith row provides the most correlated variables 
#' with the ith variable. 
#' @examples 
#' data <- genER(10, 0.1, 20, 0.5)$X
#' gEstM <- marscr(data, thres = 5)
marscr <- function(data, thres = 0.1){
  
  cord <- cor(data)
  
  if(thres < 1){
    res <- (abs(cord) >= thres) * 1
    diag(res) <- 0
  }else{
    rank_data <- apply(-abs(cord), 2, rank)
    res <- (rank_data <= (thres + 1)) * 1
    diag(res) <- 0
  }
  
  res
  
}