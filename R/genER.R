#' Generate data using the Erdos and Renyi (ER) model 
#'
#' Given the number of vertices, the probability that any two vertices connect,
#' and the sample size, generate the true graph structure and the sample data.
#'
#' @details The true graph structure is generated such that the probability of 
#' having an edge for any pair of vertices is \eqn{p_E = d_0/p}, where \eqn{d_0}
#' is average number of edges per node. Then, all graphs with p vertices and
#' d edges have probability of \eqn{p_E^d(1-p_E)^{\choose(d, 2) -d}} to be
#' generated. After we generate the adjacency matrix \eqn{A}, the sample 
#' data are simulated from the multivariate normal distribution with
#' mean 0 and covariance \deqn{ (I - B)^{-1}(I-B^T)^{-1} } where \eqn{B = } the
#' effect size \eqn{\times} the lower triangular matrix of \eqn{A}.
#' 
#' @param p number of vertices
#' @param pE probability that two vertices will connect
#' @param n sample size for data generation
#' @param effS the effect size, must be a single number - - default is 
#' \code{NULL} which generates random coefficients for each covariate from a 
#' Uniform(0,1) distribution.
#' 
#' @export
#' @return A list with the following components:
#' \describe{
#'   \item{A}{adjacency matrix of the simulated graph}
#'   \item{X}{simulated sample data}
#'   \item{G}{simulated graph as a graphNEL object}
#' } 
#' @examples
#' genER(10, 0.1, 20, 0.5)
genER <- function(p, pE, n, effS = NULL) {
  stopifnot(pE>0,pE<1)
  A = matrix(0,p,p)
  w = which(lower.tri(A))
  
  if(is.null(effS)){
    effS = runif(n=length(w),min = 0, max = 1)
  }
  
  A[w] = rbinom(n=length(w),size=1,prob=pE)*effS
  Sigh = solve(diag(p) - A)
  X = rmvnorm(n=n,sigma = Sigh %*% t(Sigh))
  
  # shuffling
  id.shuffle = sample(1:p,p,replace=FALSE)
  A = A[id.shuffle,id.shuffle]
  X = X[,id.shuffle]
  return(list(A=A,X=X, G=as(t(A),"graphNEL")))
}

