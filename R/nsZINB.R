#' Neighborhood selection for zero-inflated negative binomial regression
#'
#' Perform variable selection for ZINB regression via penalized maximum 
#' likelihood. 
#' 
#'
#' @param y zero-inflated count response
#' @param X covariate matrix. Intercept is added within the function. This 
#' could take '1' as the input which indicates an intercept-only model.
#' @param unpenalizedx,unpenalizedz Additional unpenalized covariates for 
#' negative binomial and logistic regression respectively. Default is 
#' \code{NULL}.
#' @param lambdas,taus specific tuning parameter values you want to run the 
#' model with. Default is \code{NULL} where the function will auto-generate
#' a tuning parameter search grid. If default is used, must have input for
#' nlambda and ntau.
#' @param nlambda,ntau number of unique lambda and tau values - default are 10 
#' and 3.
#' @param naPercent allowable percentage of observations with missing values - 
#' default is .4.
#' @param maxIT maximum number of EM iterations - default is 1000.
#' @param maxIT2 maximum number of iterations for updating the coefficients
#' in the regression model - default is 25.
#' @param track default is \code{NULL} (deactivated). Otherwise, it takes a 
#' single integer value which activates tracking mode for that tuning
#' parameter pair. See output change details below.
#' @param theta.st default is \code{NULL} (deactivated) where theta estimation 
#' is done using MLE. Otherwise, takes a single value for theta to hold 
#' constant for all estimation.
#' @param stepThrough default is \code{NULL} (deactivated). Otherwise needs to 
#' be a length 2 vector to activate debugging mode. The first number is
#' the theta iteration and second is the EM iteration to prompt
#' stepthrough debugger.
#' @param loud default is \code{NULL} (deactivated). Otherwise takes a 
#' positive integer x to announce at every xth iteration of EM and 
#' numerical optimization algorithm.
#' @param optimType options are "EM" and "optim". Default is "EM" which 
#' runs the EM algorithm prior to using BFGS optimization. "optim" skips the 
#' EM algorithm. 
#' @param warmStart default is TRUE, which keeps previous estimates as starting
#' points for estimation for the next tuning parameter. Other options are
#' 'cond', which resets the the starting point to the original starting point
#' when non-convergence happens. FALSE uses the same starting point for all tp.
#' @param irlsConv forces each estimate of beta and gamma to converge first if 
#' set to TRUE. Default is FALSE.
#' @param weightedPen default is TRUE. Weights the penalty using the Hessian. 
#' @param numericalDeriv default is FALSE. Calculates the Hessian numerically 
#' when set to TRUE. Otherwise, calculates it analytically.
#' @param pfactor default is 1e-2. The multiplier for the largest calculated 
#' penalty to determine smallest penalty value. Use in conjunction with
#' nlambda/ntau to control the granularity of the tp grid.
#' @param oneTheta default is FALSE (deactivated). If set to TRUE, only 
#' estimates theta once per tuning parameter pair.
#' @param maxOptimIT maximum number of iterations for numerical optimization 
#' (BFGS) after the EM algorithm. By default is set to 50. Convergence time
#' is long. 
#' @param errorDirec whenever the EM algorithm errors out, the function auto 
#' saves to this directory and moves on to the next tuning parameters. Default
#' is the R temporary directory. 
#' @param eps threshold for convergence for the EM algorithm.
#' @param convType manages the order of convergence within the EM algorithm. 
#' Options are 1 (default) and 2. Type 1 forces convergence of the binomial
#' and negative binomial parts together. Type 2 forces convergence of binomial
#' part first, then negative binomial part.
#' @param start default is \code{NULL} which sets starting coefficients values 
#' to 0. If set to 'jumpstart', then will estimate the starting coefficients
#' from penalized negative binomial estimation and logistic regression based on
#' the penalized library. Otherwise, can also take direct input for starting
#' values. Must be in the form of list(betas = v1, gammas = v2), where v1 and
#' v2 are vectors the length of the number of covariates in X.
#' @param order default is FALSE. If TRUE, then order of estimation is ordered 
#' by marginal correlation with response.
#' @param alwaysUpdateWeights default is FALSE. Updates the penalty weights 
#' after each convergence of the EM algorithm.
#' @param penType options are 1 (default) or 2. 1 is the group log penalty. 2 
#' is lasso.   
#' 
#' @details  If tracking, this function returns a nested list of all 
#' estimated with the following hierarchichy:
#' \itemize{
#'   \item Level 1 - TP (may not be included if only 1 tp pair is tracked);
#'    last two elements are lambda and tau,
#'   \item Level 2 - Theta Estimate (up to 20); last two elements are theta,
#'   \item Level 3 - EM Iteration (up tp max it), 
#' }  
#' with the following values: loglik, loglik.em, loglikZI, loglikNB, pen, 
#' betas, gammas. 
#' @return A list with each element corresponding to each tuning 
#' parameter pair. Each element contains the following components:
#' \describe{
#'   \item{X}{The design matrix used for calculations. Non-empty only for the 
#'   first tuning parameter pair.}
#'   \item{betas}{Non-zero beta coefficients corresponding to betas.w.}
#'   \item{gammas}{Non-zero gamma coefficients corresponding to gammas.w.}
#'   \item{loglik.obs}{Observed data log likelihood at convergence.}
#'   \item{pen}{Value of the penalty at convergence.}
#'   \item{theta.r}{Theta path.}
#'   \item{theta}{Theta estimate at convergence.}
#'   \item{BIC}{BIC value at convergence.}
#'   \item{extBIC}{Extended BIC value at convergence.}
#'   \item{extBICGG}{Extended BIC GG value at convergence.}
#'   \item{lambda, tau}{Tuning parameter pair used.}
#' }
#' @seealso \code{\link{gentp}} for generating tuning parameters of lambdas
#' and taus.
#' @export
nsZINB <- function(y, X, unpenalizedx = NULL, unpenalizedz=NULL, lambdas=NULL, 
                   taus=NULL, nlambda=10, ntau=3, naPercent =.4, maxIT = 1000, 
                   maxIT2 = 25, track = NULL, theta.st = NULL, 
                   stepThrough = NULL, optimType = "EM", 
                   loud = NULL, warmStart = TRUE, bicgamma = NULL,
                   irlsConv = FALSE, weightedPen = TRUE, 
                   numericalDeriv = FALSE, pfactor = 1e-2, 
                   oneTheta = FALSE, maxOptimIT = 50, 
                   errorDirec = tempdir(), eps = 1e-3, 
                   convType = 1, start = NULL, order = FALSE, 
                   alwaysUpdateWeights = FALSE, penType = 1){
  
  # Create Intermediate Variables
  ##############################
  intOnly = FALSE
  if(inherits(X, 'numeric') & all(X == 1)){
    X = as.matrix(rep(1, length(y)))
    intOnly = TRUE
    M = 1
  } else {
    M = ncol(X) + 1
  }
  N = length(y)
  Y0 <- which(y <= 0)
  Y1 <- which(y > 0)
  if(.Platform$OS.type == "unix"){
    sinkNul = "/dev/null"
  } else {
    sinkNul = 'NUL'
  }
  # Check inputs
  ##############################
  
  if(!is.numeric(y)){
    stop("y must be a numeric vector\n")
  }
  
  if(!is.matrix(X)){
    stop("X must be a matrix\n")
  }
  
  if(N <= 2){
    warning("sample size is too small\n")
    return(NULL)
  }
  
  if(nrow(X) != N){
    stop("the dimension of X and y do not match\n")
  }
  
  if(!is.null(unpenalizedx)){
    if((! is.numeric(unpenalizedx)) || nrow(unpenalizedx) != N){
      stop("unpenalizedx must be a numeric matrix of 
           the number of rows as the length as y\n")
    }
    useOffsetx = 1
  }else{
    useOffsetx = 0
    unpenalizedx = NULL
  }
  
  if(!is.null(unpenalizedz)){
    if((! is.numeric(unpenalizedz)) || nrow(unpenalizedz) != N){
      stop("unpenalizedz must be a numeric matrix of the number of rows 
           as the length as y\n")
    }
    useOffsetz = 1
  }else{
    useOffsetz = 0
    unpenalizedz = NULL
  }
  
  # Filter missing values
  ##############################
  
  isNA = apply(X, 1, function(v){ any(is.na(v)) })
  isNA = is.na(y) | isNA
  
  if(length(which(isNA))/N > naPercent){
    cat(y)
    print(X)
    stop("percent of missing data is too high\n")
  }
  
  w2kp  = which(!isNA)
  
  if(length(w2kp) == 0){
    warning("there is no non-missing values\n")
    return(NULL)
  }
  
  yR    = y[w2kp]
  XR    = X[w2kp, ,drop=FALSE]
  
  if(useOffsetx){
    unpenalizedx = unpenalizedx[w2kp, ,drop=FALSE]
  }
  
  if(useOffsetz){
    unpenalizedz = unpenalizedz[w2kp, ,drop=FALSE]
  }
  
  XR = scale(XR, center = FALSE) 
  yN    = yR
  
  
  
  if(order == TRUE){
    corr = cor(y, X)
    o = order(corr,decreasing=T)
    o.back = order(o)
    X = X[,o]
  } else {
    o.back = seq(1:ncol(X))
  }
  
  
  # Add intercept
  if(intOnly == FALSE){
    X = cbind(1, scale(X, center = FALSE))
  }
  
  # Initialize starting parameters and result 
  # container
  
  if(is.null(start)){
    start = list(betas = rep(0, M), gammas = rep(0, M))
    
    if(useOffsetx){
      fit_nb <- glm.nb(y~1+unpenalizedx)
      b0 <- coefficients(fit_nb)
      offsetx <- c(unpenalizedx %*% b0[-1])
    }else{
      offsetx <- rep(0, N)
    }
    
   if(useOffsetz){
     fit_log <- glm((y==0)~1+unpenalizedz, family = "binomial")
     g0 <- coefficients(fit_log)
     offsetz <- c(unpenalizedz %*% g0[-1])
   }else{
     offsetz <- rep(0, N)
   }
    

    
  } else if(start == 'jumpstart'){
    if(intOnly){
      start = list()
      
      if(useOffsetx){
        fit_nb <- glm.nb(y~1+unpenalizedx)
        b0 <- coefficients(fit_nb)
        b <- b0[1]
        offsetx <- c(unpenalizedx %*% b0[-1])
        start$betas = b
      }else{
        fit_nb <- glm.nb(y~1)
        b0 <- coefficients(fit_nb)
        b <- b0[1]
        offsetx <- rep(0, N)
        start$betas = b
      }
      
      if(useOffsetz){
        fit_log <- glm((y==0)~1+unpenalizedz, family = "binomial")
        g0 <- coefficients(fit_log)
        g <- g0[1]
        offsetz <- c(unpenalizedz %*% g0[-1])
        start$gammas = g 
      }else{
        fit_log <- glm((y==0)~1, family = "binomial")
        g0 <- coefficients(fit_log)
        g <- g0[1]
        offsetz <- rep(0, N)
        start$gammas = g 
      }
      
    } else {
      start = list()
      
      if(useOffsetx){
        fit_nb2 <- penalized(y, X[, -1, drop = FALSE], 
                             unpenalized = ~ 1 + unpenalizedx, 
                             model="poisson", lambda1 = 1, 
                             trace = FALSE)
        b <- c(coef(fit_nb2, "unpenalized")[1], coef(fit_nb2, "penalized"))
        offsetx <- c(unpenalizedx %*% coef(fit_nb2, "unpenalized")[-1])
        start$betas = b
      }else{
        fit_nb2 <- penalized(y, X[, -1, drop = FALSE], 
                             model="poisson", lambda1 = 1, 
                             trace = FALSE)
        b <- coef(fit_nb2, "all")
        offsetx <- rep(0, N)
        start$betas = b
      }
      
      if(useOffsetz){
        fit_log2 <- penalized(y==0, X[, -1, drop = FALSE], 
                              unpenalized = ~ 1 + unpenalizedz, 
                              model="logistic", lambda1 = 1, 
                              trace = FALSE)
        
        g <- c(coef(fit_log2, "unpenalized")[1], coef(fit_log2, "penalized"))
        offsetz <- c(unpenalizedz %*% coef(fit_log2, "unpenalized")[-1])
        start$gammas = g     
      }else{
        fit_log2 <- penalized(y==0, X[, -1, drop = FALSE], 
                              model="logistic", lambda1 = 1, 
                              trace = FALSE)
        
        g <- coef(fit_log2, "all")
        offsetz <- rep(0, N)
        start$gammas = g     
      }
      
      
      
    }
  } 
  
  betas = start$betas
  gammas = start$gammas
  
  theta = theta.st
  if(is.null(theta)){
    linkobj.z = binomial()
    linkobj.count = negative.binomial(1)
    
    etaz = X%*%gammas + offsetz
    pij = linkobj.z$linkinv(etaz)
    
    eta = X%*%betas + offsetx
    mu = linkobj.count$linkinv(eta)
    
    fnb = dnbinom(y, mu = mu, size = 1)
    
    ############################ Calculate expected value of latent var z
    zGy= (pij/(pij + fnb*(1-pij)))*(y==0)
    temp = updateTheta(betas, gammas, offsetx, offsetz, lambda, tau, theta.st, 
                       1, Y0, X, zGy, y)
    theta = temp[1]
  }
  
  
  nconv = 0
  it = 0
  
  result = list()
  
  if(weightedPen == TRUE){
    
    if(numericalDeriv == FALSE){
      
      betaWeight = hessBeta(X, rep(0,M), offsetz, rep(0,M), offsetx, y, 1)
      gammaWeight = hessGamma(X, rep(0,M), offsetz, rep(0,M), offsetx, y, 1)
      
    } else {
      betaWeight = -diag(hessian(updateObsLik.beta, rep(0,M),X= X, 
                                 gammas=rep(0,M), offsetz=offsetz, 
                                 offsetx=offsetx, y =y, Y0=Y0, 
                                 theta=1, lambda=0, tau =0.1))
      gammaWeight = -diag(hessian(updateObsLik.gamma, rep(0,M),X= X, 
                                  betas=rep(0,M), offsetz=offsetz, 
                                  offsetx=offsetx, y =y, Y0=Y0, 
                                  theta=1, lambda=0, tau =0.1))
    }
    ##############################################
    
    if(any(is.na(c(betaWeight, gammaWeight))) | 
       any(c(betaWeight, gammaWeight) < 0)){
      # browser()
      # warning(paste("Numerical problem calculating Hessian at:", 
      #               lambda, ",", tau,"// theta iteration:", it.theta))
      betaWeight = rep(1, length(betaWeight))
      gammaWeight = rep(1, length(gammaWeight))
    }
  } else {
    betaWeight = rep(1, M)
    gammaWeight = rep(1,M)
  }
  
  # Importan: Change Weights -------------
  betaWeight <- sqrt(abs(betaWeight))
  gammaWeight <- sqrt(abs(gammaWeight))
  
  tmp <- max(betaWeight)
  
  betaWeight <- betaWeight / tmp
  gammaWeight <- gammaWeight / tmp
  
  # Create tuning parameter grid if not given
  ##############################
  
  if(!is.null(lambdas)){
    nlambda = length(lambdas)
  }
  
  if(!is.null(taus)){
    ntau = length(lambdas)
  }
  
  if(is.null(lambdas) | is.null(taus)){
    tps = gentp(nlambda,ntau,y,X,unpenalizedx, unpenalizedz, offsetx,
                offsetz)
    
    lambdas = tps[,1]
    taus = tps[,2]    
  }
  
  
  ############################## Start Algorithm
  ##############################
  for(tp.iter in 1:length(lambdas)){
    
    warn = 0
    
    # Track only specific tuning parameter pair
    if(!is.null(track)){
      if(tp.iter != track){
        next 
      }
    }
    
    lambda = lambdas[tp.iter]
    tau = taus[tp.iter]
    # print(paste(Sys.time(),"TP:", tp.iter, " // lambda: ", lambda, 
    #             " // tau: ", tau))
    
    # no warm start if last iteration did not converge
    if(warmStart == FALSE){
      betas = start$betas
      gammas = start$gammas
      
      loglik = loglik.em = loglik.obs = loglik0.em = loglik1.em =
        -.Machine$double.xmax
      pen = .Machine$double.xmax
    } else if((it == maxIT | nconv > 0) & warmStart == 'cond'){
      betas = start$betas
      gammas = start$gammas
      
      loglik = loglik.em = loglik.obs = loglik0.em = loglik1.em =
        -.Machine$double.xmax
      pen = .Machine$double.xmax
    }
    
    ############################## Initialize
    
    theta.diff = .Machine$double.xmax
    it.theta = 0
    theta.record = vector()
    
    incrLikCount = 0
    nconv = 0
    
    thetamethod = 0
    
    result.theta = list()  
    
    
    ############################## We perform the EM algorithm within each 
    ############################## estimate of theta
    while(theta.diff > 1e-2 & it.theta < 1){
      ############################## Housekeeping
      theta.s = theta 
      it.theta = it.theta + 1
      # print(paste0("Theta: ", theta))
      if(oneTheta == TRUE){
        it.theta = 6
      }
      
      ############################## Initialize
      it = 0
      
      loglik = loglik.em = loglik.obs = loglik0.em = loglik1.em = BIC = 
        extBIC = extBICGG = -.Machine$double.xmax
      pen = .Machine$double.xmax
      diff.like = thresh = theta.diff = .Machine$double.xmax
      incrLikCount = 0
      nconv = 0
      thetamethod = 0
      
      result.EM = list()
      
      # no warm start if last iteration did not converge or across thetas
      if(warmStart == FALSE){
        betas = start$betas
        gammas = start$gammas
        
        loglik = loglik.em = loglik.obs = loglik0.em = loglik1.em =
          -.Machine$double.xmax
        pen = .Machine$double.xmax
      } else if((it == maxIT | nconv > 0) & warmStart == 'cond'){
        betas = start$betas
        gammas = start$gammas
        
        loglik = loglik.em = loglik.obs = loglik0.em = loglik1.em =
          -.Machine$double.xmax
        pen = .Machine$double.xmax
      }
      
      if(alwaysUpdateWeights == TRUE){
        if(weightedPen == TRUE){
          
          if(numericalDeriv == FALSE){
            betaWeight = hessBeta(X, rep(0,M), offsetz, rep(0,M), offsetx, y, 
                                  theta)
            gammaWeight = hessGamma(X, rep(0,M), offsetz, rep(0,M), offsetx, y, 
                                    theta)
          } else {
            betaWeight = -diag(hessian(updateObsLik.beta, rep(0,M),X= X, 
                                       gammas=rep(0,M), offsetz=offsetz, 
                                       offsetx=offsetx, y =y, Y0=Y0, 
                                       theta=theta, lambda=0, tau =0.1))
            gammaWeight = -diag(hessian(updateObsLik.gamma, rep(0,M),X= X, 
                                        betas=rep(0,M), offsetz=offsetz, 
                                        offsetx=offsetx, y =y, Y0=Y0, 
                                        theta=theta, lambda=0, tau =0.1))
          }
          ##############################################
          
          if(any(is.na(c(betaWeight, gammaWeight))) | 
             any(c(betaWeight, gammaWeight) < 0)){
            # browser()
            # warning(paste("Numerical problem calculating Hessian at:", lambda, 
            #               ",", tau,"// theta iteration:", it.theta))
            betaWeight = rep(1, length(betaWeight))
            gammaWeight = rep(1, length(gammaWeight))
          }
          
        } else {
          betaWeight = rep(1, M)
          gammaWeight = rep(1,M)
        }
      }
      
      
      if(optimType == "EM"){
        ############################## EM Algorithm loop
        ############################## 
        if(is.null(loud)){
          loudEM = 0
        } else {
          loudEM = loud
        } 
        if(convType == 1){
          res = updateEM(y, betas, offsetx, gammas, offsetz, X, theta, 
                              lambda, 
                              tau, irlsConv, betaWeight, gammaWeight, maxIT, 
                              loudEM, eps, 
                              model = 1, penType, maxIT2)
        } else {
          res = updateEM2(y, betas, offsetx, gammas, offsetz, X, theta, lambda, 
                          tau, irlsConv, betaWeight, gammaWeight, maxIT, loudEM,
                          eps, penType)
        }
        if(res$nconv == 1 | any(is.na(c(res$betas, res$gammas)))){
          fname = tempfile(pattern = "errorFile_", tmpdir = errorDirec, 
                           fileext = ".RData")
          save(y, betas, offsetx, gammas, offsetz, X, theta, lambda, tau, 
               irlsConv, betaWeight, gammaWeight, maxIT, loudEM, file=fname)
          # browser()
          nconv = 1
          next
        }
        
        betas = res$betas
        gammas = res$gammas
        
        temp = updateTheta(betas, gammas, offsetx, offsetz, lambda, tau, 
                           theta.st, theta, Y0, X, zGy, y)
        
        theta.record[it.theta] = theta = temp[1]
        thetamethod = temp[2]
      } 
      
      if(nconv == 1){
        next
      }
      
      # Start BFGS Optimization Coordinate Descent Loop
      it = 0
      while(abs(diff.like) > 1e-5 & it < maxOptimIT & incrLikCount < 5){
        
        ############################# Step Through
        if(!is.null(stepThrough)){
          if(stepThrough[1] == it.theta & stepThrough[2] == it){
            browser()
          }
        }
        
        # Housekeeping
        it = it + 1
        loglik.s = loglik
        
        sto = list(betas = betas, gammas = gammas, loglik.obs = loglik.obs, 
                   pen = pen, theta.r = theta.record, theta = theta,
                   BIC = BIC, extBIC = extBIC, extBICGG = extBICGG)
        
        # Calculate current mean estimates
        linkobj.z = binomial()
        linkobj.count = negative.binomial(theta = theta)
        
        etaz = X%*%gammas + offsetz
        pij = linkobj.z$linkinv(etaz)
        
        eta = X%*%betas + offsetx
        mu = linkobj.count$linkinv(eta)
        
        fnb = dnbinom(y, mu = mu, size = theta)
        
        # Maximizing over each gene
        useBFGS = 1
        for(k in 1:ncol(X)){
          res = updateOptim.k(y, betas, gammas, theta, k, offsetx, offsetz, 
                              lambda, tau, X, theta.st, betaWeight, gammaWeight,
                              useBFGS)
          nconv = nconv + res[3]
          if(nconv > 0){
            break
          }
          
          betas[k] = res[1]
          gammas[k] = res[2]
          useBFGS = res[4]
        }
        
        if(nconv > 0){
          break
        }
        
        ## Update theta
        if(oneTheta == TRUE){
          theta = optimize(optimThetaFunc, c(0.001,100), X = X, gammas = gammas, 
                           offsetz = offsetz, betas = betas, offsetx = offsetx, 
                           y = y, 
                           Y0 = Y0, lambda = lambda, tau = tau, maximum = TRUE, 
                           betaWeight = betaWeight, 
                           gammaWeight = gammaWeight)$maximum
        }
        
        # Calculate penalty
        if(penType == 1){
          pen = sum(lambda*log(abs(betaWeight[-1]*betas[-1]) +
                                 abs(gammaWeight[-1]*gammas[-1])+tau))
        } else if (penType == 2){
          pen = sum(lambda*(abs(betaWeight[-1]*betas[-1]) +
                              abs(gammaWeight[-1]*gammas[-1])+tau))
        }
        
        # Get Likelihoods
        loglik.obs = updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, 
                                  theta, lambda, tau)
        
        loglik = loglik.obs - pen
        
        diff.like = loglik - loglik.s
        
        # ############# Record EM path if tracking
        
        # if(!is.null(track)){
        #   result.EM[[it]] = list(loglik.obs = loglik.obs, pen = pen, 
        #                          betas = betas, gammas = gammas)
        # }
        
        ############# Calculate End Criteria
        if(is.na(diff.like)){
          break
        }
        if(is.null(track)){
          if(diff.like < -0.01) {
            lowestLikSto = sto
            # incrLikCount = incrLikCount + 1
          } else {
            incrLikCount = 0
            lowestLikSto = NULL
          }
        }
        
        if(diff.like < -0.01){
          # print(paste(">>>>>>>>>>ERROR<<<<<<<<<<<<<<<<", it, ":", diff.like, 
          # " // theta method:", thetamethod, "// theta:", theta))
          warn = 1
        } 
        
        if(!is.null(loud)){
          if(it %% loud == 0){
            print(paste("BFGS:", it, ":", diff.like, " // theta method:", 
                        thetamethod, "// theta:", theta))
          } 
        }
      } # end BFGS coord descent while loop
      
      
      if(optimType == "optim"){
        theta = optimize(optimThetaFunc, c(0.001,100), X = X, gammas = gammas, 
                         offsetz = offsetz, betas = betas, offsetx = offsetx, 
                         y = y, 
                         Y0 = Y0, lambda = lambda, tau = tau, maximum = TRUE, 
                         betaWeight = betaWeight, 
                         gammaWeight = gammaWeight)$maximum
      }
      
      ############# Calculate penalty
      if(penType == 1){
        pen = sum(lambda*log(abs(betaWeight[-1]*betas[-1])+
                               abs(gammaWeight[-1]*gammas[-1])+tau))
      } else if (penType == 2){
        pen = sum(lambda*(abs(betaWeight[-1]*betas[-1])+
                            abs(gammaWeight[-1]*gammas[-1])+tau))
      }
      
      ############# Get Likelihoods
      loglik.obs = updateObsLik(X, gammas, offsetz, betas, offsetx, y, Y0, 
                                theta, lambda, tau)
      
      loglik = loglik.obs - pen
      
      theta.diff = abs(theta - theta.s)
      
      ############ Record theta path if tracking
      if(!is.null(track)){
        result.theta[[it.theta]] = list(loglik.obs = loglik.obs, pen = pen, 
                                        betas = betas, gammas = gammas, 
                                        theta = theta)
      }
      
      # print(paste(it.theta, "// theta diff:", theta.diff, 
      #             " // theta:", theta))
    } # end theta while loop
    
    ############# Calculate BIC
    if(is.null(bicgamma)){
      bicgamma = log(M)/log(N)
    }
    BIC = -2.0*loglik.obs + log(N)*(sum(betas != 0) + sum(gammas != 0) - 2)
    extBIC = BIC + 2*log(M)*bicgamma*(sum(betas != 0) + sum(gammas != 0) - 2)
    extBICGG = BIC + bicgamma*2*log(choose((ncol(X)-1)*2,
                                           (sum(betas != 0) + 
                                            sum(gammas != 0) - 2)))
    
    ############ Return values for tp pair
    if(!is.null(track)){
      if(length(track) > 1 | track == "all"){
        ############ Record tp path if tracking
        if(!is.null(track)){
          result[[paste(tp.iter)]] = append(result.theta, 
                                            list(X= X, lambda = lambda, 
                                                 tau = tau))
        }
      } else {
        return(append(result.theta, list(X = X, lambda = lambda, tau = tau)))
      }
    } else if (incrLikCount > 4){
      res.tp = lowestLikSto
      result[[paste(tp.iter)]] = res.tp
    } else {
      if(tp.iter == 1){
        if(intOnly == FALSE){
          X.output = cbind(X[,1], X[,-1, drop=FALSE][,o.back])          
        } else {
          X.output = X
        }
      } else {
        X.output = NULL
      }
      if(intOnly == FALSE){
        betas.output = c(betas[1], betas[-1][o.back])
        gammas.output = c(gammas[1], gammas[-1][o.back])        
      } else {
        betas.output = betas
        gammas.output = gammas
      }
      res.tp = list(X = X.output, 
                    betas = betas.output[which(betas.output != 0)], 
                    gammas = gammas.output[which(gammas.output!=0)], 
                    betas.w = which(betas.output!=0), 
                    gammas.w = which(gammas.output!=0), 
                    loglik.obs = loglik.obs, pen = pen, theta.r = theta.record, 
                    theta = theta, BIC = BIC, extBIC = extBIC, 
                    extBICGG = extBICGG, 
                    lambda = lambda, tau = tau)
      result[[paste(tp.iter)]] = res.tp
    } 
    if(warn == 1){
      warning(paste("Decrease in EM likelihood for tuning parameter pair:", 
                    lambda,",",tau))
    }
  } # end tp for loop
  return(result)
}

