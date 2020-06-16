#' Estimate the skeleton of the graphical model using a modified PC algorithm
#' 
#'
#' After estimating moral graph using neighborhood selection, the false 
#' positive connections can be removed using a modified PC algorithm, which 
#' removes an edge 
#' between two variables if they are conditional independent given any subset 
#' of other variables.
#' 
#' @param estG the estimated moral graph model
#' @param variables the sample data
#' @param alpha the level of significance used in the conditional independence 
#' test - default is 0.01
#' @param fun the function used for the conditional independence test
#' @export
#' @return A list consisting of two components:
#' \describe{
#'   \item{res}{the removed vertex pairs and their d-separation sets}
#'   \item{X}{simulated sample data}
#'   \item{G}{the matrix adjacency matrix for the resulting graph}
#' } 
pcZINB <- function(estG, variables, alpha = 0.01, fun = cond_ZINB){
  #Implement Stable PC Algorithm
  if(.Platform$OS.type == "unix"){
    sinkNul = "/dev/null"
  } else {
    sinkNul = 'NUL'
  }
  
  l = -1
  p = ncol(estG)
  sepset = pval.stor = pl = vector("list", p)
  for (i in 1:p) sepset[[i]] = pl
  G = estG
  if(sum(G) == 0){
    warning("Nothing was selected in Step 1: Your data may not satisfy assumptions of this method.")
    
    res = list()
    output2 = list()
    output2[["first"]] = NA
    output2[["second"]] = NA
    output2[["third"]] = "Nothing was selected in Step 1."
    
    res[[1]] = output2
    
    result = list()
    result[["res"]] = res
    result[["G"]] = G
    return(result)
  }
  S = list()
  storage = NULL
  
  connect <- as.matrix(which(G == 1, arr.ind=T))
  
  connect = connect[which(connect[,1] < connect[,2]),,drop=FALSE]
  
  # ? test marginal dependence 
  if(nrow(connect) > 0){
    for(iterd in 1:nrow(connect)){
      i = connect[iterd,1]
      j = connect[iterd,2]
      
      if(estG[i,j] == 1 & (j > i)){
        #Test Marginal Dependence:
        x = variables[,i]
        y = variables[,j]
        
        # sink(sinkNul) 
        result.temp = fun(x,y,numeric(0), alpha)
        # sink()
        if(result.temp > alpha){
          G[i,j] = 0
          G[j,i] = 0
          storage = rbind(storage,c(i,j,0,result.temp,NA,NA))
        }
      } #end estG
    }
  }
  
  
  disCG = matrix(10, nrow = nrow(G), ncol=ncol(G))
  
  # connect <- matrix(1:2, nrow = 1)
  # connect2 <- matrix(1:4, nrow = 2)
  # if conditional set with size 6
  
  # while(max(disCG, na.rm= T) > l & nrow(connect) != nrow(connect2)){
  while(max(disCG, na.rm= T) > l){
    disCG = matrix(NA, nrow = nrow(G), ncol=ncol(G))
    
    l = l+1 # Filtering for conditional sizes: at most 10
    print(paste("Filtering for conditional sizes:", l))
    G.tilde = G
    
    #Build list of connections:
    connect <- as.matrix(which(G.tilde == 1, arr.ind=T))
    connect = connect[which(connect[,1] < connect[,2]),,drop=FALSE]
    
    if(nrow(connect) > 0){
      res = foreach(iterd = 1:nrow(connect), 
                    .export=c("connectedComp", "cond_ZINB", 
                              "penZINB")) %dopar% { 
                                tryCatch({
                                  if(iterd %in% floor(quantile(c(1:nrow(connect)), 
                                                               c(0.25, 0.5, 0.75)))){
                                    print(paste("Removal:", iterd,"out of",nrow(connect)))
                                  }
                                  
                                  output = NA
                                  
                                  Sij = vector()
                                  
                                  i = connect[iterd, 1]
                                  j = connect[iterd, 2]
                                  
                                  nbrsUnion = (G.tilde[, i] | G.tilde[, j])
                                  nbrsUnion[i] = FALSE
                                  nbrsUnion[j] = FALSE
                                  nbrsIsect = (G.tilde[, i] & G.tilde[, j])
                                  nbrsIsect[i] = FALSE
                                  nbrsIsect[j] = FALSE
                                  
                                  nbrsdec = nbrsIsect
                                  if (sum(nbrsIsect) != 0) {
                                    Gsub = G.tilde # stability
                                    Gsub[,c(i,j)] = Gsub[c(j,i),] = FALSE
                                    nbrsdec = (connectedComp(Gsub,
                                                             which(nbrsIsect)) & nbrsUnion)
                                  }
                                  
                                  Aij = which(nbrsUnion)
                                  Cij = which(nbrsdec)
                                  
                                  #step = 1
                                  
                                  ##Begin 3.2
                                  if(length(Cij) >= l){
                                    #Iteratively select a subset in Cij
                                    sw = 0 #switch for whether we removed the edge
                                    
                                    if(length(Cij) == l){
                                      Gamma = matrix(Cij,ncol = 1)
                                    } else {
                                      Gamma = combn(Cij, l) #all possible sets
                                    }
                                    
                                    Gamma.i = 0
                                    while(sw == 0){ 
                                      Gamma.i = Gamma.i + 1
                                      Gamma.temp = Gamma[,Gamma.i]
                                      
                                      if(length(which(Aij %in% Gamma.temp)) > 0){
                                        Kappa = Aij[-which(Aij %in% Gamma.temp)]  
                                      } else {
                                        Kappa = Aij
                                      }
                                      
                                      # step = 2
                                      x = variables[,i]
                                      y = variables[,j]
                                      
                                      if(length(Kappa) == 0){
                                        sw = 1
                                        next
                                      } else {
                                        z = variables[,Kappa]
                                      }
                                      
                                      # sink(sinkNul) 
                                      result.temp = fun(x,y,z,alpha)
                                      # sink()
                                      
                                      if(result.temp > alpha){
                                        G[i,j] = 0
                                        G[j,i] = 0
                                        Sij = Kappa
                                        sw =  1
                                        storage = rbind(storage,c(i,j,0, result.temp,
                                                                  NA,l))
                                      }
                                      
                                      if(Gamma.i == ncol(Gamma)){
                                        sw = 1
                                      }
                                    }
                                  } #End 3.2
                                  
                                  output2 = list()
                                  output2[["first"]] = c(i, j, G[i,j])
                                  output2[["second"]] = Sij
                                  output2[["third"]] = storage
                                  
                                  # res[[iterd]] = output2
                                  return(output2)
                                }, error = function(e) {
                                  output2 = list()
                                  output2[["first"]] = NA #MAYBE THIS NEEDS TO BE VECTOR OF 3?
                                  output2[["second"]] = NA
                                  output2[["third"]] = NA
                                  output2[["error"]] = paste("The", 
                                                             iterd,
                                                             "th Iteration (with i=",i," and j=",
                                                             j,") Caused Error:",e,".")
                                  print("Error in pcZINB, check output2.")
                                  return(output2)
                                }) #end tryCatch
                              } #end parallel
      
    }
    
    #CHECK RES HERE
    
    res2 = t(sapply(res, "[[", "first"))
    G = convert(res2, G)
    
    #REMEMBER STILL NEED TO FIGURE OUT A WAY TO STORE SIJ
    
    disCG = vector()
    for(i in 1:(ncol(G) - 1)){
      for(j in (i+1):ncol(G)){
        if(G[i,j] != 0){
          nbrsUnion = (G.tilde[, i] | G.tilde[, j])
          nbrsUnion[i] = FALSE
          nbrsUnion[j] = FALSE
          nbrsIsect = (G.tilde[, i] & G.tilde[, j])
          nbrsIsect[i] = FALSE
          nbrsIsect[j] = FALSE
          
          nbrsdec = nbrsIsect
          if (sum(nbrsIsect) != 0) {
            Gsub = G.tilde # stability
            Gsub[,c(i,j)] = Gsub[c(j,i),] = FALSE
            nbrsdec = (connectedComp(Gsub,which(nbrsIsect)) & nbrsUnion)
          }
          
          Cij = which(nbrsdec)
          
          disCG = c(disCG, length(Cij))
        }
      }
    }
  }
  
  result = list()
  result[["res"]] = res
  result[["G"]] = G
  return(result)
}  