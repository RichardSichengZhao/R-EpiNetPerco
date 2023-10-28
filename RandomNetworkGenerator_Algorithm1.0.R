#Test Distributions
library(pracma)
library(gsl)
rm(list=ls())



##################Notes for SimNet Version 1.0: Old original Method
# 1. Include Erdős–Gallai but only at the beginning
# 2. Without in-loop EG checking, the generation can be failure due to sample not-graphical/not-realizable residual degree
# 3. Because of 2, there is a cap for number of loops in each generation procedure
# 4. Also, the sampling procedure is fastest since for each vertex, all edges will constructed in one simple sampling.
### Note this cap should be parameter, now it is a constant=15
# 5. Pros: Fastest method for large network generation for any distribution
# 6. Cons: Random failure and the graph not generated uniformly(No rigorous proof)
# 7. Good for quick test of large network.

# Graph realization problem #
# Condition 1: total degree must be even

# Condition 2: The degree distribution can actually form a network
# Havel–Hakimi algorithm is actually performed when construct the network below: low calculation efficiency
# Erdős–Gallai theorem: Equivalent but seems more simple
EG_check_V1 <- function(DegreeDist, N){
  check_vec <- c(0)
  Cum <- cumsum(DegreeDist)
  for (k in c(1:N)) {
    RHS <- k*(k-1) + k*(length(DegreeDist[which(which(DegreeDist>=k)>k)]))+sum(DegreeDist[which(DegreeDist<k)])
    LHS <- Cum[k]    
    if (LHS<=RHS){
      check_vec[k] <- 1
    } else{
      check_vec[k] <- 0
    }
  }
  return(min(check_vec))
}

# make the edge sparse matrix to adjacency list to solve the RAM/Storage bottle neck problem for larger networks

# Network generate function
SimNet_V1 <- function(Pk, N=500){
  EG_result <- 0
  
  while(EG_result==0){
    DegreeSample <- sample(c(0:max(Pk$kvalue)),N,replace = T,prob = Pk$pk)
    DegreeSample <- sort(DegreeSample, decreasing = T)
    
    # Graph realization problem
    
    # Condition 1: total degree must be even
    while (sum(DegreeSample)%%2 == 1) {
      alt_index <- sample(N,1)
      DegreeSample[alt_index] <- sample(c(0:max(Pk$kvalue)),1,replace = T,prob = Pk$pk)
    }
    
    # Condition 2: The degree distribution can actually form a network
    EG_result <- EG_check(DegreeSample,N)
  }
  
  DegreeSample <- sort(DegreeSample, decreasing = T)
  Gmatrix <- matrix(0, nrow = N, ncol = N)
  ResDegree <- DegreeSample-rowSums(Gmatrix)
  
  temp_ResDegree <- ResDegree
  failure <- 0
  
  while (max(ResDegree)>0) {
    ResIndex <- which(ResDegree>0)
    #Always paring from max leftover degree to avoid not paired vertices at the end
    VIndex <- which.max(ResDegree)
    NIndex <- ResIndex[-which(ResIndex==VIndex)]
    
    VPrevEdge <- which(Gmatrix[VIndex, ]==1)
    NIndex <- NIndex[! NIndex %in% VPrevEdge]
    VDegree <- ResDegree[VIndex]
    
    #Edges paring can be failure due to randomness
    #while length(NIndex) < VDegree, left over vertices cannot satisfy the degree
    #Restart the paring until success
    
    if (length(NIndex)<VDegree){
      ResDegree <- temp_ResDegree
      failure <- failure+1
      if (failure > 15){
        break
      }
    } else {
      if (length(NIndex)==1){
        Nsample <- NIndex
      } else {
        Nsample <- sample(NIndex,VDegree,replace = FALSE)  
      }
    Gmatrix[VIndex,Nsample] <- 1
    Gmatrix[Nsample,VIndex] <- 1
    ResDegree <- DegreeSample-rowSums(Gmatrix)
    }
  }
  return(list(DegreeDist=DegreeSample,EdgeMatrix=Gmatrix, fail=failure))
}

