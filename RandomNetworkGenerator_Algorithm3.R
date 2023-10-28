#Test Distributions
library(pracma)
library(gsl)
rm(list=ls())



##################Notes for SimNet Version 3.0: Mixed Sequential Method
# 1. Include Erdős–Gallai at the beginning and each sampling
# 2. With in-loop EG checking, the generation is guaranteed to be successful
# 3. All residual stubs will be sampled each time and then check EG 
# 4. This will leads to large restarting for higher degree vertex
# 5. No analysis about the distribution

# 6. Pros: faster for smaller network with concentrated distribution like Mean Field ans Poisson, and with higher mean degree. Guarantee no failure
# 7. Cons: No rigorous analysis now, restarting problem.

### Maybe should do a time comparison test for different algorithm


# Graph realization problem #
# Condition 1: total degree must be even

#Adjust to even total degree by Newman's Procedure

# sum(EX_sample)%%2
# 
# while (sum(EX_sample)%%2 == 1) {
#   alt_index <- sample(N,1)
#   EX_sample[alt_index] <- sample(c(0:max(Pk_EX$kvalue)),1,replace = T,prob = Pk_EX$pk)
# }
# sum(EX_sample)

# Condition 2: The degree distribution can actually form a network
# Havel–Hakimi algorithm is actually performed when construct the network below: low calculation efficiency
# Erdős–Gallai theorem: Equivalent but seems more simple

EG_check <- function(DegreeDist){
  check_vec <- c(0)
  DegreeDist <- sort(DegreeDist,decreasing = T)
  N <- length(DegreeDist)
  Cum <- cumsum(DegreeDist)
  
  #Mark3: corrected Durfee number m 
  Dlist <- DegreeDist- c(0:(N-1)) >=0
  m <- length(which(Dlist==TRUE))
  
  for (k in c(1:m)) {
    RHS <- k*(k-1) + k*(length(DegreeDist[which(which(DegreeDist>=k)>k)]))+sum(DegreeDist[which(DegreeDist<k)])
    LHS <- Cum[k]    
    if (LHS<=RHS) {
      check_vec[k] <- 1
    } else {
      check_vec[k] <- 0
    }
  }
  return(min(check_vec))
}

#"A sequential importance sampling algorithm for generating random graphs with prescribed degrees" by JOE. Blitzstein and P. Diaconis
# Method update to assure generalization

# Maybe later
# make the edge sparse matrix to adjacency list to solve the RAM/Storage bottle neck problem for larger networks


# New Sequential Algorithm
# Network generate function
SimNet <- function(Pk, N=500){
  EG_result <- 0
  
  while (EG_result==0) {
    DegreeSample <- sample(c(0:max(Pk$kvalue)),N,replace = T,prob = Pk$pk)
    DegreeSample <- sort(DegreeSample, decreasing = T)
    
    # Step 0: Graph realization problem, make sure a graphical/realizable sequence is generated
    
    # Condition 1: total degree must be even
    while (sum(DegreeSample)%%2 == 1) {
      alt_index <- sample(N,1)
      DegreeSample[alt_index] <- sample(c(0:max(Pk$kvalue)),1,replace = T,prob = Pk$pk)
    }
    
    # Condition 2: The degree distribution can actually form a network
    EG_result <- EG_check(DegreeSample)
  }
  DegreeSample <- sort(DegreeSample, decreasing = F)
  
  # Step 1: create empty list of edges
  Gmatrix <- matrix(0, nrow = N, ncol = N)
  ResDegree <- DegreeSample-rowSums(Gmatrix)
  
  # Step 8: (recursive) return to step 2
  # Step 2: if d=0, terminate with output
  while (max(ResDegree)>0) {
   
    # Step 3: Choose the least VIndex with VDegree a minimal positive entry
    
    # Condition 1: positive
    ResIndex <- which(ResDegree>0)

    # Condition 2: minimal
    # Always paring from lower leftover degree to avoid not paired vertices at the end
    VIndex <- ResIndex[which.min(ResDegree[ResIndex])]
    VDegree <- ResDegree[VIndex]
    
    # New Algorithm Test: randomly sampling #VDegree vertices based on distribution proportional to its residual degree
    
    # Step 4: Compute candidate list
    # Condition 1: i \neq j
    NIndex <- ResIndex[-which(ResIndex==VIndex)]
    
    # Condition 2: {i,j} \notin E
    VPrevEdge <- which(Gmatrix[VIndex, ]==1)
    NIndex <- NIndex[! NIndex %in% VPrevEdge]
    
    # Condition 3: after add edges, the residual seq still realizable
    EG_Test <- 0
    
    # Then loop until the sample pass the EG test
    while (EG_Test == 0) {
      # Create sampling distribution based on residual degree of possible candidates
      NDegree <- ResDegree[NIndex]
      
      # Sampling
      if (length(NIndex)==1){
        Nsample <- NIndex
      } else {
        Nsample <- sample(NIndex, VDegree, replace = FALSE, prob = NDegree)  
      }
      
      # EG check for sample
      # No need to check even property: by add stubs to both vertices connected to an edge, we assure the residual total degree is always even
    
      # Update test edges
      TestDegree <- ResDegree
      TestDegree[VIndex] <- TestDegree[VIndex]-VDegree
      TestDegree[Nsample] <- TestDegree[Nsample]-1
      
      EG_Test <- EG_check(TestDegree)
    }
      # Step 6: add the edge {VIndex,Nsample} to E and update ResDegree
      Gmatrix[VIndex,Nsample] <- 1
      Gmatrix[Nsample,VIndex] <- 1
      
      ResDegree <- DegreeSample-rowSums(Gmatrix)
      #Prog <- max(ResDegree)
      Prog <- round(1-length(which(ResDegree>0))/N,4)*100
      print(Prog)
  }
  return(list(DegreeDist=DegreeSample,EdgeMatrix=Gmatrix))
}
