#Test Distributions
library(pracma)
library(gsl)
rm(list=ls())

##################Notes for SimNet Version 2.1: Optimized Blitzstein & Diaconis's Sequential Method
# 1. Include Erdős–Gallai at the beginning and each sampling
# 2. With in-loop EG checking, the generation is guaranteed to be successful
# 3. Each edge will have an EG checking before sampling from candidate, so with no restarting procedure
# 4. Applied Optimization ( Mark 2 & Mark 3(Durfee Number)) by their article but slightly harder to illustrate

### Note: Trying to add Moseman's Improvement later
### The Durfee Number optimization is performs better than original for EG checking thus the EG checking function should be this version
# 5. Have not introduce the importance sampling technique to guarantee uniform!!!!!!!!!!!!!
# 6. Asymptotically Uniform? (by Greenhill ariticle? need to do more research) But assure that all graph will have positive prob to appear and close to uniform

# 5. Pros: Fastest method for large network for long tailed distribution like powerlaw and exponential, and with higher mean degree. Guarantee no failure and no restart
# 6. Cons: Slow for distribution with low mean, distribution with higher degree concentration and small network size
# 7. Checked Optimization of Ver 2.0
# 8. good for long-tail distribution, high mean and large size? (need to test in future)
# 9. Try Multi-Core method?

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
    
    # Now apply Mark2(Based on Test file1): for each VIndex, candidates won't change
    
    # Step 4: Compute candidate list
    # Condition 1: i \neq j
    # Initialized NIndex without previous CIndex (Candidate list with same VIndex)
      NIndex <- ResIndex[-which(ResIndex==VIndex)]
    
    # Step 7: (Recursive) repeat step 4-6 until degree of VIndex is 0
    while (VDegree>0) {
      # Step 4: Compute candidate list
      # Condition 1: i \neq j (Did previously)
      
      # Condition 2: {i,j} \notin E
      VPrevEdge <- which(Gmatrix[VIndex, ]==1)
      NIndex <- NIndex[! NIndex %in% VPrevEdge]
      
      # Condition 3: after add edges, the residual seq still realizable
      
      # EG_condition initialize
      L <- length(NIndex)
      EG_ResSeq <- rep(0,L)
      
      # EG check for each possible index satisfy C1 and C2
      # No need to check even property: by add stubs to both vertices connected to an edge, we assure the residual total degree is always even
      for (i in c(1:L)) {
        JIndex <- NIndex[i]
        TestDegree <- ResDegree
        TestDegree[VIndex] <- TestDegree[VIndex]-1
        TestDegree[JIndex] <- TestDegree[JIndex]-1
        EG_ResSeq[i] <- EG_check(TestDegree)
      }
      
      # Candidate list
      CIndex <- NIndex[which(EG_ResSeq==1)]
      
      # Step 5: Pick Nsample with probability proportional to its current residual degree
      # Current residual degree of candidates 
      CDegree <- ResDegree[CIndex]
      
      # Degree distribution of candidates
      CTotalDegree <- sum(CDegree)
      CDist <- CDegree/CTotalDegree
      
      # Sampling
      if (length(CIndex)==1){
        Nsample <- CIndex
      } else {
        Nsample <- sample(CIndex,1,replace = FALSE, prob = CDist)  
      }
      
      # Step 6: add the edge {VIndex,Nsample} to E and update ResDegree
      Gmatrix[VIndex,Nsample] <- 1
      Gmatrix[Nsample,VIndex] <- 1
      
      ResDegree <- DegreeSample-rowSums(Gmatrix)
      VDegree <- ResDegree[VIndex]
      
      # Now apply Mark2(Based on Test file1): for each VIndex, candidates won't change in same stage
      NIndex <- CIndex
    }
    Prog <- round(1-length(which(ResDegree>0))/N,4)*100
    print(Prog)
  }
  return(list(DegreeDist=DegreeSample,EdgeMatrix=Gmatrix))
}

