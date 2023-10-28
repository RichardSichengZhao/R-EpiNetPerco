#Test Distributions
library(pracma)
library(gsl)
rm(list=ls())

GilAlgo <- function(Network, beta, gamma, MaxTime, InitInfSize=1, TrackDyn=F){
  # Since we sort the vertices with degree, we need random initial infection with size s
  s <- InitInfSize
  
  G <- Network
  N <- length(G[1,])
  g <- gamma
  b <- beta
  
  # Radomly chose s vertices to be infected
  InitIndex <- sample.int(N,s)
  
  #Status: S=0, I=1, R=2

  #Initialize status and rate
  t <- 0
  Status <- rep(0,N)
  Status[InitIndex] <- 1
  Istep <- length(which(Status==1))
  
  if (TrackDyn==T){
    NumStep <- 1
    t_vec <- c(t)
    S_vec <- c(length(which(Status==0))/N)
    I_vec <- c(length(which(Status==1))/N)
    R_vec <- c(0)
  }

  Rate <- rep(0,N)
  Rate[InitIndex] <- g

  for (i in c(1:s)) {
    Contact <- which(G[InitIndex[i],]>0 & Status<1)
    Rate[Contact] <- b
  }

  # while loop start t<tmax & Istep != 0
  while(t<MaxTime & Istep != 0){
    Sum <- sum(Rate)
    Cum <- cumsum(Rate)
    r <- runif(2, min = 0, max = 1)

    Event <- min(which(Cum>r[1]*Sum))
    Status[Event] <- Status[Event]+1
    Contact <- which(G[Event,]>0 & Status<1)

    if (Status[Event]==2){
      Rate[Event] <- 0
      Rate[Contact] <- Rate[Contact]-b
    } else if (Status[Event]==1){
      Rate[Event] <- g
      Rate[Contact] <- Rate[Contact]+b
    } else {
    }
  
    Tstep <- -log(r[2])/Sum  
    t <- t+Tstep
    Istep <- length(which(Status==1))
  
    if (TrackDyn==T){
      NumStep <- NumStep+1
      t_vec[NumStep] <- t
      S_vec[NumStep] <- length(which(Status==0))/N
      I_vec[NumStep] <- length(which(Status==1))/N
      R_vec[NumStep] <- length(which(Status==2))/N
    }
  }
  
  FinishTime <- t
  Ssize <- length(which(Status==0))/N
  Isize <- length(which(Status==1))/N
  Rsize <- length(which(Status==2))/N
  FinalStat <- data.frame(FinishTime,Ssize,Isize,Rsize)
  
  if (TrackDyn==T){
    Track <- cbind(t_vec,S_vec,I_vec,R_vec)
    return(list(FinalStat=FinalStat,Details=Track))
  } else {
    return(FinalStat) 
  }
}

