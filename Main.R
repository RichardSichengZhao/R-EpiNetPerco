library(pracma)
library(gsl)
rm(list=ls())

set.seed(12345)

#Translate per-node degree distribution column vector into data frame
#For future network simulation usage
DegreeDF <- function(dist){
  n <- length(dist)
  degree <- sort(dist)
  node <- c(1:n)
  df <- data.frame(node, degree)
  return(df)
}
#dist is the raw data vector for degree without sorting
#df is the data frame with node number with ascending ordered degree
#This is for create degree distribution from raw data


#Degree Distribution data frame
DDistPK <- function(df){
  m <- max(df[,2])
  kvalue <- c(0:m)
  pk <- rep(0,m)
  for (k in kvalue) {
    pk[k+1] <- length(c(which(df[,2]==k)))/nrow(df)
  }
  Pk <- data.frame(kvalue,pk)
  return(Pk)
}
# Creat degree distribution frame from raw data df


#PGFs and Derivatives

#G0
PGFG0 <- function(x,Pk){
  G0 <- 0
  for (k in Pk[,1]) {
    G0 <- G0+(x^k)*(Pk[k+1,2])
  }
  return(G0)
}

#G'0
PGFd1G0 <- function(x,Pk){
  d1G0 <- 0
  for (k in c(1:max(Pk[,1]))) {
    d1G0 <- d1G0+(x^(k-1))*k*(Pk[k+1,2])
  }
  return(d1G0)
}

#G''0
PGFd2G0 <- function(x,Pk){
  d2G0 <- 0
  m <- max(Pk[,1])
  for (k in c(2:m)) {
    d2G0 <- d2G0+(x^(k-2))*k*(k-1)*(Pk[k+1,2])
  }
  return(d2G0)
}

#<K^n>
Kn <- function(Pk,n){
  Knvalue <- 0
  for (k in Pk[,1]) {
    Knvalue <- Knvalue+(k^n)*(Pk[k+1,2])
  }
  return(Knvalue)
}

#G1
PGFG1 <- function(x,Pk){
  G1 <- PGFd1G0(x,Pk)/Kn(Pk,1)
  return(G1)
}

#G'1
PGFd1G1 <- function(x,Pk){
  G1 <- PGFd2G0(x,Pk)/Kn(Pk,1)
  return(G1)
}


#u_T=G_q(u_T) self contain equation
ueqn <- function(x) {
  PGFG1(1+(x-1)*Tvalue,Pk_value)-x
}


##Changing input from beta, gamma to T
##Two different type of constant T assumption
#Newman's concentration assumption
Tconst_Newman <- function(beta, gamma){
  Tvalue <-1-exp(-beta/gamma)
  return(Tvalue)
}

#SIR exponential assumption
Tconst_exp <- function(beta,gamma){
  Tvalue <-beta/(beta+gamma)
  return(Tvalue)
}

# Typical Percolation Process
TypProc <- function(Pk,Tvalue,tol=1e-3){
  Tc_value <- 1/(PGFd1G1(1,Pk))
  OBType <- ''
  s <- 0
  Rinfty <- 0
  u <- 0
  v <- 0
  ueqn <- function(x) {
    PGFG1(1+(x-1)*Tvalue,Pk)-x
  }
  if (Tvalue<Tc_value){
    OBType <- 'Limited'
    s <- 1+Tvalue*PGFd1G0(1,Pk)/(1-Tvalue*PGFd1G1(1,Pk))
    Rinfty <- 0
    u <- 1
    v <- 1
  }
  else if(Tvalue==Tc_value){
    OBType <- 'Undefined'
    s <- 0
    Rinfty <- 0
  }
  else if(Tvalue>Tc_value){
    OBType <- 'Epidemic'
    s <- 0
    #usol <- uniroot(ueqn,c(0+tol,1-tol),tol = 1e-11)
    
    LB_u <- 0
    UB_u <- 1
    u_vec <- seq(from=LB_u,to=UB_u,by=tol)
    u_mat <- matrix(0,nrow = length(u_vec),ncol = 3)
    u_mat[,1] <- u_vec
    u0 <- ueqn(0) 
    
    for (i in c(1:length(u_vec))) {
      u_mat[i,2] <- ueqn(u_vec[i])
      u_mat[i,3] <- u_mat[i,2]/u0
    }
    
    if (length(which(u_mat[,3]<0)) == 0){
      u <- 1
    }
    else{
      UB_u <- u_mat[min(which(u_mat[,3]<0)),1]
      usol <- uniroot(ueqn,c(LB_u,UB_u),tol = 1e-10)
      u <- usol$root
    }
    v <- 1-Tvalue+Tvalue*u
    Rinfty <- 1-PGFG0(v,Pk)
  }
  
  eta <- 0
  m <- max(Pk[,1])
  for (j in c(0:m)) {
    eta <- eta+Pk[,2][which(Pk[,1]==j)]*(v^j)
  }
  
  OutputDF <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,u,v,eta)
  return(OutputDF)
}


#Percolation Process for Partial Immunity

# 2 Wave Result, the first wave is always typical process with naive population
# Input is output data frame from previous round


# Polarized Immunity Assumption

# uninfected degree distribution after first outbreak
P2u_polar <- function(alpha, Pk1, Result1){
  k <- 0
  l <- 0
  u1 <- Result1$u
  plk <-0
  plku <- data.frame(k,l,plk)
  m <- max(Pk1$kvalue)
  for (i in c(0:m)) {
    for (j in c(0:i)) {
      c <- (i^2+i)/2+j+1
      plku[c,1] <- i
      plku[c,2] <- j
      plku[c,3] <- choose(i,j)*(u1+(1-u1)*alpha)^j*((1-u1)*(1-alpha))^(i-j)
    }
  }
  return(plku)
}



# re-suscepitible degree distribution after first outbreak
P2s_polar <- function(alpha, Pk1, Result1, sigma, tau){
  k <- 0
  l <- 0
  u1 <- Result1$u
  plk <- 1
  plks <- data.frame(k,l,plk)
  m <- max(Pk1$kvalue)
  for (i in c(1:m)) {
    for (j in c(0:i)) {
      c <- (i^2+i)/2+j+1
      plks[c,1] <- i
      plks[c,2] <- j
      if(j==0){
        plks[c,3] <- (1-alpha)*tau^(i-1)
      }
      else if(j==i){
        plks[c,3] <- alpha*sigma^(i-1)
      }
      else{
        plks[c,3] <- alpha*choose(i-1,j-1)*sigma^(j-1)*tau^((i-1)-(j-1))+(1-alpha)*choose(i-1,j)*sigma^j*tau^(i-1-j)
      }
    }
  }
  return(plks)
}


# Residual Network degree distribution and residual size
Pl2_polar <- function(alpha, Pk1, Result1, ShowSepDist=FALSE){
  if(Result1$OBType=='Epidemic'){
    u1 <- Result1$u
    T1 <- Result1$Tvalue
    zeta <- Result1$v
    sigma <- u1*(1-T1)*(1-alpha)+alpha
    tau <- (1-u1+T1*u1)*(1-alpha)
    ST <- data.frame(sigma,tau)
    Pu <- P2u_polar(alpha,Pk1,Result1)
    Ps <- P2s_polar(alpha,Pk1,Result1,sigma,tau)
    m <- max(Pk1$kvalue)
    lvalue <- 0
    pl <- 0
    Pl2 <- data.frame(lvalue,pl)
    denom <- 0
    eta <- 0
    for (j in c(0:m)) {
      denom <- denom+Pk1$pk[which(Pk1$kvalue==j)]*(zeta^j+alpha*(1-zeta^j))
      eta <- eta+Pk1$pk[which(Pk1$kvalue==j)]*(zeta^j)
    }
    residualfrac=eta+alpha*(1-eta)
    for (l in c(0:m)) {
      numer <- 0
      for (k in c(l:m)) {
        numer <- numer+Pk1$pk[which(Pk1$kvalue==k)]*(zeta^k*Pu$plk[which(Pu$k==k & Pu$l==l)]+alpha*(1-zeta^k)*Ps$plk[which(Ps$k==k & Ps$l==l)])
      }
      Pl2[l+1,1] <- l
      Pl2[l+1,2] <- numer/denom
    }
    if(ShowSepDist==TRUE){
      return(list(ResidualDDist=Pl2,ResidualFrac=residualfrac, UninfectDist=Pu, ResusDist=Ps, ResusParam=ST,Eta=eta))
    }
    else{
      return(list(ResidualDDist=Pl2,ResidualFrac=residualfrac))
    }
  }
  else{
    return('First outbreak must be epidemic!')
  }
}
# ResidualFrac is the fraction of residual network size after polarized immunity devided by original network size
# The second endemic outbreak proportion size S is based on residual network size=ResidualFrac*N
# the limited outbreak size <s> is number of vertices, and by definition independent with network size



# Overall proportion of infected among whole population during the second outbreak
TypProc_polar <- function(Pk1,T1,alpha,T2,tol=1e-5){
  Result1 <- TypProc(Pk1,T1,tol)
  if(Result1$OBType=='Epidemic'){
    ResNetwork <- Pl2_polar(alpha, Pk1, Result1)
    Pl2 <- ResNetwork$ResidualDDist
    ResFrac <- ResNetwork$ResidualFrac
    Result2 <- TypProc(Pl2,T2,tol)
    Result1[2,]<- Result2[1,]
    if(Result2$OBType=='Epidemic'){
      RinftyAll <- Result2$Rinfty*ResFrac
      return(list(ResultOfOutbreaks=Result1,ResidualFrac=ResFrac,OverallInfectProb=RinftyAll))
    }
    else{
      RinftyAll <- Result2$Rinfty*ResFrac
      return(list(ResultofOutbreaks=Result1,ResidualFrac=ResFrac,OverallInfectProb=RinftyAll))
    }
  }
  else{
    return(list(Resultof1stOB=Result1, Resultof2ndOB='First Outbreak is not epidemic!'))
  }
}

#Here the overall Infect prob, considering the residual size, is corresponding to original network.





# Leaky 
# Two Type Percolation

# p_ij
P2_leaky <- function(Pk1,Result1){
  m <- max(Pk1$kvalue)
  u1 <- Result1$u
  T1 <- Result1$Tvalue
  zeta <- Result1$v
  P2_matrix <- matrix(0,m+1,m+1)
  denom <- 0
  for (k in c(0:m)) {
    denom <- denom+Pk1$pk[which(Pk1$kvalue==k)]*(zeta^k)
  }
  for (i in c(0:m)) {
    for (j in c(0:m)) {
      if (i+j>m){
        P2_matrix[i+1,j+1] <- 0
      }
      else{
        kk <- i+j
        P2_matrix[i+1,j+1] <- (Pk1$pk[which(Pk1$kvalue==(kk))]*(zeta^(kk))*choose(kk,i)*(u1^i)*(1-u1)^j)/denom
      }
    }
  }
  return(P2_matrix)
}

Q2_leaky<- function(Pk1,Result1){
  m <- max(Pk1$kvalue)
  u1 <- Result1$u
  T1 <- Result1$Tvalue
  zeta <- Result1$v
  Q2_matrix <- matrix(0,m+1,m+1)
  denom <- 0
  for (k in c(0:m)) {
    denom <- denom+Pk1$pk[which(Pk1$kvalue==k)]*(1-zeta^k)
  }
  for (i in c(0:m)) {
    for (j in c(0:m)) {
      if (i+j>m){
        Q2_matrix[i+1,j+1] <- 0
      }
      else if(j==0){
        Q2_matrix[i+1,j+1] <- 0
      }
      else{
        kk <- i+j
        Q2_matrix[i+1,j+1] <- (Pk1$pk[which(Pk1$kvalue==(kk))]*(1-zeta^(kk))*choose(kk-1,i)*((u1*(1-T1))^i)*(1-u1+u1*T1)^(j-1))/denom
      }
    }
  }
  return(Q2_matrix)
}

Tmatrix <- function(alpha, T2){
  output <- matrix(0,nrow = 2, ncol = 2)
  output[1,] <- c(T2, alpha*T2)
  output[2,] <- c(alpha*T2, alpha^2*T2)
  return(output)
}

#1-4
###Bansal Meyers version: extra T

### New Version
PGFdf1_leaky <- function(P2,Q2,m){
  dfAdx <- 0
  dfAdy <- 0
  dfBdx <- 0
  dfBdy <- 0
  
  for (i in c(0:m)) {
    for (j in c(0:m)) {
      dfAdx <- dfAdx+P2[i+1,j+1]*i
      dfAdy <- dfAdy+P2[i+1,j+1]*j
      dfBdx <- dfBdx+Q2[i+1,j+1]*i
      dfBdy <- dfBdy+Q2[i+1,j+1]*j
    }
  }
  output <- c(dfAdx,dfAdy,dfBdx,dfBdy)
  return(output)
}

#5-12
### Version without T
PGFdf2_leaky <- function(P2,Q2,Tmatrix,m){
  TAA <- Tmatrix[1,1]
  TAB <- Tmatrix[1,2]
  TBA <- Tmatrix[2,1]
  TBB <- Tmatrix[2,2]
  
  dfAAdx <- 0
  dfAAdy <- 0
  dfBAdx <- 0
  dfBAdy <- 0
  dfABdx <- 0
  dfABdy <- 0
  dfBBdx <- 0
  dfBBdy <- 0
  
  denom_AA <- 0
  denom_BA <- 0
  denom_AB <- 0
  denom_BB <- 0
  
  for (i in c(0:m)) {
    for (j in c(0:m)) {
      denom_AA <- denom_AA+i*P2[i+1,j+1]
      denom_BA <- denom_BA+j*P2[i+1,j+1]
      denom_AB <- denom_AB+i*Q2[i+1,j+1]
      denom_BB <- denom_BB+j*Q2[i+1,j+1]
      
      dfAAdx <- dfAAdx+P2[i+1,j+1]*i*(i-1)
      dfAAdy <- dfAAdy+P2[i+1,j+1]*i*j
      dfBAdx <- dfBAdx+P2[i+1,j+1]*i*j
      dfBAdy <- dfBAdy+P2[i+1,j+1]*j*(j-1)
      dfABdx <- dfABdx+Q2[i+1,j+1]*i*(i-1)
      dfABdy <- dfABdy+Q2[i+1,j+1]*i*j
      dfBBdx <- dfBBdx+Q2[i+1,j+1]*i*j
      dfBBdy <- dfBBdy+Q2[i+1,j+1]*j*(j-1)
    }
  }
  
  dfAAdx <- dfAAdx/denom_AA
  dfAAdy <- dfAAdy/denom_AA
  dfBAdx <- dfBAdx/denom_BA
  dfBAdy <- dfBAdy/denom_BA
  dfABdx <- dfABdx/denom_AB
  dfABdy <- dfABdy/denom_AB
  dfBBdx <- dfBBdx/denom_BB
  dfBBdy <- dfBBdy/denom_BB
  
  output <- c(dfAAdx,dfAAdy,dfBAdx,dfBAdy,dfABdx,dfABdy,dfBBdx,dfBBdy)
  return(output)
}

#Version with correction
s_leaky <- function(alpha,T2,Pk1,Result1,tol_chi=1e-2){
  m <- max(Pk1$kvalue)
  P2 <- P2_leaky(Pk1,Result1)
  Q2 <- Q2_leaky(Pk1,Result1)
  Tmatrix <- Tmatrix(alpha,T2)
  PGFdf1 <- PGFdf1_leaky(P2,Q2,m)
  PGFdf2 <- PGFdf2_leaky(P2,Q2,Tmatrix,m)
  
  U1 <- PGFG0(Result1$v,Pk1)
  V1 <- 1-U1
  K <- PGFd1G0(1,Pk1)
  
  TAA <- Tmatrix[1,1]
  TAB <- Tmatrix[1,2]
  TBA <- Tmatrix[2,1]
  TBB <- Tmatrix[2,2]
  
  dfAdx <- PGFdf1[1]
  dfAdy <- PGFdf1[2]
  dfBdx <- PGFdf1[3]
  dfBdy <- PGFdf1[4]
  
  dfAAdx <- PGFdf2[1]
  dfAAdy <- PGFdf2[2]
  dfBAdx <- PGFdf2[3]
  dfBAdy <- PGFdf2[4]
  dfABdx <- PGFdf2[5]
  dfABdy <- PGFdf2[6]
  dfBBdx <- PGFdf2[7]
  dfBBdy <- PGFdf2[8]
  
  f1 <- function(t,a){
    (1-t*(dfAAdx-dfABdx))*(-t*a*t*a*dfAAdy)
  }
  
  f2 <- function(t,a){
    (t*a^2*dfBAdy*dfBBdx-dfBAdx*(1-t*dfBBdy))
  }
  
  f3 <- function(t,a){
    (1-t*dfAAdx)*(1-t*a^2*dfBBdy)
  }
  
  chi <- f1(T2,alpha)*f2(T2,alpha)+f3(T2,alpha)
  
  chi_eqn <- function(x) {
    f1(x,alpha)*f2(x,alpha)+f3(x,alpha)
  }
  
  LB_chi <- 0
  UB_chi <- 1
  B_vec <- seq(from=LB_chi,to=UB_chi,by=tol_chi)
  B_mat <- matrix(0,nrow = length(B_vec),ncol = 3)
  B_mat[,1] <- B_vec
  chi0 <- chi_eqn(0) 
  
  for (i in c(1:length(B_vec))) {
    B_mat[i,2] <- chi_eqn(B_vec[i])
    B_mat[i,3] <- B_mat[i,2]/chi0
  }
  
  if (length(which(B_mat[,3]<0)) == 0){
    Tc <- 1
  }
  else{
    UB_chi <- B_mat[min(which(B_mat[,3]<0)),1]
    Tc_sol <- uniroot(chi_eqn,c(LB_chi,UB_chi))
    Tc <- Tc_sol$root
  }
  
  #Tc_sol <- uniroot(chi_eqn,c(0+tol,1-tol))
  #Tc <- Tc_sol$root
  
  # Tc here is threshold for T2, which makes denom chi=0
  # But what is T2c_leaky defined in BM's article, what is threshold for individual-level immunity model
  # Seems to be some average based on edge type
  
  e_UU <- U1*dfAdx/K
  e_UI <- U1*dfAdy/K
  e_IU <- V1*dfBdx/K 
  e_II <- V1*dfBdy/K
  
  T2c_leaky <- Tc*(e_UU+e_UI/alpha+e_IU/alpha+e_II/alpha^2)
  
  DFAAdx <- (1-TBB*dfBBdy)/chi  
  DFABdx <- (1-TAA*(dfAAdx-dfABdx))*(TBB*TBA*dfBAdy*dfBBdx+TBA*dfBAdx*(1-TBB*dfBBdy))/chi
  DFBAdx <- (1-TAA*(dfAAdx-dfABdx))*(1-TBB*dfBBdy)/chi
  DFBBdx <- (1-TAA*(dfAAdx-dfABdx))*(TBA*dfBBdx)/chi
  DFAAdy <- (1-TBB*(dfBBdy-dfBAdy))*(TAB*dfAAdy)/chi
  DFABdy <- (1-TBB*(dfBBdy-dfBAdy))*(1-TAA*dfAAdx)/chi
  DFBAdy <- (1-TBB*(dfBBdy-dfBAdy))*(TAB*dfAAdy)*(1-TAA*(dfAAdx-dfABdx))/chi
  DFBBdy <- ((1-TAA*dfAAdx)-(TAB*dfAAdy)*(TBA*(dfBAdx+dfBBdx)*(1-TAA*(dfAAdx-dfAAdy))))/chi
  
  sAA <- (TAA*dfAdx*DFAAdx+TAB*dfAdy*DFABdx)+1
  sAB <- (TAA*dfAdx*DFAAdy+TAB*dfAdy*DFABdy)
  sBA <- (TBA*dfBdx*DFBAdx+TBB*dfBdy*DFBBdx)
  sBB <- (TBA*dfBdx*DFBAdy+TBB*dfBdy*DFBBdy)+1
  
  sA <- sAA+sBA
  sB <- sAB+sBB
  s <- sA+sB
  output <- data.frame(sAA,sAB,sBA,sBB,chi,sA,sB,s,Tc,T2c_leaky)
  return(output)
}


# system of epidemic outbreak

fsys_soln <- function(P2,Q2,Tmatrix,m, tol=1e-5){
  TAA <- Tmatrix[1,1]
  TAB <- Tmatrix[1,2]
  TBA <- Tmatrix[2,1]
  TBB <- Tmatrix[2,2]
  
  #system
  fAA_T <- function(x,y){
    denom <- 0
    numer <- 0
    
    X <- 1+(x-1)*TAA
    Y <- 1+(y-1)*TAB
    
    for (i in c(0:m)) {
      for (j in c(0:m)) {
        numer <- numer+P2[i+1,j+1]*i*X^(i-1)*Y^(j)
        denom <- denom+P2[i+1,j+1]*i    
      }
    }
    output <- numer/denom
    return(output)
  }
  fBA_T <- function(x,y){
    denom <- 0
    numer <- 0
    
    X <- 1+(x-1)*TAA
    Y <- 1+(y-1)*TAB
    
    for (i in c(0:m)) {
      for (j in c(0:m)) {
        numer <- numer+P2[i+1,j+1]*j*X^(i)*Y^(j-1)
        denom <- denom+P2[i+1,j+1]*j    
      }
    }
    output <- numer/denom
    return(output)
  }
  fAB_T <- function(x,y){
    denom <- 0
    numer <- 0
    
    X <- 1+(x-1)*TBA
    Y <- 1+(y-1)*TBB
    
    for (i in c(0:m)) {
      for (j in c(0:m)) {
        numer <- numer+Q2[i+1,j+1]*i*X^(i-1)*Y^(j)
        denom <- denom+Q2[i+1,j+1]*i    
      }
    }
    output <- numer/denom
    return(output)
  }
  fBB_T <- function(x,y){
    denom <- 0
    numer <- 0
    
    X <- 1+(x-1)*TBA
    Y <- 1+(y-1)*TBB
    
    for (i in c(0:m)) {
      for (j in c(0:m)) {
        numer <- numer+Q2[i+1,j+1]*j*X^(i)*Y^(j-1)
        denom <- denom+Q2[i+1,j+1]*j    
      }
    }
    output <- numer/denom
    return(output)
  }
  
  fsys <- function(x){
    c(
      fAA_T(x[1],x[3])-x[1],
      fBA_T(x[1],x[3])-x[2],
      fAB_T(x[2],x[4])-x[3],
      fBB_T(x[2],x[4])-x[4]
    )
  }
  
  x0 <- c(0,0,0,0)
  soln <- fsolve(fsys,x0,tol=tol)
  
  return(soln$x)
}

ES_leaky <- function(alpha,T2,Pk1,Result1, tol=1e-8){
  m <- max(Pk1$kvalue)
  P2 <- P2_leaky(Pk1,Result1)
  Q2 <- Q2_leaky(Pk1,Result1)
  Tmatrix <- Tmatrix(alpha,T2)
  
  TAA <- Tmatrix[1,1]
  TAB <- Tmatrix[1,2]
  TBA <- Tmatrix[2,1]
  TBB <- Tmatrix[2,2]
  
  fsol <- fsys_soln(P2,Q2,Tmatrix,m,tol = tol)
  a <- fsol[1]
  b <- fsol[2]
  c <- fsol[3]
  d <- fsol[4]
  
  FAsol <- 0
  FBsol <- 0
  
  for (i in c(0:m)) {
    for (j in c(0:m)){
      FAsol <- FAsol+P2[i+1,j+1]*((1+(a-1)*TAA)^i)*((1+(c-1)*TAB)^j)
      FBsol <- FBsol+Q2[i+1,j+1]*((1+(b-1)*TBA)^i)*((1+(d-1)*TBB)^j)
    } 
  }
  SA <- 1-FAsol
  SB <- 1-FBsol
  
  U1 <- PGFG0(Result1$v,Pk1)
  V1 <- 1-U1
  
  S2_leaky <- SA*U1+SB*V1
  return(data.frame(SA,SB,S2_leaky,FAsol,FBsol, a, b, c, d, U1))
}

TypProc_leaky <- function(Pk1,T1,alpha,T2,tol=1e-8,ShowParam=FALSE){
  Result1 <- TypProc(Pk1,T1,tol)
  if(Result1$OBType=='Epidemic'){
    sframe <- s_leaky(alpha,T2,Pk1,Result1)
    ESframe <- ES_leaky(alpha,T2,Pk1,Result1,tol)
    Tvalue <- T2
    Tc_value <- sframe$Tc
    T2c_leaky <- sframe$T2c_leaky
    if(T2<Tc_value){
      OBType <- 'Limited'
      s <- sframe$s
      Rinfty <- 0
      Result2 <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,T2c_leaky)
      sAA <- sframe$sAA
      sBA <- sframe$sBA
      sAB <- sframe$sAB
      sBB <- sframe$sBB
      sA <- sframe$sA
      sB <- sframe$sB
      sTypeDetail <- data.frame(sA,sB,sAA,sBA,sAB,sBB)
      if(ShowParam==FALSE){
        return(list(Resultof1stOB=Result1,Resultof2ndOB=Result2))
      }
      else{
        return(list(Resultof1stOB=Result1,Resultof2ndOB=Result2,TypeDetail=sTypeDetail))
      }
    }
    else if(Tvalue==Tc_value){
      OBType <- 'Undefined'
      s <- 0
      Rinfty <- 0
      Result2 <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,T2c_leaky)
      return(list(Resultof1stOB=Result1,Resultof2ndOB=Result2))
    }
    else if(Tvalue>Tc_value){
      OBType <- 'Epidemic'
      s <- 0
      Rinfty <- ESframe$S2_leaky
      Result2 <- data.frame(Tvalue,Tc_value,OBType,s,Rinfty,T2c_leaky)
      uA <- ESframe$SA
      uB <- ESframe$SB
      uAA <- ESframe$a
      uBA <- ESframe$b
      uAB <- ESframe$c
      uBB <- ESframe$d
      uTypeDetail <- data.frame(uA,uB,uAA,uBA,uAB,uBB)
      if(ShowParam==FALSE){
        return(list(Resultof1stOB=Result1,Resultof2ndOB=Result2))
      }
      else{
        return(list(Resultof1stOB=Result1,Resultof2ndOB=Result2,TypeDetail=uTypeDetail))
      }
    }
  }
  else{
    return(list(Resultof1stOB=Result1, Resultof2ndOB='First Outbreak is not epidemic!'))
  }
}



#Miller Slim and Voltz
#Configuration model
library(deSolve)

ModProc_CM <- function(Pk, beta, gamma, init_theta=1e-3, ODEmaxTime=50, ODEstep=1e-2,ThetaTol=1e-9, TrackDyn=T){
  if (TrackDyn==T){
  Sys <- function(t, y, parms){
    with(as.list(c(parms,y)),{
      dtheta <- (-b)*theta+b*PGFd1G0(theta,Pk)/PGFd1G0(1,Pk)+g*(1-theta)
      dR <- g*(1-PGFG0(theta,Pk)-R)
      return(list(c(dtheta,dR))) 
    }) 
  }
  parms <- c(b=beta,g=gamma)
  times <- seq(0,ODEmaxTime,by=ODEstep)
  y <- c(theta=1-init_theta,R=0)
  
  Sys_out <- ode(y,times,Sys,parms)
  S_out <- PGFG0(Sys_out[,2],Pk)
  I_out <- 1-S_out-Sys_out[,3]
  Sys_out <- as.matrix(cbind(Sys_out,S_out,I_out))
  }
  
  g <- gamma
  b <- beta
  
  thetaEqn<- function(x) {
    g/(b+g)+b/(b+g)*PGFd1G0(x,Pk)/PGFd1G0(1,Pk)-x
  }
  
  LB_theta <- 0
  UB_theta <- 1
  step_theta <- 1e-2
  
  Btheta_vec <- seq(from=LB_theta,to=UB_theta,by=step_theta)
  Btheta_mat <- matrix(0,nrow = length(Btheta_vec),ncol = 3)
  Btheta_mat[,1] <- Btheta_vec
  theta0 <- thetaEqn(0) 
  
  for (i in c(1:length(Btheta_vec))) {
    Btheta_mat[i,2] <- thetaEqn(Btheta_vec[i])
    Btheta_mat[i,3] <- Btheta_mat[i,2]/theta0
  }
  
  if (length(which(Btheta_mat[,3]<0)) == 0){
    thetaInf <- 1
  }else{
    UB_theta <- Btheta_mat[min(which(Btheta_mat[,3]<0)),1]
    theta_sol <- uniroot(thetaEqn,c(LB_theta,UB_theta),tol = 1e-9)
    thetaInf <- theta_sol$root
  }
  
  R0 <- b/(b+g)*PGFd2G0(1,Pk)/PGFd1G0(1,Pk)
  RInf <- 1-PGFG0(thetaInf,Pk)
  
  if (TrackDyn==T){
    return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf, Dynamic=Sys_out))
  } else {
    return(list(R0=R0,RInfinity=RInf, ThetaInfinity=thetaInf))
  }
}
