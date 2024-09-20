rm(list=ls())
source("generate_data.R")
source("L3Ofn.R")

# Demo for judge example without covariates 
load("df_nocov.RData")
# MX is the residual from regressing X on Z
# Me is the residual from regressing e on Z
# MY is the residual from regressing Y on Z
LMnum <- GetLM_nocov(df=dnc,X,e,groupZ=group) #numerator
LMVar <- L3Ovar_gloop_nocov(df=dnc,group,X,e,MX,Me) #L3O variance
LMt <- LMnum/sqrt(LMVar) #t statistic for LM test to compare with std normal
CI <- GetCIvals(GetCIcoef_nocov(df=dnc,groupZ=group,X,Y,MX,MY))
LMt; CI

# Demo for QOB example with covariates 
load("df_cov.RData")
# groupW are covariate groups (i.e. state indicators)
# group indexes every QOB and groupW combination
LMnum <- GetLM(df=dc,X,e,groupW,group) #numerator
LMVar <- L3Ovar_gloop_cov(df=dc,group,groupW,X,e,MX,Me) #L3O variance
LMt <- LMnum/sqrt(LMVar) #t statistic for LM test to compare with std normal
CI <- GetCIvals(GetCIcoef(df=dc,groupW,group,X,Y,MX,MY))
LMt; CI

# Demo for QOB example using stored G and P. Slower but allows continuous W, Z
# Construct G and P
GetGP <- function(group,groupW,n){
  groupZ <- groupW*(group%%2)
  # Create Z matrix of indicators
  Z <- matrix(0, nrow=length(groupZ), ncol=length(unique(groupZ))-1)
  Z[cbind(seq_along(groupZ), groupZ)] <- 1
  
  ## Calculations related to Z
  ZZ <- t(Z) %*% Z; ZZ_inv <- solve(ZZ)
  ZZ_inv2 <- chol(ZZ_inv); ZZZ_inv <- Z%*% ZZ_inv
  ZZZ_inv2 <- Z%*%ZZ_inv2 #n x k mx
  HZ <- Z%*% ZZ_inv %*% t(Z)
  
  ## Calculations related to W
  Wd <- matrix(0, nrow=length(groupW), ncol=length(unique(groupW)))
  Wd[cbind(seq_along(groupW), groupW)] <- 1
  #W <- cbind(Wd,Wcfix)
  W <- Wd
  WW <- t(W) %*% W; WW_inv <- solve(WW)
  WW_inv2 <- chol(WW_inv); WWW_inv <- W%*% WW_inv
  WWW_inv2 <- W%*%WW_inv2 
  HW <- W %*% WW_inv %*% t(W)
  MW <- diag(n) - HW
  
  ## Combine Z and W
  Q <- cbind(Z,W)
  QQ <- t(Q) %*% Q; QQ_inv <- solve(QQ)
  QQ_inv2 <- chol(QQ_inv); QQQ_inv <- Q%*% QQ_inv
  QQQ_inv2 <- Q%*%QQ_inv2 
  
  # P for full projection on Z and W
  P <- Q%*% QQ_inv %*% t(Q)
  M <- diag(n) - P
  
  ## G for UJIVE
  G <- solve(diag(n)-diag(diag(P))) %*% (P - diag(diag(P))) - 
    solve(diag(n) - diag(diag(HW))) %*% (HW - diag(diag(HW)))
  
  list(G=G,P=P,Z=Z,W=W)
}
GP <- GetGP(dc$group,dc$groupW,n=nrow(dc))
# Run the functions as before
LMnum <- dc$e %*% GP$G %*% dc$X
LMVar <- L3Ovar_iloop_cov(dc$X,dc$e,GP$P,GP$G) #L3O variance
LMt <- LMnum/sqrt(LMVar) #t statistic for LM test to compare with std normal
CI <- GetCIvals(GetCIcoef_iloop(df=dc,GP$P,GP$G,X,Y,MX,MY,GP$Z,GP$W,noisy=FALSE))
LMt; CI #Numerically identical to previous approach
