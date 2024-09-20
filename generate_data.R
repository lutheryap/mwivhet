rm(list=ls())

## Function to generate judge data without covariates
# This DGP is used in Section 2 of the paper
# sig objects are covariance objects, S is the IV strength and Het the heterogeneity in TE
# K is the number of instruments (K+1 judges), c is the number of cases per judge
GenData_nocov <- function(S=3,Het=3,sigee=1,sigvv=1,sigexi=0,sigev=0,ConTE=FALSE,
                          beta=0,beta0=0,K=40,c=5) {
  # Calculations as function of parameters
  J <- K+1
  n <- J*c
  group <- as.vector(seq(1,J) %x% rep(1,c))
  
  # Create Z matrix of indicators
  Z <- matrix(0, nrow=length(group), ncol=length(unique(group)))
  Z[cbind(seq_along(group), group)] <- 1
  
  ## Calculations related to the projection matrix
  ZZ <- t(Z) %*% Z; ZZ_inv <- solve(ZZ)
  ZZ_inv2 <- chol(ZZ_inv); ZZZ_inv <- Z%*% ZZ_inv
  ZZZ_inv2 <- Z%*%ZZ_inv2 #n x k mx
  H <- rowSums(ZZZ_inv2^2) #n x 1 vector of P_ii
  P <- Z %*% ZZ_inv %*% t(Z); PP <- P*P
  
  h <- sqrt(Het/sqrt(K)/(c-1))
  s <- sqrt(S/sqrt(K)/(c-1))
  pi <- c(0,rep(-s,K/2),rep(s,K/2)) #balanced
  
  # Construct DGP Sigma
  sigxixi <- h/sigvv
  
  if (ConTE == TRUE) {
    Sigma <- matrix(c(sigee,sigev,sigev,sigvv),nrow=2)
  } else {
    sigxiv01 <- h
    sigxiv23 <- - h
    Sigma0 <- matrix(c(sigee,sigexi,sigev,sigexi,sigxixi,0,sigev,0,sigvv),nrow=3)
    #matrixcalc::is.positive.definite(Sigma0)
    Sigma01 <- matrix(c(sigee,sigexi,sigev,sigexi,sigxixi,sigxiv01,sigev,sigxiv01,sigvv),nrow=3)
    Sigma23 <- matrix(c(sigee,sigexi,sigev,sigexi,sigxixi,sigxiv23,sigev,sigxiv23,sigvv),nrow=3)
  }
  
  ## Construct dataframe
  df <- as.data.frame(group)
  df$pi <- as.vector(Z%*%pi)
  
  # Make random draws. 
  if (ConTE == TRUE) {
    vmat <- mvtnorm::rmvnorm(n,mean=rep(0,2),sigma=Sigma)
    df$eps <- vmat[,1]
    df$xi <- 0
    df$v <- vmat[,2]
  } else {
    vmat0 <- mvtnorm::rmvnorm(c,mean=rep(0,3),sigma=Sigma0)
    vmat01 <- mvtnorm::rmvnorm(K*c/2,mean=rep(0,3),sigma=Sigma01)
    vmat23 <- mvtnorm::rmvnorm(K*c/2,mean=rep(0,3),sigma=Sigma23)
    vmat <- rbind(vmat0,vmat01[1:(K*c/4),],vmat23,vmat01[((K*c/4)+1):(K*c/2),])
    df$eps <- vmat[,1]
    df$xi <- vmat[,2]
    df$v <- vmat[,3]
  }
  
  df$X <- df$pi + df$v
  df$Y <- df$X*(beta+df$xi) + df$eps
  df$e <- df$Y - df$X*beta0 
  M <- diag(n) - P
  df$MX <- c(M %*% df$X); df$Me <- c(M %*% df$e); df$MY <- M %*% df$Y
  
  df
}

## Function to generate data with covariates. Interact QOB with state.
# This DGP is used in Section 6.3 of the paper
# sig objects are covariance objects, S is the IV strength and Het the heterogeneity in TE
# K is the number of states, c is the number of cases per QOB times state
GenData_cov <- function(S=5,p1=7/8,Het=3,sigeps=.5,sigev=.1,beta=0,beta0=0,omega=.1,
                        K=20,c=5) {
  
  # Calculations as function of parameters
  J <- 2*(K+1)
  n <- J*c
  group <- as.vector(seq(1,J) %x% rep(1,c))
  
  # Create discrete W 
  groupW <- as.vector(seq(1,K+1) %x% rep(1,2*c))
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
  
  PP_off <- G*G
  
  h <- sqrt(Het/sqrt(K)/(c-1))
  s <- sqrt(S/sqrt(K)/(c-1))
  pi <- c(0,rep(s,K/2),rep(-s,K/2))
  
  # Generate gammad 
  gammad <- c(0,rep(omega,K/2),rep(-omega,K/2))
  
  ## Construct dataframe
  df <- as.data.frame(group)
  df$groupZ <- groupZ
  df$groupW <- groupW
  df$pi <- as.vector(Z%*%pi)
  df$gammad <- as.vector(Wd%*%gammad)
  dfbase <- df # This is fixed across simulations
  
  # Make random draws. This varies across sims
  df <- dfbase
  df$v <- runif(n,-1,1)
  df$hetdraw <- runif(n,0,1)
  df$vpos <- ifelse(df$v>0,1,-1)
  df$ppos <- p1; df$ppos[df$vpos <0] <- 1-p1
  
  df$eps <- sigev*df$vpos+rnorm(n,0,sigeps)
  df$xi <- c(rep(0,2*c),rep(c(rep(1,c),rep(0,c),rep(0,c),rep(0,c)),K/2))*ifelse(df$hetdraw<df$ppos,h,-h) +
    c(rep(0,2*c),rep(c(rep(0,c),rep(0,c),rep(1,c),rep(0,c)),K/2))*ifelse(df$hetdraw<df$ppos,-h,h)
  df$X <- ifelse(df$v < df$pi+df$gammad,1,-1) 
  df$Y <- df$X*(beta+df$xi)+df$gammad+df$eps
  df$e <- df$Y - df$X*beta0
  df$MX <- c(M %*% df$X)
  df$Me <- c(M %*% df$e)
  df$MY <- c(M %*% df$Y)
  df$dM <- diag(M)
  
  df
}



# Generate data
set.seed(2024)

# Judge data without covariates
df <- GenData_nocov()
# Keep only observable variables for dataset with no covariates (dnc)
dnc <- data.frame(group=df$group,X=df$X,Y=df$Y,e=df$e,MX=df$MX,Me=df$Me,MY=df$MY)
save(dnc,file="df_nocov.RData")

# Interacted data with covariates
df <- GenData_cov()
# Keep only observable variables for dataset with no covariates (dc)
dc <- data.frame(group=df$group,groupW=df$groupW,X=df$X,Y=df$Y,e=df$e,MX=df$MX,Me=df$Me,MY=df$MY)
save(dc,file="df_cov.RData")
