L3Ovar_iloop_cov <- function(X,e,P,G,noisy=FALSE) {
  n <- length(X)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Me <- M%*%e
  MX <- M%*%X
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Gicol <- matrix(G[,i],ncol=1) #force column vector
    Girow <- matrix(G[i,],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Gofficol <- Goff[,i]
    Goffirow <- matrix(Goff[i,],ncol=1)
    
    A11i <- t(X*Girow)%*%D2D3i%*%(Girow*X)*(t(Mi)%*%e)
    A12i <- t((Me)*X*Girow*Mi)%*%recD3i%*%(Girow*X*dM)-
      t(X*Girow*Mi)%*%(recD3i*M)%*%(Girow*X*(Me))
    A13i <- t(dM*X*Girow)%*%((onesN%x%t(Me))*recD3i)%*%(Girow*Mi*X)
    A14i <- t((Me)*X*Girow)%*%(M*recD3i)%*%(Gicol*Mi*X)
    A15i <- (t(Mi)%*%e)*(t(Goffirow^2*recD2i)%*%(dM*X^2))-
      (t(Goffirow^2*Mi*recD2i)%*%(Me*X^2))
    
    A21i <- t(X*Girow)%*%D2D3i%*%(Gicol*e)*(t(Mi)%*%X)
    A22i <- t((MX)*X*Girow*Mi)%*%recD3i%*%(Gicol*e*dM)-
      t(X*Girow*Mi)%*%(recD3i*M)%*%(Gicol*e*(MX))
    A23i <- t(dM*X*Girow)%*%((onesN%x%t(MX))*recD3i)%*%(Gicol*Mi*e)
    A24i <- t((MX)*X*Girow)%*%(M*recD3i)%*%(Gicol*Mi*e)
    A25i <- (t(Mi)%*%X)*(t(Goffirow*Gofficol*recD2i)%*%(dM*X*e))-
      (t(Goffirow*Gofficol*Mi*recD2i)%*%(MX*X*e))
    
    A31i <- t(e*Gicol)%*%D2D3i%*%(Gicol*e)*(t(Mi)%*%X)
    A32i <- t((MX)*e*Gicol*Mi)%*%recD3i%*%(Gicol*e*dM)-
      t(e*Gicol*Mi)%*%(recD3i*M)%*%(Gicol*e*(MX))
    A33i <- t(dM*e*Gicol)%*%((onesN%x%t(MX))*recD3i)%*%(Gicol*Mi*e)
    A34i <- t((MX)*e*Gicol)%*%(M*recD3i)%*%(Gicol*Mi*e)
    A35i <- (t(Mi)%*%X)*(t(Gofficol^2*recD2i)%*%(dM*e*e))-
      (t(Gofficol^2*Mi*recD2i)%*%(MX*e*e))
    
    A41i <- t(Gicol^2*e*dM*recD2i*Me)%*%(D2D3i)%*%(Mi*X) -
      t(Gicol^2*e*recD2i*Me)%*%(D2D3i*M*(onesN%x%t(X)))%*%(Mi)
    A42i <- t(Gicol^2*e*dM*recD2i)%*%(recD3i*M)%*%(Mi^2*X) -
      t(Gicol^2*e*Mi*recD2i)%*%(recD3i*M*M)%*%(Mi*X) -
      t(Gicol^2*e*dM*Mi*recD2i)%*%((onesN%x%t(dM))*recD3i)%*%(Mi*X) +
      t(Gicol^2*e*Mi^2*recD2i)%*%(recD3i*M)%*%(dM*X)
    A43i <- t(Gicol^2*e*dM*Mi*recD2i)%*%(recD3i)%*%(Mi^2*X*Me) -
      t(Gicol^2*e*Mi^2*recD2i)%*%(recD3i*M)%*%(Mi*X*Me) -
      M[i,i]*t(Gicol^2*e*dM*recD2i)%*%(recD3i*M)%*%(Mi*X*Me) +
      M[i,i]*t(Gicol^2*e*Mi*recD2i)%*%(recD3i*M*M)%*%(Me*X)
    A44i <- X[i]*M[i,i]*(t(Gofficol^2*recD2i)%*%(Me*e)) -
      X[i]*(t(Mi)%*%e)*(t(Gofficol^2*recD2i)%*%(Mi*e))
    
    A51i <- t(Girow*Gicol*e*dM*recD2i*MX)%*%(D2D3i)%*%(Mi*X) -
      t(Girow*Gicol*e*recD2i*MX)%*%(D2D3i*M*(onesN%x%t(X)))%*%(Mi)
    A52i <- t(Girow*Gicol*e*dM*recD2i)%*%(recD3i*M)%*%(Mi^2*X) -
      t(Girow*Gicol*e*Mi*recD2i)%*%(recD3i*M*M)%*%(Mi*X) -
      t(Girow*Gicol*e*dM*Mi*recD2i)%*%((onesN%x%t(dM))*recD3i)%*%(Mi*X) +
      t(Girow*Gicol*e*Mi^2*recD2i)%*%(recD3i*M)%*%(dM*X)
    A53i <- t(Girow*Gicol*e*dM*Mi*recD2i)%*%(recD3i)%*%(Mi^2*X*MX) -
      t(Girow*Gicol*e*Mi^2*recD2i)%*%(recD3i*M)%*%(Mi*X*MX) -
      M[i,i]*t(Girow*Gicol*e*dM*recD2i)%*%(recD3i*M)%*%(Mi*X*MX) +
      M[i,i]*t(Girow*Gicol*e*Mi*recD2i)%*%(recD3i*M*M)%*%(MX*X)
    A54i <- X[i]*M[i,i]*(t(Goffirow*Gofficol*recD2i)%*%(MX*e)) -
      X[i]*(t(Mi)%*%X)*(t(Goffirow*Gofficol*recD2i)%*%(Mi*e))
    
    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i*(M[i,]%*%e) + A43i + A44i
    A5vec[i] <- A51i + A52i*(M[i,]%*%X) + A53i + A54i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  
  sum(A1vec*e + 2*A2vec*e +A3vec*X -A4vec*X-A5vec*e)
}


bindrowvecs <- function(X,n) matrix(rep(X,n),nrow=n,byrow=TRUE)
bindcolvecs <- function(X,n) matrix(rep(X,n),nrow=n,byrow=FALSE)

ivectomats <- function(ds,X,g) {
  ds$group <- eval(substitute(group),ds)
  gidx <- which(ds$group==g)
  sumq <- sum(X[gidx])
  subtractq <- matrix(rep(0,nrow(ds)^2),nrow=nrow(ds))
  subtractq[gidx,] <- X[gidx]
  subtractq[,gidx] <-  subtractq[,gidx] + t(subtractq[gidx,])
  out <- (matrix(rep(1,nrow(ds)^2),nrow=nrow(ds)) - diag(nrow(ds)))*sumq -subtractq
  out
}

Getgroupindex <- function(ds,group) {
  ds$group <- eval(substitute(group),ds)
  ds$groupidx <- 0
  gvals <- unique(ds$group)
  for (g in unique(ds$group)) {
    ds$groupidx[ds$group==g] <- which(unique(ds$group)==g)
  }
  ds$groupidx
}

A1type_iloop_sum <- function(df,P,G,ipos,jpos,kpos,lpos,noisy=FALSE) {
  #A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,max(df$group))
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A1vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Gi <- matrix(G[i,],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Goffi <- Goff[i,]
    lposmx <- matrix(df$lpos,ncol=1)
    
    A11i <- t(df$jpos*Gi)%*%D2D3i%*%(Gi*df$kpos)*(df$lpos[i])
    A12i <- t(df$lpos*df$jpos*Gi*Mi)%*%recD3i%*%(Goffi*df$kpos*dM)-
      t(df$jpos*Gi*Mi)%*%(recD3i*M)%*%(Goffi*df$kpos*df$lpos)
    A13i <- t(dM*df$jpos*Gi)%*%((onesN%x%t(lposmx))*recD3i)%*%(Gi*Mi*df$kpos)
    A14i <- t((df$lpos)*df$jpos*Gi)%*%(M*recD3i)%*%(Gi*Mi*df$kpos)
    A15i <- (df$lpos[i])*(t(Goffi^2*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Goffi^2*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
   
    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  
  sum(A1vec*df$ipos)
}

A2type_iloop_sum <- function(df,P,G,ipos,jpos,kpos,lpos,noisy=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A2vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Girow <- matrix(G[i,],ncol=1) #force column vector
    Gicol <- matrix(G[,i],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Gofficol <- Goff[,i]
    Goffirow <- matrix(Goff[i,],ncol=1)
    lposmx <- matrix(df$lpos,ncol=1)
    
    A21i <- t(df$jpos*Girow)%*%D2D3i%*%(Gicol*df$kpos)*(df$lpos[i])
    A22i <- t(df$lpos*df$jpos*Girow*Mi)%*%recD3i%*%(Gicol*df$kpos*dM)-
      t(df$jpos*Girow*Mi)%*%(recD3i*M)%*%(Gicol*df$kpos*df$lpos)
    A23i <- t(dM*df$jpos*Girow)%*%((onesN%x%t(lposmx))*recD3i)%*%(Gicol*Mi*df$kpos)
    A24i <- t((df$lpos)*df$jpos*Girow)%*%(M*recD3i)%*%(Gicol*Mi*df$kpos)
    A25i <- (df$lpos[i])*(t(Gofficol*Goffirow*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Gofficol*Goffirow*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
    
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  
  sum(A2vec*df$ipos)
}

A3type_iloop_sum <- function(df,P,G,ipos,jpos,kpos,lpos,noisy=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A3vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Gi <- matrix(G[,i],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Goffi <- Goff[,i]
    lposmx <- matrix(df$lpos,ncol=1)
    
    A31i <- t(df$jpos*Gi)%*%D2D3i%*%(Gi*df$kpos)*(df$lpos[i])
    A32i <- t(df$lpos*df$jpos*Gi*Mi)%*%recD3i%*%(Goffi*df$kpos*dM)-
      t(df$jpos*Gi*Mi)%*%(recD3i*M)%*%(Goffi*df$kpos*df$lpos)
    A33i <- t(dM*df$jpos*Gi)%*%((onesN%x%t(lposmx))*recD3i)%*%(Gi*Mi*df$kpos)
    A34i <- t((df$lpos)*df$jpos*Gi)%*%(M*recD3i)%*%(Gi*Mi*df$kpos)
    A35i <- (df$lpos[i])*(t(Goffi^2*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Goffi^2*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
    
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  
  sum(A3vec*df$ipos)
}

A4type_iloop_sum <- function(df,P,G,ipos,jpos,kpos,lpos,noisy=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A4vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Gi <- matrix(G[,i],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Goffi <- Goff[,i]
    lposmx <- matrix(df$lpos,ncol=1)
    kposmx <- matrix(df$kpos,ncol=1)
    
    A41i <- t(Gi^2*df$jpos*dM*recD2i*df$lpos)%*%(D2D3i)%*%(Mi*df$kpos) -
      t(Gi^2*df$jpos*recD2i*df$lpos)%*%(D2D3i*M*(onesN%x%t(kposmx)))%*%(Mi)
    A42i <- t(Gi^2*df$jpos*dM*recD2i)%*%(recD3i*M)%*%(Mi^2*df$kpos) -
      t(Gi^2*df$jpos*Mi*recD2i)%*%(recD3i*M*M)%*%(Mi*df$kpos) -
      t(Gi^2*df$jpos*dM*Mi*recD2i)%*%((onesN%x%t(dM))*recD3i)%*%(Mi*df$kpos) +
      t(Gi^2*df$jpos*Mi^2*recD2i)%*%(recD3i*M)%*%(dM*df$kpos)
    A43i <- t(Gi^2*df$jpos*dM*Mi*recD2i)%*%(recD3i)%*%(Mi^2*df$kpos*df$lpos) -
      t(Gi^2*df$jpos*Mi^2*recD2i)%*%(recD3i*M)%*%(Mi*df$kpos*df$lpos) -
      M[i,i]*t(Gi^2*df$jpos*dM*recD2i)%*%(recD3i*M)%*%(Mi*df$kpos*df$lpos) +
      M[i,i]*t(Gi^2*df$jpos*Mi*recD2i)%*%(recD3i*M*M)%*%(df$lpos*df$kpos)
    A44i <- df$kpos[i]*M[i,i]*(t(Goffi^2*recD2i)%*%(df$lpos*df$jpos)) -
      df$kpos[i]*(df$lpos[i])*(t(Goffi^2*recD2i)%*%(Mi*df$jpos))
    
    
    A4vec[i] <- A41i + A42i*df$lpos[i] + A43i + A44i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  ret <- (A4vec*df$ipos) 
  
  sum(ret)
}

A5type_iloop_sum <- function(df,P,G,ipos,jpos,kpos,lpos,noisy=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))
  
  A5vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    Girow <- matrix(G[i,],ncol=1) #force column vector
    Gicol <- matrix(G[,i],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    Gofficol <- Goff[,i]
    Goffirow <- matrix(Goff[i,],ncol=1)
    lposmx <- matrix(df$lpos,ncol=1)
    kposmx <- matrix(df$kpos,ncol=1)
    
    A51i <- t(Girow*Gicol*df$jpos*dM*recD2i*df$lpos)%*%(D2D3i)%*%(Mi*df$kpos) -
      t(Girow*Gicol*df$jpos*recD2i*df$lpos)%*%(D2D3i*M*(onesN%x%t(kposmx)))%*%(Mi)
    A52i <- t(Girow*Gicol*df$jpos*dM*recD2i)%*%(recD3i*M)%*%(Mi^2*df$kpos) -
      t(Girow*Gicol*df$jpos*Mi*recD2i)%*%(recD3i*M*M)%*%(Mi*df$kpos) -
      t(Girow*Gicol*df$jpos*dM*Mi*recD2i)%*%((onesN%x%t(dM))*recD3i)%*%(Mi*df$kpos) +
      t(Girow*Gicol*df$jpos*Mi^2*recD2i)%*%(recD3i*M)%*%(dM*df$kpos)
    A53i <- t(Girow*Gicol*df$jpos*dM*Mi*recD2i)%*%(recD3i)%*%(Mi^2*df$kpos*df$lpos) -
      t(Girow*Gicol*df$jpos*Mi^2*recD2i)%*%(recD3i*M)%*%(Mi*df$kpos*df$lpos) -
      M[i,i]*t(Girow*Gicol*df$jpos*dM*recD2i)%*%(recD3i*M)%*%(Mi*df$kpos*df$lpos) +
      M[i,i]*t(Girow*Gicol*df$jpos*Mi*recD2i)%*%(recD3i*M*M)%*%(df$lpos*df$kpos)
    A54i <- df$kpos[i]*M[i,i]*(t(Gofficol*Goffirow*recD2i)%*%(df$lpos*df$jpos)) -
      df$kpos[i]*(df$lpos[i])*(t(Gofficol*Goffirow*recD2i)%*%(Mi*df$jpos))
    
    
    A5vec[i] <- A51i + A52i*df$lpos[i] + A53i + A54i
    
    if(noisy) {
      if (i%%10 == 0) cat(i/n," ")
    }
  }
  ret <- (A5vec*df$ipos) 
  
  sum(ret)
}

A1type_sum <- function(df,group,groupW,ipos,jpos,kpos,lpos,noisy=FALSE) {
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,max(df$group))
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW ==s,]

    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group==j] <- sum(ds$group==j)
    }
    ds <- ds[ds$numingrp>=3,]
    if (nrow(ds)==0) {
      for (g in unique(ds$group))  {
        A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow=length(ds$group), ncol=length(unique(ds$group)))
      ds$groupidx <- Getgroupindex(ds,group)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
      #ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))
      
      PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
      #PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)

      ZWmat <- matrix(1, nrow=length(ds$groupW), ncol=length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))
      
      # calculate values specific to this subset
      Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) -
            diag(1/(diag(diag(nrow(ds))) -pmin(diag(PW),.99)))%*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms),ncol=1)
      D2s <- dMs %*% t(dMs) - Ms*Ms
      recD2s <- 1/D2s; diag(recD2s) <- 0
      
      for (g in unique(ds$group)) {
        if (nrow(ds[ds$group==g,]) <= 3) {
          A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
        } else {
          repidx <- min(which(ds$group==g)) #representative index
          Pis <- matrix(Ps[,repidx],ncol=1)
          Pgs <- ifelse(Pis==0,0,1)%*%matrix(1,ncol=length(Pis),nrow=1)
          Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
          Mis <- matrix(-Ps[,repidx],ncol=1)
          recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
          D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
          D2D3is <- D2s/D3is; diag(D2D3is) <- 0
          recD3is <- 1/D3is; diag(recD3is) <- 0
          ones <- matrix(rep(1,nrow(ds)),ncol=1)
          Mes <- matrix(ds$lpos,ncol=1)
          
          A11vecs[g] <- t(ds$jpos*Gis)%*%(D2D3is*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(ds$kpos*Gis)
          A12vecs[g] <- t(ds$jpos*ds$lpos*Gis*Mis)%*%(recD3is*ivectomats(ds,ds$ipos,g))%*%(Gis*ds$kpos*dMs) -
            t(ds$jpos*Gis*Mis)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(Gis*ds$kpos*ds$lpos)
          A13vecs[g] <- t(ds$jpos*dMs*Gis)%*%(recD3is*(ones %*% t(Mes))*ivectomats(ds,ds$ipos,g))%*%(ds$kpos*Gis*Mis)
          A14vecs[g] <- t(ds$jpos*ds$lpos*Gis)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(ds$kpos*Gis*Mis)
          A15vecs[g] <- t(ds$lpos*ds$ipos)%*%(recD2s*Gs^2*Pgs)%*%(ds$jpos*ds$kpos*dMs) -
            t(ds$ipos)%*% (Gs^2*Ms*recD2s*Pgs) %*% (ds$lpos*ds$jpos*ds$kpos)
        }
      }
      
    }
    if(noisy) {
      cat(s, "of", max(df$groupW), "done. ")
    }
  }
  
  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) 
  
  sum(ret)
}

A4type_sum <- function(df,group,groupW,ipos,jpos,kpos,lpos,noisy=FALSE) {
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0,max(df$group))
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW ==s,]
    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group==j] <- sum(ds$group==j)
    }
    ds <- ds[ds$numingrp>=3,]
    if (nrow(ds)==0) {
      for (g in unique(ds$group))  {
        A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
      }
    }  else {
      ZQ <- matrix(0, nrow=length(ds$group), ncol=length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds,group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    #ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    #PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)
    
    ZWmat <- matrix(1, nrow=length(ds$groupW), ncol=length(ds$groupW))
    PW <- ZWmat / (length(ds$groupW))
    
    # calculate values specific to this subset
    Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) -
          diag(1/(diag(diag(nrow(ds))) -pmin(diag(PW),.99)))%*% (PW - diag(diag(PW)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms),ncol=1)
    D2s <- dMs %*% t(dMs) - Ms*Ms
    recD2s <- 1/D2s; diag(recD2s) <- 0
    
    for (g in unique(ds$group)) {
      if  (nrow(ds[ds$group==g,]) <= 3) {
        A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
      } else {
        repidx <- min(which(ds$group==g)) #representative index
        Pis <- matrix(Ps[,repidx],ncol=1)
        Pgs <- ifelse(Pis==0,0,1)%*%matrix(1,ncol=length(Pis),nrow=1)
        Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
        Mis <- matrix(-Ps[,repidx],ncol=1)
        recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
        D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
        D2D3is <- D2s/D3is; diag(D2D3is) <- 0
        recD3is <- 1/D3is; diag(recD3is) <- 0
        ones <- matrix(rep(1,nrow(ds)),ncol=1)
        Mes <- matrix(ds$lpos,ncol=1)
        
        A41vecs[g] <- t(Gis^2*ds$jpos*dMs*ds$lpos*recD2is)%*%(D2D3is*ivectomats(ds,ds$ipos,g))%*%(Mis*ds$kpos) -
          t(Gis^2*ds$jpos*ds$lpos*recD2is)%*%(D2D3is*Ms*(ones%x%t(ds$kpos))*ivectomats(ds,ds$ipos,g))%*%(Mis)
        A42vecs[g] <- t(Gis^2*ds$jpos*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis^2*ds$kpos) -
          t(Gis^2*ds$jpos*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis*ds$kpos) -
          t(Gis^2*ds$jpos*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis*ds$kpos) + 
          t(Gis^2*ds$jpos*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(dMs*ds$kpos)
        A43vecs[g] <- t(Gis^2*ds$jpos*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$ipos,g))%*%(Mis^2*ds$kpos*ds$lpos) -
          t(Gis^2*ds$jpos*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(Mis*ds$kpos*ds$lpos) -
          t(Gis^2*ds$jpos*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*dMs,g))%*%(Mis*ds$kpos*ds$lpos) +
          t(Gis^2*ds$jpos*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$ipos*dMs,g))%*%(ds$kpos*ds$lpos)
        A44vecs[g] <- t(ds$ipos*ds$kpos*dMs)%*%(recD2s*Gs^2*Pgs)%*%(ds$lpos*ds$jpos) -
          t(ds$ipos*ds$kpos*ds$lpos)%*%(recD2s*Gs^2*Pgs)%*%(Mis*ds$jpos)
      }
      
    }
    }
    if(noisy) {
      cat(s, "of", max(df$groupW), "done. ")
    }
  }
  
  ret <- (A41vecs + A42vecs +A43vecs +A44vecs) 
  
  sum(ret)
}


A1type_ijloop_sum <- function(df,ipos,jpos,kpos,lpos,IdPQ,IdPW,noisyi=FALSE,noisyj=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  #M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1,n),ncol=1)
  QQQinv <- Q%*%QQinv
  WWWinv <- W%*%WWinv
  
  
  A1vec <- rep(0,n)
  for (i in 1:n) {
    dPQi <- rep(0,n); dPQi[i] <- dPQ[i]
    dPWi <- rep(0,n); dPWi[i] <- dPW[i]
    
    Pi <- QQQinv %*% Q[i,]
    Gi <- (QQQinv %*% Q[i,] - dPQi)/dM[i] - 
      (WWWinv %*% W[i,] - dPWi)/dMW[i]
    Mi <- -Pi; Mi[i] <- 1-Pi[i]
    D2i <- dM*dM[i] - Mi^2
    recD2i <- 1/D2i; recD2i[i] <- 0
    
    A11i <- A12i <- A13i <- A14i <- rep(0,n)
    
    for (j in (1:n)[-i]) {
      dPQj <- rep(0,n); dPQj[j] <- dPQ[j]
      dPWj <- rep(0,n); dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j,]
      Mj <- -Pj; Mj[j] <- 1-Pj[j]
      D2j <- dM*dM[j] - Mj^2
      
      D3ij <- Mi[i]*D2j-(dM*Mi[j]^2+Mi^2*dM[j] - 2*Mj*Mi*Mi[j])
      recD3ij <- 1/D3ij; recD3ij[i] <- 0; recD3ij[j] <- 0
      D2D3ij <- D2j/D3ij; D2D3ij[i] <- 0; D2D3ij[j] <- 0
      
      A11ij <- (t(df$jpos*Gi)%*%D2D3ij)*(Gi[j]*df$kpos[j])*(df$lpos[i])
      A12ij <- (t(df$lpos*df$jpos*Gi*Mi)%*%recD3ij)*(Gi[j]*df$kpos[j]*dM[j])-
        (t(df$jpos*Gi*Mi)%*%(recD3ij*Mj))%*%(Gi[j]*df$kpos[j]*df$lpos[j])
      A13ij <- (t(dM*df$jpos*Gi)%*%(onesN*df$lpos[j]*recD3ij))*(Gi[j]*Mi[j]*df$kpos[j])
      A14ij <- (t((df$lpos)*df$jpos*Gi)%*%(Mj*recD3ij))*(Gi[j]*Mi[j]*df$kpos[j])
      
      A11i[j] <- A11ij; A12i[j] <- A12ij; A13i[j] <- A13ij; A14i[j] <- A14ij
      
      if (noisyj==TRUE & j%%10 == 0) cat("j:", j/n," ")
    }
    
    A1vec[i] <- sum(A11i - A12i - A13i + A14i)
    A15i <- (df$lpos[i])*(t(Gi^2*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Gi^2*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
    A1vec[i] <-  A1vec[i] + A15i
    
    if (noisyi==TRUE) cat("i:", i/n," ")
  }
  
  sum(A1vec*df$ipos)
}

A2type_ijloop_sum <- function(df,ipos,jpos,kpos,lpos,IdPQ,IdPW,noisyi=FALSE,noisyj=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  #M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1,n),ncol=1)
  QQQinv <- Q%*%QQinv
  WWWinv <- W%*%WWinv
  
  
  A2vec <- rep(0,n)
  for (i in 1:n) {
    dPQi <- rep(0,n); dPQi[i] <- dPQ[i]
    dPWi <- rep(0,n); dPWi[i] <- dPW[i]
    
    Pi <- QQQinv %*% Q[i,]
    Girow <- (QQQinv %*% Q[i,] - dPQi)/dM[i] - 
      (WWWinv %*% W[i,] - dPWi)/dMW[i]
    Gicol <- t(QQQinv %*% Q[i,] - dPQi)/dM - 
      t(WWWinv %*% W[i,] - dPWi)/dMW
    
    Mi <- -Pi; Mi[i] <- 1-Pi[i]
    D2i <- dM*dM[i] - Mi^2
    recD2i <- 1/D2i; recD2i[i] <- 0
    
    A21i <- A22i <- A23i <- A24i <- rep(0,n)
    
    for (j in (1:n)[-i]) {
      dPQj <- rep(0,n); dPQj[j] <- dPQ[j]
      dPWj <- rep(0,n); dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j,]
      Mj <- -Pj; Mj[j] <- 1-Pj[j]
      D2j <- dM*dM[j] - Mj^2
      
      D3ij <- Mi[i]*D2j-(dM*Mi[j]^2+Mi^2*dM[j] - 2*Mj*Mi*Mi[j])
      recD3ij <- 1/D3ij; recD3ij[i] <- 0; recD3ij[j] <- 0
      D2D3ij <- D2j/D3ij; D2D3ij[i] <- 0; D2D3ij[j] <- 0
      
      A21ij <- (t(df$jpos*Girow)%*%D2D3ij)*(Gicol[j]*df$kpos[j])*(df$lpos[i])
      A22ij <- (t(df$lpos*df$jpos*Girow*Mi)%*%recD3ij)*(Gicol[j]*df$kpos[j]*dM[j])-
        (t(df$jpos*Girow*Mi)%*%(recD3ij*Mj))%*%(Gicol[j]*df$kpos[j]*df$lpos[j])
      A23ij <- (t(dM*df$jpos*Girow)%*%(onesN*df$lpos[j]*recD3ij))*(Gicol[j]*Mi[j]*df$kpos[j])
      A24ij <- (t((df$lpos)*df$jpos*Girow)%*%(Mj*recD3ij))*(Gicol[j]*Mi[j]*df$kpos[j])
      
      A21i[j] <- A21ij; A22i[j] <- A22ij; A23i[j] <- A23ij; A24i[j] <- A24ij
      
      if (noisyj==TRUE & j%%10 == 0) cat("j:", j/n," ")
    }
    
    A2vec[i] <- sum(A21i - A22i - A23i + A24i)
    A25i <- (df$lpos[i])*(t(Girow*Gicol*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Girow*Gicol*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
    A2vec[i] <-  A2vec[i] + A25i
    
    if (noisyi==TRUE) cat("i:", i/n," ")
  }
  
  sum(A2vec*df$ipos)
}

A3type_ijloop_sum <- function(df,ipos,jpos,kpos,lpos,IdPQ,IdPW,noisyi=FALSE,noisyj=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  #M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1,n),ncol=1)
  QQQinv <- Q%*%QQinv
  WWWinv <- W%*%WWinv
  
  
  A3vec <- rep(0,n)
  for (i in 1:n) {
    dPQi <- rep(0,n); dPQi[i] <- dPQ[i]
    dPWi <- rep(0,n); dPWi[i] <- dPW[i]
    
    Pi <- QQQinv %*% Q[i,]
    Gi <- (QQQinv %*% Q[i,] - dPQi)/dM - 
      (WWWinv %*% W[i,] - dPWi)/dMW
    Mi <- -Pi; Mi[i] <- 1-Pi[i]
    D2i <- dM*dM[i] - Mi^2
    recD2i <- 1/D2i; recD2i[i] <- 0
    
    A31i <- A32i <- A33i <- A34i <- rep(0,n)
    
    for (j in (1:n)[-i]) {
      dPQj <- rep(0,n); dPQj[j] <- dPQ[j]
      dPWj <- rep(0,n); dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j,]
      Mj <- -Pj; Mj[j] <- 1-Pj[j]
      D2j <- dM*dM[j] - Mj^2
      
      D3ij <- Mi[i]*D2j-(dM*Mi[j]^2+Mi^2*dM[j] - 2*Mj*Mi*Mi[j])
      recD3ij <- 1/D3ij; recD3ij[i] <- 0; recD3ij[j] <- 0
      D2D3ij <- D2j/D3ij; D2D3ij[i] <- 0; D2D3ij[j] <- 0
      
      A31ij <- (t(df$jpos*Gi)%*%D2D3ij)*(Gi[j]*df$kpos[j])*(df$lpos[i])
      A32ij <- (t(df$lpos*df$jpos*Gi*Mi)%*%recD3ij)*(Gi[j]*df$kpos[j]*dM[j])-
        (t(df$jpos*Gi*Mi)%*%(recD3ij*Mj))%*%(Gi[j]*df$kpos[j]*df$lpos[j])
      A33ij <- (t(dM*df$jpos*Gi)%*%(onesN*df$lpos[j]*recD3ij))*(Gi[j]*Mi[j]*df$kpos[j])
      A34ij <- (t((df$lpos)*df$jpos*Gi)%*%(Mj*recD3ij))*(Gi[j]*Mi[j]*df$kpos[j])
      
      A31i[j] <- A31ij; A32i[j] <- A32ij; A33i[j] <- A33ij; A34i[j] <- A34ij
      
      if (noisyj==TRUE & j%%10 == 0) cat("j:", j/n," ")
    }
    
    A3vec[i] <- sum(A31i - A32i - A33i + A34i)
    A35i <- (df$lpos[i])*(t(Gi^2*recD2i)%*%(dM*df$jpos*df$kpos))-
      (t(Gi^2*Mi*recD2i)%*%(df$lpos*df$jpos*df$kpos))
    
    A3vec[i] <-  A3vec[i] + A35i
    
    if (noisyi==TRUE) cat("i:", i/n," ")
  }
  
  sum(A3vec*df$ipos)
}

A4type_ijloop_sum <- function(df,ipos,jpos,kpos,lpos,IdPQ,IdPW,noisyi=FALSE,noisyj=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  
  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1,n),ncol=1)
  QQQinv <- Q%*%QQinv
  WWWinv <- W%*%WWinv
  
  A4vec <- rep(0,n)
  for (i in 1:n) {
    dPQi <- rep(0,n); dPQi[i] <- dPQ[i]
    dPWi <- rep(0,n); dPWi[i] <- dPW[i]
    Pi <- QQQinv %*% Q[i,]
    Gi <- (QQQinv %*% Q[i,] - dPQi)/dM - 
      (WWWinv %*% W[i,] - dPWi)/dMW
    Mi <- -Pi; Mi[i] <- 1-Pi[i]
    D2i <- dM*dM[i] - Mi^2
    recD2i <- 1/D2i; recD2i[i] <- 0
    
    A41i <- A42i <- A43i <- rep(0,n)
    
    for (j in (1:n)[-i]) {
      dPQj <- rep(0,n); dPQj[j] <- dPQ[j]
      dPWj <- rep(0,n); dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j,]
      Mj <- -Pj; Mj[j] <- 1-Pj[j]
      D2j <- dM*dM[j] - Mj^2
      
      D3ij <- Mi[i]*D2j-(dM*Mi[j]^2+Mi^2*dM[j] - 2*Mj*Mi*Mi[j])
      recD3ij <- 1/D3ij; recD3ij[i] <- 0; recD3ij[j] <- 0
      D2D3ij <- D2j/D3ij; D2D3ij[i] <- 0; D2D3ij[j] <- 0
      
      A41ij <- (t(Gi^2*df$jpos*dM*recD2i*df$lpos)%*%(D2D3ij))*(Mi[j]*df$kpos[j]) -
        (t(Gi^2*df$jpos*recD2i*df$lpos)%*%(D2D3ij*Mj*(onesN*df$kpos[j])))*Mi[j]
      A42ij <- (t(Gi^2*df$jpos*dM*recD2i)%*%(recD3ij*Mj))*(Mi[j]^2*df$kpos[j]) -
        (t(Gi^2*df$jpos*Mi*recD2i)%*%(recD3ij*Mj^2))*(Mi[j]*df$kpos[j]) -
        (t(Gi^2*df$jpos*dM*Mi*recD2i)%*%((onesN*dM[j])*recD3ij))*(Mi[j]*df$kpos[j]) +
        (t(Gi^2*df$jpos*Mi^2*recD2i)%*%(recD3ij*Mj))%*%(dM[j]*df$kpos[j])
      A43ij <- (t(Gi^2*df$jpos*dM*Mi*recD2i)%*%(recD3ij))*(Mi[j]^2*df$kpos[j]*df$lpos[j]) -
        (t(Gi^2*df$jpos*Mi^2*recD2i)%*%(recD3ij*Mj))*(Mi[j]*df$kpos[j]*df$lpos[j]) -
        Mi[i]*(t(Gi^2*df$jpos*dM*recD2i)%*%(recD3ij*Mj))*(Mi[j]*df$kpos[j]*df$lpos[j]) +
        Mi[i]*(t(Gi^2*df$jpos*Mi*recD2i)%*%(recD3ij*Mj^2))*(df$lpos[j]*df$kpos[j])
      
      A41i[j] <- A41ij; A42i[j] <- A42ij; A43i[j] <- A43ij
      
      if (noisyj==TRUE & j%%10 == 0) cat("j:", j/n," ")
    }
    
    A44i <- df$kpos[i]*Mi[i]*(t(Gi^2*recD2i)%*%(df$lpos*df$jpos)) -
      df$kpos[i]*(df$lpos[i])*(t(Gi^2*recD2i)%*%(Mi*df$jpos))
    A4vec[i] <- sum(A41i + A42i*df$lpos + A43i) + A44i
    
    if (noisyi==TRUE) cat("i:", i/n," ")
  }
  
  sum(A4vec*df$ipos)
}

A5type_ijloop_sum <- function(df,ipos,jpos,kpos,lpos,IdPQ,IdPW,noisyi=FALSE,noisyj=FALSE) {
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  n <- nrow(df)
  
  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1,n),ncol=1)
  QQQinv <- Q%*%QQinv
  WWWinv <- W%*%WWinv
  
  A5vec <- rep(0,n)
  for (i in 1:n) {
    dPQi <- rep(0,n); dPQi[i] <- dPQ[i]
    dPWi <- rep(0,n); dPWi[i] <- dPW[i]
    Pi <- QQQinv %*% Q[i,]
    Girow <- (QQQinv %*% Q[i,] - dPQi)/dM[i] - 
      (WWWinv %*% W[i,] - dPWi)/dMW[i]
    Gicol <- t(QQQinv %*% Q[i,] - dPQi)/dM - 
      t(WWWinv %*% W[i,] - dPWi)/dMW
    Mi <- -Pi; Mi[i] <- 1-Pi[i]
    D2i <- dM*dM[i] - Mi^2
    recD2i <- 1/D2i; recD2i[i] <- 0
    
    A51i <- A52i <- A53i <- rep(0,n)
    
    for (j in (1:n)[-i]) {
      dPQj <- rep(0,n); dPQj[j] <- dPQ[j]
      dPWj <- rep(0,n); dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j,]
      Mj <- -Pj; Mj[j] <- 1-Pj[j]
      D2j <- dM*dM[j] - Mj^2
      
      D3ij <- Mi[i]*D2j-(dM*Mi[j]^2+Mi^2*dM[j] - 2*Mj*Mi*Mi[j])
      recD3ij <- 1/D3ij; recD3ij[i] <- 0; recD3ij[j] <- 0
      D2D3ij <- D2j/D3ij; D2D3ij[i] <- 0; D2D3ij[j] <- 0
      
      A51ij <- (t(Girow*Gicol*df$jpos*dM*recD2i*df$lpos)%*%(D2D3ij))*(Mi[j]*df$kpos[j]) -
        (t(Girow*Gicol*df$jpos*recD2i*df$lpos)%*%(D2D3ij*Mj*(onesN*df$kpos[j])))*Mi[j]
      A52ij <- (t(Girow*Gicol*df$jpos*dM*recD2i)%*%(recD3ij*Mj))*(Mi[j]^2*df$kpos[j]) -
        (t(Girow*Gicol*df$jpos*Mi*recD2i)%*%(recD3ij*Mj^2))*(Mi[j]*df$kpos[j]) -
        (t(Girow*Gicol*df$jpos*dM*Mi*recD2i)%*%((onesN*dM[j])*recD3ij))*(Mi[j]*df$kpos[j]) +
        (t(Girow*Gicol*df$jpos*Mi^2*recD2i)%*%(recD3ij*Mj))%*%(dM[j]*df$kpos[j])
      A53ij <- (t(Girow*Gicol*df$jpos*dM*Mi*recD2i)%*%(recD3ij))*(Mi[j]^2*df$kpos[j]*df$lpos[j]) -
        (t(Girow*Gicol*df$jpos*Mi^2*recD2i)%*%(recD3ij*Mj))*(Mi[j]*df$kpos[j]*df$lpos[j]) -
        Mi[i]*(t(Girow*Gicol*df$jpos*dM*recD2i)%*%(recD3ij*Mj))*(Mi[j]*df$kpos[j]*df$lpos[j]) +
        Mi[i]*(t(Girow*Gicol*df$jpos*Mi*recD2i)%*%(recD3ij*Mj^2))*(df$lpos[j]*df$kpos[j])
      
      A51i[j] <- A51ij; A52i[j] <- A52ij; A53i[j] <- A53ij
      
      if (noisyj==TRUE & j%%10 == 0) cat("j:", j/n," ")
    }
    
    A54i <- df$kpos[i]*Mi[i]*(t(Girow*Gicol*recD2i)%*%(df$lpos*df$jpos)) -
      df$kpos[i]*(df$lpos[i])*(t(Girow*Gicol*recD2i)%*%(Mi*df$jpos))
    A5vec[i] <- sum(A51i + A52i*df$lpos + A53i) + A54i
    
    if (noisyi==TRUE) cat("i:", i/n," ")
  }
  
  sum(A5vec*df$ipos)
}



A1type_sum_nocov <- function(df,groupZ,ipos,jpos,kpos,lpos,noisy=FALSE) {
  df$groupZ <- eval(substitute(groupZ),df)
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,max(df$groupZ))
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ ==s,]
    ZQ <- matrix(0, nrow=length(ds$groupZ), ncol=length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds,groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    
    # calculate values specific to this subset
    Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) 
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms),ncol=1)
    D2s <- dMs %*% t(dMs) - Ms*Ms
    recD2s <- 1/D2s; diag(recD2s) <- 0
    
    g=s
    
    repidx <- min(which(ds$groupZ==g)) #representative index
    Pis <- matrix(Ps[,repidx],ncol=1)
    Pgs <- ifelse(Pis==0,0,1)%*%matrix(1,ncol=length(Pis),nrow=1)
    Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
    Mis <- matrix(-Ps[,repidx],ncol=1)
    recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
    D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
    D2D3is <- D2s/D3is; diag(D2D3is) <- 0
    recD3is <- 1/D3is; diag(recD3is) <- 0
    ones <- matrix(rep(1,nrow(ds)),ncol=1)
    Mes <- matrix(ds$lpos,ncol=1)
    
    A11vecs[g] <- t(ds$jpos*Gis)%*%(D2D3is*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(ds$kpos*Gis)
    A12vecs[g] <- t(ds$jpos*ds$lpos*Gis*Mis)%*%(recD3is*ivectomats(ds,ds$ipos,g))%*%(Gis*ds$kpos*dMs) -
      t(ds$jpos*Gis*Mis)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(Gis*ds$kpos*ds$lpos)
    A13vecs[g] <- t(ds$jpos*dMs*Gis)%*%(recD3is*(ones %*% t(Mes))*ivectomats(ds,ds$ipos,g))%*%(ds$kpos*Gis*Mis)
    A14vecs[g] <- t(ds$jpos*ds$lpos*Gis)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(ds$kpos*Gis*Mis)
    A15vecs[g] <- t(ds$lpos*ds$ipos)%*%(recD2s*Gs^2*Pgs)%*%(ds$jpos*ds$kpos*dMs) -
      t(ds$ipos)%*% (Gs^2*Ms*recD2s*Pgs) %*% (ds$lpos*ds$jpos*ds$kpos)
    if(noisy) {
      cat(s, "of", max(df$groupZ), "done. ")
    }
  }
  
  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) 
  
  sum(ret)
}

A4type_sum_nocov <- function(df,groupZ,ipos,jpos,kpos,lpos,noisy=FALSE) {
  df$groupZ <- eval(substitute(groupZ),df)
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0,max(df$groupZ))
  
  df$ipos <- eval(substitute(ipos),df)
  df$jpos <- eval(substitute(jpos),df)
  df$kpos <- eval(substitute(kpos),df)
  df$lpos <- eval(substitute(lpos),df)
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ ==s,]
    ZQ <- matrix(0, nrow=length(ds$groupZ), ncol=length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds,groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    
    # calculate values specific to this subset
    Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) 
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms),ncol=1)
    D2s <- dMs %*% t(dMs) - Ms*Ms
    recD2s <- 1/D2s; diag(recD2s) <- 0
    
    g=s
    
    repidx <- min(which(ds$groupZ==g)) #representative index
    Pis <- matrix(Ps[,repidx],ncol=1)
    Pgs <- ifelse(Pis==0,0,1)%*%matrix(1,ncol=length(Pis),nrow=1)
    Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
    Mis <- matrix(-Ps[,repidx],ncol=1)
    recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
    D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
    D2D3is <- D2s/D3is; diag(D2D3is) <- 0
    recD3is <- 1/D3is; diag(recD3is) <- 0
    ones <- matrix(rep(1,nrow(ds)),ncol=1)
    Mes <- matrix(ds$lpos,ncol=1)
    
    A41vecs[g] <- t(Gis^2*ds$jpos*dMs*ds$lpos*recD2is)%*%(D2D3is*ivectomats(ds,ds$ipos,g))%*%(Mis*ds$kpos) -
      t(Gis^2*ds$jpos*ds$lpos*recD2is)%*%(D2D3is*Ms*(ones%x%t(ds$kpos))*ivectomats(ds,ds$ipos,g))%*%(Mis)
    A42vecs[g] <- t(Gis^2*ds$jpos*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis^2*ds$kpos) -
      t(Gis^2*ds$jpos*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis*ds$kpos) -
      t(Gis^2*ds$jpos*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(Mis*ds$kpos) + 
      t(Gis^2*ds$jpos*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*ds$lpos,g))%*%(dMs*ds$kpos)
    A43vecs[g] <- t(Gis^2*ds$jpos*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$ipos,g))%*%(Mis^2*ds$kpos*ds$lpos) -
      t(Gis^2*ds$jpos*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos,g))%*%(Mis*ds$kpos*ds$lpos) -
      t(Gis^2*ds$jpos*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$ipos*dMs,g))%*%(Mis*ds$kpos*ds$lpos) +
      t(Gis^2*ds$jpos*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$ipos*dMs,g))%*%(ds$kpos*ds$lpos)
    A44vecs[g] <- t(ds$ipos*ds$kpos*dMs)%*%(recD2s*Gs^2*Pgs)%*%(ds$lpos*ds$jpos) -
      t(ds$ipos*ds$kpos*ds$lpos)%*%(recD2s*Gs^2*Pgs)%*%(Mis*ds$jpos)
    if(noisy) {
      cat(s, "of", max(df$groupZ), "done. ")
    }
  }
  
  ret <- (A41vecs + A42vecs +A43vecs +A44vecs) 
  
  sum(ret)
}


GetLM <- function(df,X,e,groupW,group,noisy=FALSE) {
  LMvecs <- rep(0,max(df$groupW))
  df$Xpos <- eval(substitute(X),df)
  df$epos <- eval(substitute(e),df)
  df$groupW <- eval(substitute(groupW),df)
  df$group <- eval(substitute(group),df)
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW ==s,]
    ZQ <- matrix(0, nrow=length(ds$group), ncol=length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds,group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    #ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))
    ZWmat <- matrix(1, nrow=length(ds$groupW), ncol=length(ds$groupW))
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    #PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)
    PW <- ZWmat / (length(ds$groupW))
    
    # calculate values specific to this subset
    Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) -
          diag(1/(diag(diag(nrow(ds))) -pmin(diag(PW),.99)))%*% (PW - diag(diag(PW)))
    LMvecs[s] <- t(ds$epos) %*% Gs %*% ds$Xpos
    
    if(noisy) {
      cat(s, "of", max(df$groupW), "done. ")
    }
  }
  sum(LMvecs)
}

GetLM_nocov <- function(df,X,e,groupZ,noisy=FALSE) {
  df$groupZ <- eval(substitute(groupZ),df)
  LMvecs <- rep(0,max(df$groupZ))
  df$Xpos <- eval(substitute(X),df)
  df$epos <- eval(substitute(e),df)
  
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ ==s,]
    ZQ <- matrix(0, nrow=length(ds$groupZ), ncol=length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds,groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    
    # calculate values specific to this subset
    Gs <-  diag(1/(diag(diag(nrow(ds))) -pmin(diag(PQ),.99))) %*% (PQ - diag(diag(PQ))) 
    LMvecs[s] <- t(ds$epos) %*% Gs %*% ds$Xpos
    
    if(noisy) {
      cat(s, "of", max(df$groupZ), "done. ")
    }
  }
  sum(LMvecs)
}


GetCIcoef <- function(df,groupW,group,X,Y,MX,MY,q =qnorm(.975)^2,noisy=FALSE) {
  df$X <- eval(substitute(X),df)
  df$Y <- eval(substitute(Y),df)
  df$MX <- eval(substitute(MX),df)
  df$MY <- eval(substitute(MY),df)
  ## Calculate Components
  C0 <- A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
    2*A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    A1type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=Y,lpos=MX,noisy=noisy) -
    A4type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=X,lpos=MY,noisy=noisy) -
    A4type_sum(df,group,groupW,ipos=Y,jpos=Y,kpos=X,lpos=MX,noisy=noisy)
  C1 <- -(A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
            2*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
            A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) -
            A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) -
            A4type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) + 
            A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            2*A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            A1type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A4type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A4type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy))
  C2 <- A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    2*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy)
  PXY <- GetLM(df,X,Y,groupW,group)
  PXX <- GetLM(df,X,X,groupW,group)
  
  
  acon <- PXX^2- q*C2
  bcon <- -2*PXY*PXX-q*C1
  ccon <- PXY^2-q*C0
  
  c(acon,bcon,ccon)
}

GetSigMx <- function(df,groupW,group,X,Y,MX,MY,noisy=FALSE) {
  df$X <- eval(substitute(X),df)
  df$Y <- eval(substitute(Y),df)
  df$MX <- eval(substitute(MX),df)
  df$MY <- eval(substitute(MY),df)
  ## Calculate Components
  sig11 <- 4*A1type_sum(df,group,groupW,ipos=Y,jpos=Y,kpos=Y,lpos=MY,noisy=noisy) +
    2*A4type_sum(df,group,groupW,ipos=Y,jpos=Y,kpos=Y,lpos=MY,noisy=noisy) 
  sig22 <- A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
    2*A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    A1type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=Y,lpos=MX,noisy=noisy) +
    A4type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    A4type_sum(df,group,groupW,ipos=Y,jpos=Y,kpos=X,lpos=MX,noisy=noisy)
  sig33 <- 4*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    2*A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) 
  sig12 <- 2*A1type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=Y,lpos=MY,noisy=noisy) +
    2*A1type_sum(df,group,groupW,ipos=Y,jpos=X,kpos=Y,lpos=MY,noisy=noisy) +
    2*A4type_sum(df,group,groupW,ipos=X,jpos=Y,kpos=Y,lpos=MY,noisy=noisy) 
  sig23 <- 2*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
    2*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    2*A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy)
  sig13 <- 4*A1type_sum(df,group,groupW,ipos=X,jpos=X,kpos=Y,lpos=MY,noisy=noisy) +
    2*A4type_sum(df,group,groupW,ipos=X,jpos=X,kpos=Y,lpos=MY,noisy=noisy)
  
  c(sig11,sig22,sig33,sig12,sig23,sig13)
}

GetCIcoef_nocov <- function(df,groupZ,X,Y,MX,MY,q =qnorm(.975)^2,noisy=FALSE) {
  df$groupZ <- eval(substitute(groupZ),df)
  df$X <- eval(substitute(X),df)
  df$Y <- eval(substitute(Y),df)
  df$MX <- eval(substitute(MX),df)
  df$MY <- eval(substitute(MY),df)
  ## Calculate Components
  C0 <- A1type_sum_nocov(df,groupZ,ipos=Y,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
    2*A1type_sum_nocov(df,groupZ,ipos=Y,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    A1type_sum_nocov(df,groupZ,ipos=X,jpos=Y,kpos=Y,lpos=MX,noisy=noisy) -
    A4type_sum_nocov(df,groupZ,ipos=X,jpos=Y,kpos=X,lpos=MY,noisy=noisy) -
    A4type_sum_nocov(df,groupZ,ipos=Y,jpos=Y,kpos=X,lpos=MX,noisy=noisy)
  C1 <- -(A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
            2*A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
            A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) -
            A4type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) -
            A4type_sum_nocov(df,groupZ,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) + 
            A1type_sum_nocov(df,groupZ,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            2*A1type_sum_nocov(df,groupZ,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            A1type_sum_nocov(df,groupZ,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A4type_sum_nocov(df,groupZ,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A4type_sum_nocov(df,groupZ,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy))
  C2 <- A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    2*A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    A1type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A4type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A4type_sum_nocov(df,groupZ,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy)
  PXY <- GetLM_nocov(df,X,Y,groupZ)
  PXX <- GetLM_nocov(df,X,X,groupZ)
  
  
  acon <- PXX^2- q*C2
  bcon <- -2*PXY*PXX-q*C1
  ccon <- PXY^2-q*C0
  
  c(acon,bcon,ccon)
}

GetLM_WQ <- function(df,IdPW,IdPQ,dPW,dPQ,W,Q,X,Y) {
  df$Xpos <- eval(substitute(X),df)
  df$Ypos <- eval(substitute(Y),df)
  
  WW <- t(W) %*% W; WWinv <- solve(WW)
  QQ <- t(Q) %*% Q; QQinv <- solve(QQ) 
  QX <- t(Q) %*% df$Xpos; QY <- t(Q) %*% df$Ypos
  WX <- t(W) %*% df$Xpos; WY <- t(W) %*% df$Ypos
  
  QXdQ <- t(Q) %*% (df$Xpos/IdPQ); QYdQ <- t(Q) %*% (df$Ypos/IdPQ)
  WXdW <- t(W) %*% (df$Xpos/IdPW); WYdW <- t(W) %*% (df$Ypos/IdPW)
  
  
  (t(QXdQ) %*% QQinv %*% QY - sum(df$Xpos*dPQ/IdPQ*df$Ypos)) -
    (t(WXdW) %*% WWinv %*% WY - sum(df$Xpos*dPW/IdPW*df$Ypos))
}

GetCIcoef_iloop <- function(df,P,G,X,Y,MX,MY,Z,W,q =qnorm(.975)^2,noisy=FALSE) {
  df$X <- eval(substitute(X),df)
  df$Y <- eval(substitute(Y),df)
  df$MX <- eval(substitute(MX),df)
  df$MY <- eval(substitute(MY),df)
  ## Calculate Components
  C0 <- A1type_iloop_sum(df,P,G,ipos=Y,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
    2*A2type_iloop_sum(df,P,G,ipos=Y,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
    A3type_iloop_sum(df,P,G,ipos=X,jpos=Y,kpos=Y,lpos=MX,noisy=noisy) -
    A4type_iloop_sum(df,P,G,ipos=X,jpos=Y,kpos=X,lpos=MY,noisy=noisy) -
    A5type_iloop_sum(df,P,G,ipos=Y,jpos=Y,kpos=X,lpos=MX,noisy=noisy)
  C1 <- -(A1type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) +
            2*A2type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) +
            A3type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=Y,lpos=MX,noisy=noisy) -
            A4type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MY,noisy=noisy) -
            A5type_iloop_sum(df,P,G,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) + 
            A1type_iloop_sum(df,P,G,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            2*A2type_iloop_sum(df,P,G,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
            A3type_iloop_sum(df,P,G,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A4type_iloop_sum(df,P,G,ipos=X,jpos=Y,kpos=X,lpos=MX,noisy=noisy) -
            A5type_iloop_sum(df,P,G,ipos=Y,jpos=X,kpos=X,lpos=MX,noisy=noisy))
  C2 <- A1type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    2*A2type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) +
    A3type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A4type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy) -
    A5type_iloop_sum(df,P,G,ipos=X,jpos=X,kpos=X,lpos=MX,noisy=noisy)
  
  ## Generate dPQ
  n <- nrow(df)
  Q <- cbind(Z,W)
  dPQ <- dPW <- rep(0,n)
  WW <- t(W) %*% W; WWinv <- solve(WW)
  QQ <- t(Q) %*% Q; QQinv <- solve(QQ)
  for (i in 1:n) {
    dPW[i] <- W[i,] %*% WWinv %*% W[i,]
    dPQ[i] <- Q[i,] %*% QQinv %*% Q[i,]
  }
  IdPW <- pmax(1-dPW,.01); IdPQ <- pmax(1-dPQ,.01)
  
  PXY <- GetLM_WQ(df,IdPW,IdPQ,dPW,dPQ,W,Q,X,Y)
  PXX <- GetLM_WQ(df,IdPW,IdPQ,dPW,dPQ,W,Q,X,X)
  
  
  acon <- PXX^2- q*C2
  bcon <- -2*PXY*PXX-q*C1
  ccon <- PXY^2-q*C0
  
  c(acon,bcon,ccon)
}

GetCIcoef_ijloop <- function(df,IdPW,IdPQ,dPW,dPQ,X,Y,MX,MY,q =qnorm(.975)^2,noisyi=FALSE,noisyj=FALSE) {
  df$X <- eval(substitute(X),df)
  df$Y <- eval(substitute(Y),df)
  df$MX <- eval(substitute(MX),df)
  df$MY <- eval(substitute(MY),df)
  ## Calculate Components
  C0 <- A1type_ijloop_sum(df,ipos=Y,jpos=X,kpos=X,lpos=MY,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
    2*A2type_ijloop_sum(df,ipos=Y,jpos=X,kpos=Y,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
    A3type_ijloop_sum(df,ipos=X,jpos=Y,kpos=Y,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
    A4type_ijloop_sum(df,ipos=X,jpos=Y,kpos=X,lpos=MY,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
    A5type_ijloop_sum(df,ipos=Y,jpos=Y,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj)
  C1 <- -(A1type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MY,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
            2*A3type_ijloop_sum(df,ipos=X,jpos=X,kpos=Y,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
            A4type_ijloop_sum(df,ipos=X,jpos=X,kpos=Y,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
            A4type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MY,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
            A5type_ijloop_sum(df,ipos=X,jpos=Y,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) + 
            A1type_ijloop_sum(df,ipos=Y,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
            2*A2type_ijloop_sum(df,ipos=Y,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
            A3type_ijloop_sum(df,ipos=X,jpos=Y,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
            A4type_ijloop_sum(df,ipos=X,jpos=Y,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
            A5type_ijloop_sum(df,ipos=Y,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj))
  C2 <- A1type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
    2*A2type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) +
    A3type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
    A4type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj) -
    A5type_ijloop_sum(df,ipos=X,jpos=X,kpos=X,lpos=MX,IdPQ,IdPW,noisyi=noisyi,noisyj=noisyj)
  PXY <- GetLM_WQ(df,IdPW,IdPQ,dPW,dPQ,W,Q,X,Y)
  PXX <- GetLM_WQ(df,IdPW,IdPQ,dPW,dPQ,W,Q,X,X)
  
  
  acon <- PXX^2- q*C2
  bcon <- -2*PXY*PXX-q*C1
  ccon <- PXY^2-q*C0
  
  c(acon,bcon,ccon)
}


GetCIvals <- function(CIcoef) {
  acon <- CIcoef[1]; bcon <- CIcoef[2]; ccon <- CIcoef[3]
  det <- CIcoef[2]^2-4*CIcoef[1]*CIcoef[3]
  CILB <- (-bcon-sqrt(det))/(2*acon); CIUB <- (-bcon+sqrt(det))/(2*acon)
  c(CILB,CIUB)
}

GetCItypebd <- function(CIcoef) {
  det <- CIcoef[2]^2-4*CIcoef[1]*CIcoef[3]
  if (CIcoef[1]>=0 & det >=0) {
    CItype <- 1 #convex interval
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1]<0 & det >=0) {
    CItype <- 2 #donut
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1]<0 & det <0) {
    CItype <- 3 #accept everything
    CIbounds <- c(-100,100)
  } else {
    CItype <- 4 #reject everything
    CIbounds <- c(-100,100)
  }
  c(CItype,CIbounds)
}


L3Ovar_gloop_cov <- function(df,group,groupW,X,e,MX,Me) {
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0,length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0,length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0,length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0,length(unique(df$group)))
  for (s in 1:length(unique(df$groupW))) {
    ds <- df[df$groupW ==s,]
    ZQ <- matrix(0, nrow=length(ds$group), ncol=length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds,group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)
    
    # calculate values specific to this subset
    Gs <-  solve(diag(nrow(ds)) - diag(diag(PQ))) %*% (PQ - diag(diag(PQ))) -
      solve(diag(nrow(ds)) - diag(diag(PW)))%*% (PW - diag(diag(PW)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms),ncol=1)
    D2s <- dMs %*% t(dMs) - Ms*Ms
    recD2s <- 1/D2s; diag(recD2s) <- 0
    #diag0 <- matrix(1,nrow=nrow(ds),ncol=nrow(ds)); diag(diag0) <- 0
    
    for (g in unique(ds$group)) {
      repidx <- min(which(ds$group==g)) #representative index
      Pis <- matrix(Ps[,repidx],ncol=1)
      Pgs <- ifelse(Pis==0,0,1)%*%matrix(1,ncol=length(Pis),nrow=1)
      Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
      Mis <- matrix(-Ps[,repidx],ncol=1)
      recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
      D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
      D2D3is <- D2s/D3is; diag(D2D3is) <- 0
      recD3is <- 1/D3is; diag(recD3is) <- 0
      ones <- matrix(rep(1,nrow(ds)),ncol=1)
      Mes <- matrix(ds$Me,ncol=1)
      MXs <- matrix(ds$MX,ncol=1)
      es <- matrix(ds$e,ncol=1)
      
      A11vecs[g] <- t(ds$X*Gis)%*%(D2D3is*ivectomats(ds,ds$e*ds$Me,g))%*%(ds$X*Gis)
      A12vecs[g] <- t(ds$X*ds$Me*Gis*Mis)%*%(recD3is*ivectomats(ds,ds$e,g))%*%(Gis*ds$X*dMs) -
        t(ds$X*Gis*Mis)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(Gis*ds$X*ds$Me)
      A13vecs[g] <- t(ds$X*dMs*Gis)%*%(recD3is*(ones %*% t(Mes))*ivectomats(ds,ds$e,g))%*%(ds$X*Gis*Mis)
      A14vecs[g] <- t(ds$X*ds$Me*Gis)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(ds$X*Gis*Mis)
      A15vecs[g] <- t(ds$Me*ds$e)%*%(recD2s*Gs^2*Pgs)%*%(ds$X^2*dMs) -
        t(ds$e)%*% (Gs^2*Ms*recD2s*Pgs) %*% (ds$Me*ds$X^2)
      
      A21vecs[g] <- t(ds$X*Gis)%*%(D2D3is*ivectomats(ds,ds$e*ds$MX,g))%*%(ds$e*Gis)
      A22vecs[g] <- t(ds$X*ds$MX*Gis*Mis)%*%(recD3is*ivectomats(ds,ds$e,g))%*%(Gis*ds$e*dMs) -
        t(ds$X*Gis*Mis)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(Gis*ds$e*ds$MX)
      A23vecs[g] <- t(ds$X*dMs*Gis)%*%(recD3is*(ones %*% t(MXs))*ivectomats(ds,ds$e,g))%*%(ds$e*Gis*Mis)
      A24vecs[g] <- t(ds$X*ds$MX*Gis)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(ds$e*Gis*Mis)
      A25vecs[g] <- t(ds$MX*ds$e)%*%(recD2s*Gs^2*Pgs)%*%(ds$X*ds$e*dMs) -
        t(ds$e)%*% (Gs^2*Ms*recD2s*Pgs) %*% (ds$MX*ds$X*ds$e)
      
      A31vecs[g] <- t(ds$e*Gis)%*%(D2D3is*ivectomats(ds,ds$X*ds$MX,g))%*%(ds$e*Gis)
      A32vecs[g] <- t(ds$e*ds$MX*Gis*Mis)%*%(recD3is*ivectomats(ds,ds$X,g))%*%(Gis*ds$e*dMs) -
        t(ds$e*Gis*Mis)%*%(recD3is*Ms*ivectomats(ds,ds$X,g))%*%(Gis*ds$e*ds$MX)
      A33vecs[g] <- t(ds$e*dMs*Gis)%*%(recD3is*(ones %*% t(MXs))*ivectomats(ds,ds$X,g))%*%(ds$e*Gis*Mis)
      A34vecs[g] <- t(ds$e*ds$MX*Gis)%*%(recD3is*Ms*ivectomats(ds,ds$X,g))%*%(ds$e*Gis*Mis)
      A35vecs[g] <- t(ds$MX*ds$X)%*%(recD2s*Gs^2*Pgs)%*%(ds$e^2*dMs) -
        t(ds$X)%*% (Gs^2*Ms*recD2s*Pgs) %*% (ds$MX*ds$e^2)
      
      A41vecs[g] <- t(Gis^2*ds$e*dMs*ds$Me*recD2is)%*%(D2D3is*ivectomats(ds,ds$X,g))%*%(Mis*ds$X) -
        t(Gis^2*ds$e*ds$Me*recD2is)%*%(D2D3is*Ms*(ones%x%t(ds$X))*ivectomats(ds,ds$X,g))%*%(Mis)
      A42vecs[g] <- t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis^2*ds$X) -
        t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis*ds$X) -
        t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis*ds$X) + 
        t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*ds$Me,g))%*%(dMs*ds$X)
      A43vecs[g] <- t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$X,g))%*%(Mis^2*ds$X*ds$Me) -
        t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X,g))%*%(Mis*ds$X*ds$Me) -
        t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*dMs,g))%*%(Mis*ds$X*ds$Me) +
        t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$X*dMs,g))%*%(ds$X*ds$Me)
      A44vecs[g] <- t(ds$X^2*dMs)%*%(recD2s*Gs^2*Pgs)%*%(ds$Me*ds$e) -
        t(ds$X^2*ds$Me)%*%(recD2s*Gs^2*Pgs)%*%(Mis*ds$e)
      
      A51vecs[g] <- t(Gis^2*ds$e*dMs*ds$MX*recD2is)%*%(D2D3is*ivectomats(ds,ds$e,g))%*%(Mis*ds$X) -
        t(Gis^2*ds$e*ds$MX*recD2is)%*%(D2D3is*Ms*(ones%x%t(ds$X))*ivectomats(ds,ds$e,g))%*%(Mis)
      A52vecs[g] <- t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis^2*ds$X) -
        t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis*ds$X) -
        t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis*ds$X) + 
        t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*ds$MX,g))%*%(dMs*ds$X)
      A53vecs[g] <- t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$e,g))%*%(Mis^2*ds$X*ds$MX) -
        t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(Mis*ds$X*ds$MX) -
        t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*dMs,g))%*%(Mis*ds$X*ds$MX) +
        t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$e*dMs,g))%*%(ds$X*ds$MX)
      A54vecs[g] <- t(ds$X*ds$e*dMs)%*%(recD2s*Gs^2*Pgs)%*%(ds$MX*ds$e) -
        t(ds$X*ds$e*ds$MX)%*%(recD2s*Gs^2*Pgs)%*%(Mis*ds$e)
      
      #if (g%%10 == 0) cat(g/length(unique(df$group))," ")
    }
    
  }
  
  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) + 
    2*(A21vecs - A22vecs - A23vecs +A24vecs +A25vecs) +
    (A31vecs - A32vecs - A33vecs +A34vecs +A35vecs) -
    (A41vecs + A42vecs +A43vecs +A44vecs) -
    (A51vecs + A52vecs +A53vecs +A54vecs)
  
  sum(ret)
}

L3Ovar_gloop_nocov <- function(df,group,X,e,MX,Me) {
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0,length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0,length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0,length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0,length(unique(df$group)))
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0,length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0,length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0,length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0,length(unique(df$group)))
  for (g in 1:length(unique(df$group))) {
    ds <- df[df$group ==g,]
    ZQ <- matrix(0, nrow=length(ds$group), ncol=length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds,group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    
    PQ <- ZQ %*% solve(t(ZQ)%*% ZQ) %*% t(ZQ)
    
    # calculate values specific to this subset
    Gs <-  PQ
    Ps <- PQ
    Psoff <- Ps -diag(nrow(ds))*diag(Ps)
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms),ncol=1)
    D2s <- dMs %*% t(dMs) - Ms*Ms
    recD2s <- 1/D2s; diag(recD2s) <- 0
    #diag0 <- matrix(1,nrow=nrow(ds),ncol=nrow(ds)); diag(diag0) <- 0
    
    repidx <- min(which(ds$group==g)) #representative index
    Gis <- matrix(Gs[,repidx],ncol=1); Gis[repidx,1] <- Gis[repidx+1,1] #Put the G back
    Mis <- matrix(-Ps[,repidx],ncol=1)
    recD2is <- matrix(recD2s[,repidx],ncol=1); recD2is[repidx,1] <- recD2is[repidx+1,1] #Put the G back
    D3is <- Ms[repidx,repidx]*D2s-(dMs%*%t(Mis)^2+Mis^2%*%t(dMs)-2*Ms*(Mis%*%t(Mis)))
    D2D3is <- D2s/D3is; diag(D2D3is) <- 0
    recD3is <- 1/D3is; diag(recD3is) <- 0
    ones <- matrix(rep(1,nrow(ds)),ncol=1)
    Mes <- matrix(ds$Me,ncol=1)
    MXs <- matrix(ds$MX,ncol=1)
    es <- matrix(ds$e,ncol=1)
    
    A11vecs[g] <- t(ds$X*Gis)%*%(D2D3is*ivectomats(ds,ds$e*ds$Me,g))%*%(ds$X*Gis)
    A12vecs[g] <- t(ds$X*ds$Me*Gis)%*%(Ms*Psoff*recD3is*ivectomats(ds,ds$e,g))%*%(ds$X*dMs) -
      t(ds$X*Gis)%*%(recD3is*Ms*Ms*Psoff*ivectomats(ds,ds$e,g))%*%(ds$X*ds$Me)
    A13vecs[g] <- t(ds$X*dMs*Gis)%*%(recD3is*(ones %*% t(Mes))*ivectomats(ds,ds$e,g))%*%(ds$X*Gis*Mis)
    A14vecs[g] <- t(ds$X*ds$Me*Gis)%*%(recD3is*Ms*Ms*ivectomats(ds,ds$e,g))%*%(ds$X*Gis)
    A15vecs[g] <- t(ds$Me*ds$e)%*%(recD2s*Psoff^2)%*%(ds$X^2*dMs) -
      t(ds$e)%*% (Ms*recD2s*Psoff^2) %*% (ds$Me*ds$X^2)
    
    A21vecs[g] <- t(ds$X*Gis)%*%(D2D3is*ivectomats(ds,ds$e*ds$MX,g))%*%(ds$e*Gis)
    A22vecs[g] <- t(ds$X*ds$MX*Gis)%*%(Ms*Psoff*recD3is*ivectomats(ds,ds$e,g))%*%(ds$e*dMs) -
      t(ds$X*Gis)%*%(recD3is*Ms*Ms*Psoff*ivectomats(ds,ds$e,g))%*%(ds$e*ds$MX)
    A23vecs[g] <- t(ds$X*dMs*Gis)%*%(recD3is*(ones %*% t(MXs))*ivectomats(ds,ds$e,g))%*%(ds$e*Gis*Mis)
    A24vecs[g] <- t(ds$X*ds$MX*Gis)%*%(recD3is*Ms*Ms*ivectomats(ds,ds$e,g))%*%(ds$e*Gis)
    A25vecs[g] <- t(ds$MX*ds$e)%*%(recD2s*Psoff^2)%*%(ds$X*ds$e*dMs) -
      t(ds$e)%*% (Ms*recD2s*Psoff^2) %*% (ds$MX*ds$X*ds$e)
    
    A31vecs[g] <- t(ds$e*Gis)%*%(D2D3is*ivectomats(ds,ds$X*ds$MX,g))%*%(ds$e*Gis)
    A32vecs[g] <- t(ds$e*ds$MX*Gis)%*%(Ms*Psoff*recD3is*ivectomats(ds,ds$X,g))%*%(ds$e*dMs) -
      t(ds$e*Gis)%*%(recD3is*Ms*Ms*Psoff*ivectomats(ds,ds$X,g))%*%(ds$e*ds$MX)
    A33vecs[g] <- t(ds$e*dMs*Gis)%*%(recD3is*(ones %*% t(MXs))*ivectomats(ds,ds$X,g))%*%(ds$e*Gis*Mis)
    A34vecs[g] <- t(ds$e*ds$MX*Gis)%*%(recD3is*Ms*Ms*ivectomats(ds,ds$X,g))%*%(ds$e*Gis)
    A35vecs[g] <- t(ds$MX*ds$X)%*%(recD2s*Psoff^2)%*%(ds$e^2*dMs) -
      t(ds$X)%*% (Ms*recD2s*Psoff^2) %*% (ds$MX*ds$e^2)
    
    A41vecs[g] <- t(Gis^2*ds$e*dMs*ds$Me)%*%(recD3is*ivectomats(ds,ds$X,g))%*%(Mis*ds$X) -
      t(Gis^2*ds$e*ds$Me)%*%(recD3is*Ms*(ones%x%t(ds$X))*ivectomats(ds,ds$X,g))%*%(Mis)
    A42vecs[g] <- t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis^2*ds$X) -
      t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis*ds$X) -
      t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$X*ds$Me,g))%*%(Mis*ds$X) + 
      t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*ds$Me,g))%*%(dMs*ds$X)
    A43vecs[g] <- t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$X,g))%*%(Mis^2*ds$X*ds$Me) -
      t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X,g))%*%(Mis*ds$X*ds$Me) -
      t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$X*dMs,g))%*%(Mis*ds$X*ds$Me) +
      t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$X*dMs,g))%*%(ds$X*ds$Me)
    A44vecs[g] <- t(ds$X^2*dMs)%*%(recD2s*Psoff^2)%*%(ds$Me*ds$e) -
      t(ds$X^2*ds$Me)%*%(recD2s*Psoff^2)%*%(Mis*ds$e)
    
    A51vecs[g] <- t(Gis^2*ds$e*dMs*ds$MX)%*%(recD3is*ivectomats(ds,ds$e,g))%*%(Mis*ds$X) -
      t(Gis^2*ds$e*ds$MX)%*%(recD3is*Ms*(ones%x%t(ds$X))*ivectomats(ds,ds$e,g))%*%(Mis)
    A52vecs[g] <- t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis^2*ds$X) -
      t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis*ds$X) -
      t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*(ones %*% t(dMs))*ivectomats(ds,ds$e*ds$MX,g))%*%(Mis*ds$X) + 
      t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*ds$MX,g))%*%(dMs*ds$X)
    A53vecs[g] <- t(Gis^2*ds$e*dMs*Mis*recD2is)%*%(recD3is*ivectomats(ds,ds$e,g))%*%(Mis^2*ds$X*ds$MX) -
      t(Gis^2*ds$e*Mis^2*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e,g))%*%(Mis*ds$X*ds$MX) -
      t(Gis^2*ds$e*dMs*recD2is)%*%(recD3is*Ms*ivectomats(ds,ds$e*dMs,g))%*%(Mis*ds$X*ds$MX) +
      t(Gis^2*ds$e*Mis*recD2is)%*%(recD3is*Ms^2*ivectomats(ds,ds$e*dMs,g))%*%(ds$X*ds$MX)
    A54vecs[g] <- t(ds$X*ds$e*dMs)%*%(recD2s*Psoff^2)%*%(ds$MX*ds$e) -
      t(ds$X*ds$e*ds$MX)%*%(recD2s*Psoff^2)%*%(Mis*ds$e)
    
  }
  
  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) + 
    2*(A21vecs - A22vecs - A23vecs +A24vecs +A25vecs) +
    (A31vecs - A32vecs - A33vecs +A34vecs +A35vecs) -
    (A41vecs + A42vecs +A43vecs +A44vecs) -
    (A51vecs + A52vecs +A53vecs +A54vecs)
  
  sum(ret)
}

L3Ovar_block <- function(X,e,P,c) {
  n <- length(X)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Me <- M%*%e
  MX <- M%*%X
  Poff <- P - diag(diag(P))
  
  D3block <- ((c-1)^3-(c-1)*3-2)/c^3*(matrix(rep(1,c^2),nrow=c)-diag(c))
  recD3block <- c^3/((c-1)^3-(c-1)*3-2)*(matrix(rep(1,c^2),nrow=c)-diag(c))
  D2D3block <- c*((c-1)^2-1)/((c-1)^3-(c-1)*3-2)*(matrix(rep(1,c^2),nrow=c)-diag(c))
  recD2block <- c^2/((c-1)^2-1)*(matrix(rep(1,c^2),nrow=c)-diag(c))
  
  Pvec <- rep(1/c,n)
  Mvec <- rep(-1/c,n)
  
  D2D3blockrep <- as.matrix(Matrix::bdiag(replicate(J,D2D3block,simplify=FALSE)))
  recD3blockrep <- as.matrix(Matrix::bdiag(replicate(J,recD3block,simplify=FALSE)))
  recD2blockrep <- as.matrix(Matrix::bdiag(replicate(J,recD2block,simplify=FALSE)))
  blockrep <- Z %*% t(Z) -diag(n)
  ZtZ <- Z %*% t(Z)
  
  bindrowvec <- function(X) matrix(rep(X,n),nrow=n,byrow=TRUE)
  bindcolvec <- function(X) matrix(rep(X,n),nrow=n,byrow=FALSE)
  ivectomat <- function(X) (bindcolvec(ZtZ %*% X)*blockrep)-((bindrowvec(X)+bindcolvec(X))*blockrep)
  
  ## A1 
  A11 <- t(X*Pvec) %*%(D2D3blockrep*ivectomat(Me*e)) %*% (Pvec*X)
  A12 <- t(Me*X*Pvec)%*%(recD3blockrep*Poff*M*ivectomat(e))%*%(X*dM)-
    t(X*Pvec)%*%(recD3blockrep*M*M*Poff*ivectomat(e))%*%(X*(Me))
  A13 <- t(dM*X*Pvec)%*%(recD3blockrep*(onesN%x%t(Me))*M*ivectomat(e))%*%(X*Pvec)
  A14 <- t(Me*X*Pvec)%*%(recD3blockrep*M*M*ivectomat(e))%*%(X*Pvec)
  A15 <- t(Me*e)%*%((Poff^2*recD2)%*% (dM*X^2)) - t((Poff^2*M*recD2) %*% (Me*X^2))%*%e
  ## A2 
  A21 <- t(X*Pvec)%*%(D2D3blockrep*ivectomat(MX*e))%*%(Pvec*e)
  A22 <- t(MX*X*Pvec)%*%(M*recD3blockrep*ivectomat(e))%*%(Pvec*e*dM) -
    t(X*Pvec)%*%(M*M*recD3blockrep*ivectomat(e))%*%(Pvec*e*MX)
  A23 <- t(dM*X*Pvec)%*%(recD3blockrep*(onesN%x%t(MX))*M*ivectomat(e))%*%(e*Pvec)
  A24 <- t(MX*X*Pvec)%*%(recD3blockrep*M*M*ivectomat(e))%*%(e*Pvec)
  A25 <- t(MX*e)%*%((Poff^2*recD2)%*% (dM*X*e)) - t((Poff^2*M*recD2) %*% (MX*X*e))%*%e
  ## A3 
  A31 <- t(e*Pvec)%*%(D2D3blockrep*ivectomat(MX*X))%*%(Pvec*e)
  A32 <- t(MX*e*Pvec)%*%(M*recD3blockrep*ivectomat(X))%*%(Pvec*e*dM) -
    t(e*Pvec)%*%(M*M*recD3blockrep*ivectomat(X))%*%(Pvec*e*MX)
  A33 <- t(dM*e*Pvec)%*%(recD3blockrep*(onesN%x%t(MX))*M*ivectomat(X))%*%(e*Pvec)
  A34 <- t(MX*e*Pvec)%*%(recD3blockrep*M*M*ivectomat(X))%*%(e*Pvec)
  A35 <- t(MX*X)%*%((Poff^2*recD2)%*% (dM*e*e)) - t((Poff^2*M*recD2) %*% (MX*e*e))%*%X
  
  ## A4 
  A41 <- t(Pvec^2*e*dM*Me)%*%(recD3blockrep*ivectomat(X)*M)%*%(X) -
    t(Pvec^2*e*Me)%*%(recD3blockrep*M^2*(onesN%x%t(X))*ivectomat(X))%*%(onesN)
  A42 <- t(Pvec^2*e*dM)%*%(recD2*recD3blockrep*M^3*ivectomat(X*Me+e*MX))%*%(X) -
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^4*ivectomat(X*Me+e*MX))%*%(X) -
    t(Pvec^2*e*dM)%*%(M*recD2*recD3blockrep*(onesN%x%t(dM))*M*ivectomat(X*Me+e*MX))%*%(X) +
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^3*ivectomat(X*Me+e*MX))%*%(dM*X)
  A43 <- t(Pvec^2*e*dM)%*%(recD2*recD3blockrep*M^3*ivectomat(X))%*%(X*Me) -
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^4*ivectomat(X))%*%(X*Me) -
    t(Pvec^2*e*dM)%*%(recD2*recD3blockrep*M^2*ivectomat(X*dM))%*%(X*Me) +
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^3*ivectomat(X*dM))%*%(Me*X)
  A44 <- t(dM*X^2)%*%((Poff^2*recD2)%*% (Me*e)) - t((Poff^2*recD2*M) %*% (e))%*% (X^2*Me)
  
  ## A5
  A51 <- t(Pvec^2*e*dM*MX)%*%(recD3blockrep*ivectomat(e)*M)%*%(X) -
    t(Pvec^2*e*MX)%*%(recD3blockrep*M^2*(onesN%x%t(X))*ivectomat(e))%*%(onesN)
  A53 <- t(Pvec^2*e*dM)%*%(recD2*recD3blockrep*M^3*ivectomat(e))%*%(X*MX) -
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^4*ivectomat(e))%*%(X*MX) -
    t(Pvec^2*e*dM)%*%(recD2*recD3blockrep*M^2*ivectomat(e*dM))%*%(X*MX) +
    t(Pvec^2*e)%*%(recD2*recD3blockrep*M^3*ivectomat(e*dM))%*%(MX*X)
  A54 <- t(dM*X*e)%*%((Poff^2*recD2)%*% (MX*e)) - t((Poff^2*recD2*M) %*% (e))%*% (X*e*MX)
  
  
  ret <- (A11 - A12 - A13 + A14 + A15) + 2*(A21 - A22 - A23 +A24 +A25) +
    (A31 - A32 - A33 +A34 +A35) -(A41 +A42 +A43 +A44)-(A51 +A53+A54)
  
  ret
}

L3Ovar_iloop_nocov <- function(X,e,P) {
  n <- length(X)
  M <- diag(n) - P
  
  # Stuff that can be calculated once
  dM <- matrix(diag(M),ncol=1) # force column vector
  D2 <- dM %*% t(dM) - M*M
  onesN <- matrix(rep(1,n),ncol=1)
  recD2 <- 1/D2; diag(recD2) <- 0
  Me <- M%*%e
  MX <- M%*%X
  Poff <- P - diag(diag(P))
  
  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0,n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[,i],ncol=1) #force column vector
    Pi <- matrix(P[,i],ncol=1) #force column vector
    D3i <- M[i,i]*D2-(dM%*%t(Mi)^2+Mi^2%*%t(dM)-2*M*(Mi%*%t(Mi)))
    Di <- matrix(D2[,i],ncol=1)
    D2D3i <- D2/D3i; D2D3i[i,] <- 0; D2D3i[,i] <- 0; diag(D2D3i) <- 0
    recD3i <- 1/D3i; recD3i[i,] <- 0; recD3i[,i] <- 0; diag(recD3i) <- 0
    recD2i <- matrix(1/D2[,i],ncol=1);recD2i[i] <- 0
    Poffi <- Poff[,i]
    tMie <- t(Mi)%*%e; tMiX <- t(Mi)%*%X
    Pie <- Pi*e; PiX <- Pi*X
    MXe <- MX*e; Mie <- Mi*e; MiX <- Mi*X
    MrecD3iPiMie <- (M*recD3i)%*%(Pie*Mi)
    recD3iMMiXMX <- (recD3i*M)%*%(MiX*MX)
    recD3iMMiXMe <- (recD3i*M)%*%(MiX*Me)
    recD3iPiedM <- recD3i%*%(Pie*dM)
    Pi2eMirecD2irecD3iMM <- t(Pi^2*e*Mi*recD2i)%*%(recD3i*M*M)
    Pi2edMMirecD2irecD3i <- t(Pi^2*e*dM*Mi*recD2i)%*%(recD3i)
    PiXD2D3i <- t(PiX)%*%D2D3i
    
    
    A11i <- PiXD2D3i%*%(PiX)*(tMie)
    A12i <- t((Me)*PiX*Mi)%*%recD3i%*%(Poffi*X*dM)-
      t(PiX*Mi)%*%(recD3i*M)%*%(Poffi*X*(Me))
    A13i <- t(dM*PiX)%*%((onesN%x%t(Me))*recD3i)%*%(Pi*MiX)
    A14i <- t((Me)*PiX)%*%(M*recD3i)%*%(Pi*MiX)
    A15i <- (tMie)*(t(Poffi^2*recD2i)%*%(dM*X^2))-
      (t(Poffi^2*Mi*recD2i)%*%(Me*X^2))
    
    A21i <- PiXD2D3i%*%(Pie)*(tMiX)
    A22i <- t((MX)*PiX*Mi)%*%recD3iPiedM-
      t(PiX*Mi)%*%(recD3i*M)%*%(Pie*(MX))
    A23i <- t(dM*PiX)%*%((onesN%x%t(MX))*recD3i)%*%(Pie*Mi)
    A24i <- t((MX)*PiX)%*%MrecD3iPiMie
    A25i <- (tMiX)*(t(Poffi^2*recD2i)%*%(dM*X*e))-
      (t(Poffi^2*Mi*recD2i)%*%(MXe*X))
    
    A31i <- t(Pie)%*%D2D3i%*%(Pie)*(tMiX)
    A32i <- t((MX)*Pie*Mi)%*%recD3iPiedM-
      t(Pie*Mi)%*%(recD3i*M)%*%(Pie*(MX))
    A33i <- t(dM*Pie)%*%((onesN%x%t(MX))*recD3i)%*%(Pie*Mi)
    A34i <- t((MX)*Pie)%*%MrecD3iPiMie
    A35i <- (tMiX)*(t(Poffi^2*recD2i)%*%(dM*e*e))-
      (t(Poffi^2*Mi*recD2i)%*%(MXe*e))
    
    A41i <- t(Pi^2*e*dM*recD2i*Me)%*%((onesN%x%t(Di))*recD3i)%*%(MiX) -
      t(Pi^2*e*recD2i*Me)%*%((onesN%x%t(Di))*recD3i*M*(onesN%x%t(X)))%*%(Mi)
    A42i <- t(Pi^2*e*dM*recD2i)%*%(recD3i*M)%*%(Mi^2*X) -
      Pi2eMirecD2irecD3iMM%*%(MiX) -
      t(Pi^2*e*dM*Mi*recD2i)%*%((onesN%x%t(dM))*recD3i)%*%(MiX) +
      t(Pi^2*e*Mi^2*recD2i)%*%(recD3i*M)%*%(dM*X)
    A43i <- Pi2edMMirecD2irecD3i%*%(Mi^2*X*Me) -
      t(Pi^2*e*Mi^2*recD2i)%*%recD3iMMiXMe -
      M[i,i]*t(Pi^2*e*dM*recD2i)%*%recD3iMMiXMe +
      M[i,i]*Pi2eMirecD2irecD3iMM%*%(Me*X)
    A44i <- X[i]*M[i,i]*(t(Poffi^2*recD2i)%*%(Me*e)) -
      X[i]*(tMie)*(t(Poffi^2*recD2i)%*%(Mie))
    
    A51i <- t(Pi^2*e*dM*recD2i*MX)%*%((onesN%x%t(Di))*recD3i)%*%(Mi*X) -
      t(Pi^2*e*recD2i*MX)%*%((onesN%x%t(Di))*recD3i*M*(onesN%x%t(X)))%*%(Mi)
    A53i <- Pi2edMMirecD2irecD3i%*%(Mi^2*X*MX) -
      t(Pi^2*e*Mi^2*recD2i)%*%recD3iMMiXMX -
      M[i,i]*t(Pi^2*e*dM*recD2i)%*%recD3iMMiXMX +
      M[i,i]*Pi2eMirecD2irecD3iMM%*%(MX*X)
    A54i <- X[i]*M[i,i]*(t(Poffi^2*recD2i)%*%(MXe)) -
      X[i]*(tMiX)*(t(Poffi^2*recD2i)%*%(Mie))
    
    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i*(M[i,]%*%e) + A43i + A44i
    A5vec[i] <- A51i + A42i*(M[i,]%*%X) + A53i + A54i
    
    if (i%%10 == 0) cat(i/n," ")
  }
  
  sum(A1vec*e + 2*A2vec*e +A3vec*X -A4vec*X-A5vec*e)
}

