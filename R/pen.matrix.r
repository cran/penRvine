pen.matrix <- function(penden.env) {
#if(get("pen",penden.env)==3) {
#  bspl.base1<-create.bspline.basis(rangeval=c(0,1),nbasis=dim(get("tilde.Psi.d",penden.env)[,,1])[2],norder=get("q",penden.env),basisvalues=get#("tilde.Psi.d", penden.env)[, , 1],breaks=get("knots1",penden.env))
#  L1<-bsplinepen(bspl.base1,get("m",penden.env))
#  bspl.base2<-create.bspline.basis(rangeval=c(0,1),nbasis=dim(get("tilde.Psi.d",penden.env)[,,2])[2],norder=get("q",penden.env),basisvalues=get("tilde.Psi.d", penden.env)[, , 2],breaks=get("knots2",penden.env))
#  L2<-bsplinepen(bspl.base2,get("m",penden.env))
#}
if(get("pen",penden.env)==2) {
  k.dim <- get("ddb",penden.env)
  k <- k.dim-1
  
  c2 <- factorial(k+1)/factorial(k-2)
  A <- matrix(0,k.dim-2,k.dim)
  diag(A) <- 1
  diag(A[,-1]) <- -2
  diag(A[,-c(1,2)]) <- 1
  
  A.hat <- matrix(NA,k.dim-2,k.dim-2)
  for(i in 0:(k-2)) {
    i.z <- i+1
    for(j in 0:(k-2)) {
      j.z <- j+1
      A.hat[i.z,j.z] <- choose(k-2,j)*choose(k-2,i)*beta(i+j+1,2*k-i-j-3)
    }
  }  
  A2 <- matrix(NA,k.dim,k.dim)
  for(i in 0:k) {
    i.z <- i+1
    for(j in 0:k) {
      j.z <- j+1
      A2[i.z,j.z] <- choose(k,j)*choose(k,i)*beta(i+j+1,2*k-i-j+1)
    }
  }
  assign("DDD.sum",mat<-(kronecker(c2^2*(t(A)%*%A.hat%*%A),A2)+kronecker(A2,c2^2*(t(A)%*%A.hat%*%A))),penden.env)
  assign("eigen.pen.mat",help2 <- eigen(mat),penden.env)
  assign("index.eigen.pen.mat",index <- which(help2$values>1e-08),penden.env)
  assign("Utilde.eigen.pen.mat",help2$vectors[,index],penden.env)
}

if(get("pen",penden.env)==1) {
  m <- get("m",penden.env)
  K <- get("ddb",penden.env)
  q <- get("q",penden.env)
  #if(get("base",penden.env)=="B-spline"&get("q",penden.env)==1) K <- K-1
  if(m==1) {
    L <- diag(1,K)
    L.1 <- diag(-1,K,K-1)
    L.2 <- matrix(0,K,1)
    L1 <- cbind(L.2,L.1)
    L <- L+L1
    L <- L[1:(K-1),]
  }
  if(m==2) {
    L <- diag(1,K,K)
    L1 <- diag(-2,K,(K-1))
    L2 <- diag(1,K,(K-2))
    L.1 <- matrix(0,K,1)
    L1 <- cbind(L.1,L1)
    L2 <- cbind(L.1,L.1,L2)
    L <- L+L1+L2
    L <- L[1:(K-2),]
  }
  if(m==3) {
    L <- diag(1,(K-3),(K-3))
    L.help <- matrix(0,(K-3),1)
    L1 <- diag(-3,(K-3),(K-3))
    M1 <- cbind(L,L.help,L.help,L.help)
    M2 <- cbind(L.help,L1,L.help,L.help)
    M3 <- cbind(L.help,L.help,-L1,L.help)
    M4 <- cbind(L.help,L.help,L.help,-L)
    L <- (M1+M2+M3+M4)
  }
  if(m==4) {
    L <- diag(1,(K-4),(K-4))
    L.help <- matrix(0,(K-4),1)
    L1 <- diag(-4,(K-4),(K-4))
    L2 <- diag(6,(K-4),(K-4))
    M1 <- cbind(L,L.help,L.help,L.help,L.help)
    M2 <- cbind(L.help,L1,L.help,L.help,L.help)
    M3 <- cbind(L.help,L.help,L2,L.help,L.help)
    M4 <- cbind(L.help,L.help,L.help,L1,L.help)
    M5 <- cbind(L.help,L.help,L.help,L.help,L)
    L <- (M1+M2+M3+M4+M5)
}  
}
#if (q > 2) {
#  dd.A<-dim(L)[1]
#  {    for (j in 3:q)
#      {
#        {
#          L <-  L[1:(dd.A-j),] - L[2:(dd.A-j+1),]
#        }
#      }
#  }
#}
  Index.basis.D <- get("Index.basis.D",penden.env)
  #if(get("base",penden.env)=="B-spline"&get("knots.place",penden.env)=="equi") {
  #  C1 <- get("C",penden.env)
  #  ADA <- t(C1) %*% crossprod(L) %*% (C1)
  #  AIA <- crossprod(C1)
  #  DD <- get("DD",penden.env)
  #  DDD3 <- array(NA, c(DD,DD, 2))
  #  for (k in 1:2)
  #    {
  #      l.ind <- (1:2)[-k]
  #      i <- 1:DD
  #      j <- 1:DD
  #      DDD3[i,j,k] <- ADA[Index.basis.D[i,k], Index.basis.D[j,k]]
  #      for (l in l.ind)   DDD3[i,j,k] <- DDD3[i,j,k] * AIA[Index.basis.D[i,l], Index.basis.D[j,l]]
  #    }
  #  mat <- DDD3[,,1]+DDD3[,,2]
  #}
  if(get("base",penden.env)=="B-spline") {
    pen<-get("pen",penden.env)
    ADA<-AIA<-list()
    C1 <- get("C1",penden.env)
    if(pen==1|pen==2) ADA[[1]] <- t(C1) %*% crossprod(L) %*% (C1)
    if(pen==3) ADA[[1]] <- t(C1) %*% crossprod(L1) %*% (C1)
    AIA[[1]] <- crossprod(C1)
    C2<- get("C2",penden.env)
    if(pen==1|pen==2) ADA[[2]] <- t(C2) %*% crossprod(L) %*% (C2)
    if(pen==3) ADA[[2]] <- t(C2) %*% crossprod(L2) %*% (C2)
    AIA[[2]] <- crossprod(C2)
    DD <- get("DD",penden.env)
    DDD3 <- array(NA, c(DD,DD, 2))
    for (k in 1:2)
      {
        l.ind <- (1:2)[-k]
        i <- 1:DD
        j <- 1:DD
        DDD3[i,j,k] <- ADA[[k]][Index.basis.D[i,k], Index.basis.D[j,k]]
        for (l in l.ind)   DDD3[i,j,k] <- DDD3[i,j,k] * AIA[[l]][Index.basis.D[i,l], Index.basis.D[j,l]]
      }
    mat <- DDD3[,,1]+DDD3[,,2]
  }

  if(get("base",penden.env)=="Bernstein"&get("pen",penden.env)==1)   mat <- kronecker(diag(1,K),crossprod(L))+kronecker(crossprod(L),diag(1,K))
  assign("eigen.pen.mat",help2 <- eigen(mat),penden.env)
  assign("index.eigen.pen.mat",index <- which(help2$values>1e-08),penden.env)
  assign("Utilde.eigen.pen.mat",help2$vectors[,index],penden.env)
  assign("DDD.sum",mat,penden.env)
1
}
 
