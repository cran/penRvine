cal.vine <- function(obj,val,cores) {
  Dvine<-obj$vine
  registerDoParallel(cores=cores)
  K <- obj$K
  S<-obj$S
  N <- obj$N
  base <- obj$base
  help.env <- obj$env
  q<-obj$q
  level <-  2
  Tree.l <- Dvine[[level]]
  if(base=="Bernstein") ddb <- K
  #if(base=="B-spline" & q==2) ddb <- K
  #if(base=="B-spline" & q==3) ddb <- K+1
  #ddb<-get("ddb",help.env)
  #Index.basis.D<-get("Index.basis.D",help.env)
  #Index.basis.D <- matrix(NA,(ddb)^2,2)
  #Index.basis.D[,1] <- rep(seq(1,(ddb)),(ddb))
  #Index.basis.D[,2] <- sort(Index.basis.D[,1])

log.like <- AIC <- 0
  Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list, .multicombine=TRUE) %do%   {
    index.j1 <- Tree.l[[i]]$j1
    index.j2 <- Tree.l[[i]]$j2
    U.hat <- val[,c(index.j1,index.j2)]
    U <- eval.paircopula(U.hat,K=Tree.l[[i]]$K,int=FALSE,Index.basis.D=Tree.l[[i]]$Index.basis.D,ck.val=obj$vine[[level]][[i]]$v,base=base,q=q,kn1=Tree.l[[i]]$knots1,kn2=Tree.l[[i]]$knots2,ddb=Tree.l[[i]]$ddb,int.base1=Tree.l[[i]]$int.base1,int.base2=Tree.l[[i]]$int.base2)
    ll <- sum(sapply(U,log))
    list(j1=index.j1,j2=index.j2,D=NULL,v=obj$vine[[level]][[i]]$v,U.hat=U.hat,U=U,log.like=ll)#,cop=obj$vine[[level]][[i]]$cop) 
  }
  for(j in 1:length(Tree.l.temp)) {
    log.like <- log.like+Tree.l.temp[[j]]$log.like
  }

  Dvine[[level]] <- Tree.l.temp
  while( length(Dvine[[level]])>2) # loop over the levels of the Dvine. We stop if the tree in the Dvine has just one knot
  {
    level <-  level+1
    Tree.l <- Dvine[[level]] # current tree in vine
    Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
    Tree.l.temp <- foreach(i=1:length(Tree.l),.combine=list, .multicombine=TRUE) %do% {
       U.hat <- c()
       for(j in 1:2)
          {
              {
                index.ancestor <- Tree.l[[i]]$before[j]
              }
            ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}  # for j=1 and {j2,D} for j=2, respectively.
            if(Tree.l[[i]]$j1%in%c(ancestor.knot$j1,ancestor.knot$j2)) val.j<-Tree.l[[i]]$j1
            if(Tree.l[[i]]$j2%in%c(ancestor.knot$j1,ancestor.knot$j2)) val.j<-Tree.l[[i]]$j2
            diff.help <- c("u1","u2")[!(c(ancestor.knot$j1,ancestor.knot$j2)%in%val.j)]      
            U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=obj$vine[[level]][[i]]$K,diff=diff.help,Index.basis.D=obj$vine[[level]][[i]]$Index.basis.D,base=base,q=q,env=help.env,ddb=obj$vine[[level]][[i]]$ddb,int.base1=obj$vine[[level]][[i]]$int.base1,int.base2=obj$vine[[level]][[i]]$int.base2,kn1=obj$vine[[level]][[i]]$knots1,kn2=obj$vine[[level]][[i]]$knots2))
          } 
        U.hat[U.hat<0]<-0
        U.hat[U.hat>1]<-1
        U <- eval.paircopula(val=U.hat,K=obj$vine[[level]][[i]]$K,int=FALSE,Index.basis.D=obj$vine[[level]][[i]]$Index.basis.D,ck.val=obj$vine[[level]][[i]]$v,base=base,q=q,kn1=obj$vine[[level]][[i]]$knots1,kn2=obj$vine[[level]][[i]]$knots2,ddb=obj$vine[[level]][[i]]$ddb,int.base1=obj$vine[[level]][[i]]$int.base1,int.base2=obj$vine[[level]][[i]]$int.base2)
        ll <- sum(sapply(U,log))
        list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U.hat=U.hat,U=U,v=obj$vine[[level]][[i]]$v,log.like=ll)#,cop=obj$vine[[level]][[i]]$cop)
    }
    
    if(length(Dvine[[level]])==1)  {
      Tree.l.temp <- list(Tree.l.temp)
    }
    Dvine[[level]] <- Tree.l.temp
    
    for(i in 1:length(Tree.l.temp)) {
      log.like <- log.like+Tree.l.temp[[i]]$log.like
    }
  }
  level <-  level+1
  Tree.l <- Dvine[[level]] # current tree in vine
  Tree.l1 <- Dvine[[level-1]] # previous tree in vine, one level up
  i <- 1
  U.hat <- c()
  index <- list(c(Tree.l[[i]]$j1,Tree.l[[i]]$D),c(Tree.l[[i]]$j2,Tree.l[[i]]$D))
  for(j in 1:2)
    {
      index.ancestor <- c()
      for (m in 1:length(Tree.l1))
        {
          index.ancestor <- Tree.l[[i]]$before[j]
        }
      ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D}
      if(Tree.l[[i]]$j1%in%c(ancestor.knot$j1,ancestor.knot$j2)) val.j<-Tree.l[[i]]$j1
      if(Tree.l[[i]]$j2%in%c(ancestor.knot$j1,ancestor.knot$j2)) val.j<-Tree.l[[i]]$j2
      diff.help <- c("u1","u2")[!(c(ancestor.knot$j1,ancestor.knot$j2)%in%val.j)]
      U.hat <-cbind(U.hat,cond.cop(data=ancestor.knot$U.hat,coef=ancestor.knot$v,K=obj$vine[[level]][[i]]$K,diff=diff.help,Index.basis.D=obj$vine[[level]][[i]]$Index.basis.D,base=base,q=q,env=help.env,ddb=obj$vine[[level]][[i]]$ddb,int.base1=obj$vine[[level]][[i]]$int.base1,int.base2=obj$vine[[level]][[i]]$int.base2,kn1=obj$vine[[level]][[i]]$knots1,kn2=obj$vine[[level]][[i]]$knots2))
     }
  U.hat[U.hat<0]<-0
  U.hat[U.hat>1]<-1

  U <- eval.paircopula(val=U.hat,K=obj$vine[[level]][[i]]$K,int=FALSE,Index.basis.D=obj$vine[[level]][[i]]$Index.basis.D,ck.val=obj$vine[[level]][[i]]$v,base,q=q,kn1=obj$vine[[level]][[i]]$knots1,kn2=obj$vine[[level]][[i]]$knots2,ddb=obj$vine[[level]][[i]]$ddb,int.base1=obj$vine[[level]][[i]]$int.base1,int.base2=obj$vine[[level]][[i]]$int.base2)
 
  ll <- sum(sapply(U,log))
  Tree.l.temp <- list(j1=Tree.l[[i]]$j1,j2=Tree.l[[i]]$j2,D=Tree.l[[i]]$D,U=U,v=obj$vine[[level]][[i]]$v,log.like=ll)#,cop=obj$vine[[level]][[i]]$cop)
  Dvine[[level]][[1]] <- Tree.l.temp
  log.like <- log.like+ll
  U.result <- rep(1,dim(val)[1])
  for(i in 2:length(S))  {
    for(j in 1:length(Dvine[[i]]))  {
      U.result <- U.result*Dvine[[i]][[j]]$U
    }
  }
  return(list(cal=U.result,log.like=log.like,Dvine=Dvine))
}
