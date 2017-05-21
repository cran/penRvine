order.vine <- function(help.env,test.ind=FALSE) {
  
  RVM<-get("RVM",help.env)
  len <- get("p",help.env)
  no.pairs <- choose(len,2)
  order.stat <- get("order.stat",help.env)
  pairs.new <- matrix(NA,no.pairs,2)
  count <- 1
  ddb<-get("ddb",help.env)
  if(is.null(RVM)) {
    for(i in 1:(len-1)) {
      for(j in (i+1):len) {
        pairs.new[count,] <- c(i,j)
        count <- count+1
      }
    }
  }
  else {
    pairs.new<-pairs.fit<-get("cops",help.env)[[2]]
    no.pairs<-dim(pairs.new)[1]
  }
  l.search<-get("l.search",help.env)
  id<-get("id",help.env)
  if(get("selec",help.env)=="cAIC") {
  mcoptions <- list(preschedule=FALSE)
  #print(paste("level=2 ,number of pairs=",no.pairs,sep=""))
  h1 <- foreach(i=1:no.pairs,.options.multicore=mcoptions) %do% {
      if(test.ind) {
        test.UU<-independ.test(get("U",help.env)[,pairs.new[i,]])
        if(test.UU<=0.05) {
          model.l<-paircopula(data=get("U",help.env)[,pairs.new[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=get("q",help.env),max.iter=get("max.iter",help.env),fix.lambda=get("fix.lambda",help.env))
          val <- c(v=get("ck.val",model.l),log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.est=get("f.hat.val",model.l),indep=FALSE)
        }
        if(test.UU>0.05) {
          val<- c(v=get("ck.val.start",get("indep.model",help.env)),log.like = 0, AIC=0,cAIC=-0.00000000001,lambda=0,f.est=get("f.hat.val.start",get("indep.model",help.env)),indep=TRUE)
        }
      }
      else {
        if(l.search) model.l<-lam.search(data=get("U",help.env)[,pairs.new[i,]],K=get("K",help.env),lam=get("lam.vec",help.env),m=get("m",help.env),max.iter=get("max.iter",help.env),q=get("q",help.env),base=get("base",help.env),fix.lambda=get("fix.lambda",help.env),pen=get("pen",help.env),id=id)
        if(!l.search) model.l<-paircopula(data=get("U",help.env)[,pairs.new[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=get("q",help.env),max.iter=get("max.iter",help.env),fix.lambda=get("fix.lambda",help.env))
        val <- list(v=get("ck.val",model.l),log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.hat=get("f.hat.val",model.l),indep=FALSE,knots1=get("knots1",model.l),knots2=get("knots2",model.l),Index.basis.D=get("Index.basis.D",model.l),ddb=get("ddb",model.l),int.base1=get("int.base1",model.l),int.base2=get("int.base2",model.l),K=get("K",model.l))
     }
  val
  }

  h <- foreach(i=1:no.pairs,.combine=rbind) %do% {
     if(order.stat=="cAIC") obj<-c(pairs.new[i,],h1[[i]]$cAIC)
     if(order.stat=="AIC") obj<-c(pairs.new[i,],h1[[i]]$AIC)
     obj
  }

  colnames(h) <- c("i","j","cAIC")
  if(is.null(RVM)) {
    mat <- matrix(NA,len,len)
    diag(mat) <- rep(0,len)
    for(i in 1:(len-1)) {
      for(j in (i+1):len) {
        mat[i,j] <- mat[j,i] <- h[which(h[,1]==i & h[,2]==j),3]
      }
    }
  }
  assign("pairs.new",pairs.new,help.env)
  assign("fit.level2",h1,help.env)
  assign("fit.results",h,help.env)
  if(is.null(RVM)) {
    obj <- minimum.spanning.tree(graph.adjacency(mat,diag=FALSE,mode="lower",weighted=TRUE),algorithm="prim")
    assign("min.sp.tree.2",obj,help.env)
    pairs.fit <- get.edgelist(obj, names=TRUE)
    rm(obj)
    pairs.fit <- pairs.fit[order(pairs.fit[,1]),]
    VineMatrix <- matrix(NA,dim(pairs.fit)[1],len+1)
    VineMatrix[,c(1,2)] <- pairs.fit
  }
  else  {
    VineMatrix <- matrix(NA,dim(pairs.fit)[1],len+1)
    VineMatrix[,c(1,2)] <- pairs.fit
  }
  numbers <- c()
  if(dim(pairs.fit)[1]<10) {
    assign("dez",10,help.env)
    numbers[1:dim(pairs.fit)[1]] <- as.numeric(paste("2.",seq(1,dim(pairs.fit)[1]),sep=""))
  }
  if(dim(pairs.fit)[1]>9&dim(pairs.fit)[1]<100) {
    assign("dez",100,help.env)
    numbers[1:9] <- as.numeric(paste("2.0",seq(1,9),sep=""))
    numbers[10:dim(pairs.fit)[1]]  <- as.numeric(paste("2.",seq(10,dim(pairs.fit)[1]),sep=""))
  }
  VineMatrix[,(len+1)] <- numbers
  assign("VineMatrix",VineMatrix,help.env)
  assign("pairs.fit",pairs.fit,help.env)
}


if(get("selec",help.env)=="ken.tau") {
  mcoptions <- list(preschedule=FALSE)
  h <- foreach(i=1:no.pairs,.combine=rbind,.options.multicore=mcoptions) %do% {
    c(pairs.new[i,],cor(x=get("U",help.env)[,pairs.new[i,1]], y = get("U",help.env)[,pairs.new[i,2]], use = "everything", method = c("kendall")),weight=1-(cor(x=get("U",help.env)[,pairs.new[i,1]], y = get("U",help.env)[,pairs.new[i,2]], use = "everything", method = c("kendall"))))
  }
  if(is.null(RVM)) {
    colnames(h) <- c("i","j","ken.tau","weights")
    mat <-matrix(NA,len,len)
    diag(mat) <- rep(0,len)
    for(i in 1:(len-1)) {
      for(j in (i+1):len) {
        mat[i,j] <- mat[j,i] <- h[which(h[,1]==i & h[,2]==j),3]
      }
    }
     assign("pairs.new",pairs.new,help.env)
     graph.obj<-graph.adjacency(mat,diag=FALSE,mode="upper",weighted=TRUE)
     ed.list<-as_edgelist(graph.obj)
     ord<-c()
     for(jj in 1:dim(ed.list)[1]) ord[jj]<- which(ed.list[jj,1]==pairs.new[,1]&ed.list[jj,2]==pairs.new[,2])
     obj <- minimum.spanning.tree(graph.obj,algorithm="prim",weights=h[ord,4])
     assign("min.sp.tree.2",obj,help.env)
     pairs.fit <- get.edgelist(obj, names=TRUE)
     VineMatrix <- matrix(NA,dim(pairs.fit)[1],len+1)
     VineMatrix[,c(1,2)] <- pairs.fit
  }
  else {
    VineMatrix <- matrix(NA,dim(pairs.fit)[1],len+1)
    VineMatrix[,c(1,2)] <- pairs.fit
  }
  numbers <- c()
  if(dim(pairs.fit)[1]<10) {
    assign("dez",10,help.env)
    numbers[1:dim(pairs.fit)[1]] <- as.numeric(paste("2.",seq(1,dim(pairs.fit)[1]),sep=""))
  }
  if(dim(pairs.fit)[1]>9&dim(pairs.fit)[1]<100) {
    assign("dez",100,help.env)
    numbers[1:9] <- as.numeric(paste("2.0",seq(1,9),sep=""))
    numbers[10:dim(pairs.fit)[1]]  <- as.numeric(paste("2.",seq(10,dim(pairs.fit)[1]),sep=""))
  }
  VineMatrix[,(len+1)] <- numbers
  assign("VineMatrix",VineMatrix,help.env)
  fit.level2 <- foreach(i=1:dim(pairs.fit)[1],.options.multicore=mcoptions) %do% {
      if(test.ind) {
        test.UU<-independ.test(get("U",help.env)[,pairs.new[i,]])
        if(test.UU<=0.05) {
          model.l<-paircopula(data=get("U",help.env)[,pairs.fit[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=get("q",help.env),fix.lambda=get("fix.lambda",help.env))
          val <- c(v=get("ck.val",model.l),log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.est=get("f.hat.val",model.l),indep=FALSE)
        }
        if(test.UU>0.05) {
          val<- c(v=get("ck.val.start",get("indep.model",help.env)),log.like = 0, AIC=0,cAIC=-0.00000000001,lambda=0,f.est=get("f.hat.val.start",get("indep.model",help.env)),indep=TRUE)
        }
      }
      else {
        if(l.search) model.l<-lam.search(data=get("U",help.env)[,pairs.fit[i,]],K=get("K",help.env),lam=get("lam.vec",help.env),m=get("m",help.env),max.iter=get("max.iter",help.env),q=get("q",help.env),base=get("base",help.env),fix.lambda=get("fix.lambda",help.env),pen=get("pen",help.env),id=id)
        if(!l.search) model.l<-paircopula(data=get("U",help.env)[,pairs.fit[i,]],K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=get("q",help.env),fix.lambda=get("fix.lambda",help.env))
        val <- list(v=get("ck.val",model.l),log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.hat=get("f.hat.val",model.l),indep=FALSE,knots1=get("knots1",model.l),knots2=get("knots2",model.l),Index.basis.D=get("Index.basis.D",model.l),ddb=get("ddb",model.l),int.base1=get("int.base1",model.l),int.base2=get("int.base2",model.l),K=get("K",model.l))
     }
  val
  }
  assign("fit.level2",fit.level2,help.env)
  assign("pairs.fit",pairs.fit,help.env)
  assign("pairs.new",pairs.new,help.env)
}
}
