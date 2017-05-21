order.vine.level <- function(tree,help.env) {
  l.search<-get("l.search",help.env)
  id<-get("id",help.env)
  RVM<-get("RVM",help.env)
  len <- length(tree)
  base <- get("base",help.env)
  K <- get("K",help.env)
  order.stat <- get("order.stat",help.env)
  pairs.old <- get("pairs.fit",help.env)
  q <- get("q",help.env)
  base <- get("base",help.env)
  val <- unique(c(pairs.old))
  count <- 1
  pairs.new <- matrix(NA,1,2)
    for(i in 1:length(val)) {
      val.temp <- c()
      for(j in 1:dim(pairs.old)[1]) {
        if(any(pairs.old[j,]%in%val[i])) val.temp <- c(val.temp,j)
      }
      if(length(val.temp)==2) {
        if(count==1) {
          pairs.new[count,] <- val.temp
          count <- count+1
        }
        else {
          pairs.new <- rbind(pairs.new,val.temp)
          count <- count+1
        }
      }
      if(length(val.temp)>2) {
        for(j in 1:(length(val.temp)-1)) {
          for(k in (j+1):length(val.temp)) {
	    if(count==1) {
              pairs.new[count,] <- c(val.temp[j],val.temp[k])
              count <- count+1
            }
            else {
              pairs.new <- rbind(pairs.new,c(val.temp[j],val.temp[k]))
              count <- count+1
            }
          }
        }
      }
    }
  no.pairs <- dim(pairs.new)[1]
  mcoptions <- list(preschedule=FALSE)
  if(!is.null(RVM)) {
    cops<-get("cops",help.env)[[get("level",help.env)]]
    h.help <- foreach(i=1:no.pairs,.options.multicore=mcoptions,.combine=rbind) %do% {
      UU <- c()
      help.j1 <- c(tree[[pairs.new[i,1]]]$j1,tree[[pairs.new[i,1]]]$j2,tree[[pairs.new[i,1]]]$D)
      help.j2 <- c(tree[[pairs.new[i,2]]]$j1,tree[[pairs.new[i,2]]]$j2,tree[[pairs.new[i,2]]]$D)
      j1 <- help.j1[!help.j1%in%help.j2]
      j2 <- help.j2[!help.j2%in%help.j1]
      D.help <- c(tree[[pairs.new[i,1]]]$D,tree[[pairs.new[i,2]]]$D)
      D <- sort(unique(c(help.j1[help.j1%in%help.j2],D.help)))
      c(j1,j2,D)
    }
    ind<-rep(NA,dim(cops)[1])
    ind2<-c()
    h.help2<-cbind(h.help[,2],h.help[,1],h.help[,-c(1,2)])
    for(k in 1:dim(h.help)[1]) {
      for(ll in 1:dim(cops)[1]) {
        if(identical(h.help[k,],cops[ll,])) ind[ll]<-k
        if(identical(h.help2[k,],cops[ll,])) { 
          ind[ll]<-k
          ind2<-c(ind2,k)
        }
      }
    }
    if(!is.null(ind2)) pairs.new[ind2,]<-c(pairs.new[ind2,2],pairs.new[ind2,1])
    pairs.new<-pairs.new[ind,]
  }

  no.pairs <- dim(pairs.new)[1]

if(get("selec",help.env)=="cAIC") {
  h1 <- foreach(i=1:no.pairs,.combine=c,.options.multicore=mcoptions) %do%{
      UU <- c()
      help.j1 <- c(tree[[pairs.new[i,1]]]$j1,tree[[pairs.new[i,1]]]$j2,tree[[pairs.new[i,1]]]$D)
      help.j2 <- c(tree[[pairs.new[i,2]]]$j1,tree[[pairs.new[i,2]]]$j2,tree[[pairs.new[i,2]]]$D)
      j1 <- help.j1[!help.j1%in%help.j2]
      j2 <- help.j2[!help.j2%in%help.j1]
      D.help <- c(tree[[pairs.new[i,1]]]$D,tree[[pairs.new[i,2]]]$D)
      D <- sort(unique(c(help.j1[help.j1%in%help.j2],D.help)))
      index <- list(c(j1,D),c(j2,D)) # list of involved indices in current knot, i.e {j1,D} and {j2,D}
      for(j in 1:2) # Loop over both index sets, i.e. {j1,D} and {j2,D}
        {
          indexi <- c(tree[[pairs.new[i,j]]]$j1,tree[[pairs.new[i,j]]]$j2,tree[[pairs.new[i,j]]]$D)
          index.ancestor <- c()#index.indep<- c()
          for (ml in 1:length(tree))
            {
              index.ancestor <- c(index.ancestor, all(indexi==(c(tree[[ml]]$j1,tree[[ml]]$j2,tree[[ml]]$D))))
            }
          ancestor.knot <- tree[index.ancestor][[1]] 
          if(j==1) {
            diff.help <- c("u1","u2")[!(c(tree[[pairs.new[i,1]]]$j1,tree[[pairs.new[i,1]]]$j2)%in%j1)]
          }
          if(j==2) {
            diff.help <- c("u1","u2")[!(c(tree[[pairs.new[i,2]]]$j1,tree[[pairs.new[i,2]]]$j2)%in%j2)]
          }
          UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=ancestor.knot$K,diff=diff.help,ancestor.knot$Index.basis.D,base=base,q=q,env=help.env,kn1=ancestor.knot$knots1,kn2=ancestor.knot$knots2,int.base1=ancestor.knot$int.base1,int.base2=ancestor.knot$int.base2,ddb=ancestor.knot$ddb))
        }
      if(any(UU>1)) UU[which(UU>1)] <-1
      if(any(UU<0)) UU[which(UU<0)] <-0
          if(l.search) model.l<-lam.search(data=UU,K=get("K",help.env),lam=get("lam.vec",help.env),m=get("m",help.env),max.iter=get("max.iter",help.env),q=get("q",help.env),base=get("base",help.env),fix.lambda=get("fix.lambda",help.env),pen=get("pen",help.env),id=id)
          if(!l.search) model.l <- paircopula(data=UU,K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=q,max.iter=get("max.iter",help.env),fix.lambda=get("fix.lambda",help.env))
          assign("indep",FALSE,model.l)
  model.l    
  }
  if(is.null(RVM)) {
    h <- foreach(i=1:no.pairs,.combine=rbind) %do% {
      c(pairs.new[i,],get(order.stat,h1[[i]]))
    }
    colnames(h) <- c("i","j","log.like")
    mat <- matrix(0,len,len)
    diag(mat) <- rep(0,len)

    for(i in 1:dim(pairs.new)[1]) {
      mat[pairs.new[i,1],pairs.new[i,2]] <- mat[pairs.new[i,2],pairs.new[i,1]] <- h[which(h[,1]==pairs.new[i,1] & h[,2]==pairs.new[i,2]),3]
    }
  }
  assign("pairs.new",pairs.new,help.env)
  assign(paste("fit.level",get("level",help.env),sep=""),h1,help.env)
  #assign("fit.results",h,help.env)
  if(is.null(RVM)) {
    obj <- minimum.spanning.tree(graph.adjacency(mat,diag=FALSE,mode="lower",weighted=TRUE),algorithm="prim")
    pairs.fit <- get.edgelist(obj, names=TRUE)
    pairs.fit <- pairs.fit[order(pairs.fit[,1]),]
  }
  else pairs.fit<-pairs.new
  assign("pairs.fit",pairs.fit,help.env)

  #assign(paste("min.sp.tree.",get("level",help.env),sep=""),obj,help.env)
}
if(get("selec",help.env)=="ken.tau") {
 h1 <- foreach(i=1:no.pairs) %do% {
      UU <- c()
      help.j1 <- c(tree[[pairs.new[i,1]]]$j1,tree[[pairs.new[i,1]]]$j2,tree[[pairs.new[i,1]]]$D)
      help.j2 <- c(tree[[pairs.new[i,2]]]$j1,tree[[pairs.new[i,2]]]$j2,tree[[pairs.new[i,2]]]$D)
      j1 <- help.j1[!help.j1%in%help.j2]
      j2 <- help.j2[!help.j2%in%help.j1]
      D.help <- c(tree[[pairs.new[i,1]]]$D,tree[[pairs.new[i,2]]]$D)
      D <- sort(unique(c(help.j1[help.j1%in%help.j2],D.help)))
      index <- list(c(j1,D),c(j2,D)) # list of involved indices in current knot, i.e {j1,D} and {j2,D}
      for(j in 1:2) # Loop over both index sets, i.e. {j1,D} and {j2,D}
        {
          indexi <- c(tree[[pairs.new[i,j]]]$j1,tree[[pairs.new[i,j]]]$j2,tree[[pairs.new[i,j]]]$D)
          index.ancestor <- c()#index.indep<- c()
          for (ml in 1:length(tree))
            {
              index.ancestor <- c(index.ancestor, all(indexi==(c(tree[[ml]]$j1,tree[[ml]]$j2,tree[[ml]]$D))))
            }
          ancestor.knot <- tree[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D} for j=1 and {j2,D} for j=2, respectively.
          #index.indep<-c(tree[index.ancestor][[1]]$indep,index.indep)
          
          if(j==1) {
            diff.help <- c("u1","u2")[!(c(tree[[pairs.new[i,1]]]$j1,tree[[pairs.new[i,1]]]$j2)%in%j1)]
          }
          if(j==2) {
            diff.help <- c("u1","u2")[!(c(tree[[pairs.new[i,2]]]$j1,tree[[pairs.new[i,2]]]$j2)%in%j2)]
          }
          UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=ancestor.knot$K,diff=diff.help,ancestor.knot$Index.basis.D,base=base,q=q,env=help.env,kn1=ancestor.knot$knots1,kn2=ancestor.knot$knots2,int.base1=ancestor.knot$int.base1,int.base2=ancestor.knot$int.base2,ddb=ancestor.knot$ddb))
        }
      if(any(UU>1)) UU[which(UU>1)] <-1
      if(any(UU<0)) UU[which(UU<0)] <-0
 list(stat=c(pairs=pairs.new[i,],kenpairs.tau=cor(x=UU[,1], y =UU[,2], use = "everything", method = c("kendall")),weights=1-abs(cor(x=UU[,1], y =UU[,2], use = "everything", method = c("kendall")))),UU=UU)
 }
  h <- foreach(i=1:no.pairs,.combine=rbind) %do% {
    h1[[i]]$stat
  }
  if(is.null(RVM)) {
    mat <- matrix(0,len,len)
    diag(mat) <- rep(0,len)
    for(i in 1:dim(pairs.new)[1]) {
      mat[pairs.new[i,1],pairs.new[i,2]] <- mat[pairs.new[i,2],pairs.new[i,1]] <- h[which(h[,1]==pairs.new[i,1] & h[,2]==pairs.new[i,2]),3]
    }
    assign("pairs.new",pairs.new,help.env)
    graph.obj<-graph.adjacency(mat,diag=FALSE,mode="upper",weighted=TRUE)
    ed.list<-as_edgelist(graph.obj)
    ord<-c()
    for(jj in 1:dim(ed.list)[1]) ord[jj]<- which(ed.list[jj,1]==pairs.new[,1]&ed.list[jj,2]==pairs.new[,2])
    obj <- minimum.spanning.tree(graph.obj,algorithm="prim",weights=h[ord,4])
    pairs.fit <- get.edgelist(obj, names=TRUE)
    sel<-c()
    for(i in 1:dim(pairs.fit)[1]) {
      for(kk in 1:dim(pairs.new)[1]) {
        if(all(pairs.new[kk,]==pairs.fit[i,])) sel<-c(sel,kk)
      }  
    }
    h1<-h1[sel]
    assign("pairs.fit",pairs.fit,help.env)
  }
  else {
    pairs.fit<-pairs.new
    assign("pairs.fit",pairs.fit,help.env)
  }
  h1.new <- foreach(i=1:dim(pairs.fit)[1],.options.multicore=mcoptions) %do% {
        if(get("test.ind",help.env)){
          test.UU<-independ.test(h1[[i]]$UU)
          if(test.UU<=0.05) {
            model.l <- paircopula(data=h1[[i]]$UU,K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=q,max.iter=get("max.iter",help.env))
            assign("indep",FALSE,model.l)
          }
          if(test.UU>0.05) {
            model.l<-get("indep.model",help.env)
            assign("indep",TRUE,model.l)
          }
        }
        if(!get("test.ind",help.env)) {
          if(!l.search) model.l <- paircopula(data=h1[[i]]$UU,K=get("K",help.env),lambda=get("lambda",help.env),pen=get("pen",help.env),base=get("base",help.env),m=get("m",help.env),q=q,max.iter=get("max.iter",help.env),fix.lambda=get("fix.lambda",help.env))
          if(l.search) model.l<-lam.search(data=h1[[i]]$UU,K=get("K",help.env),lam=get("lam.vec",help.env),m=get("m",help.env),max.iter=get("max.iter",help.env),q=get("q",help.env),base=get("base",help.env),fix.lambda=get("fix.lambda",help.env),id=id,pen=get("pen",help.env))
          assign("indep",FALSE,model.l)
        }
  model.l    
  }
assign(paste("fit.level",get("level",help.env),sep=""),h1.new,help.env)
}
}
