vine <- function(data,K=8,lambda=100,pen=1,base="B-spline",m=2,cores=NULL,q=2,test.ind=FALSE,test.ind1=FALSE,selec="cAIC",max.iter=51,RVM=NULL,lam.vec=NULL,l.search=FALSE,fix.lambda=FALSE,id=NULL) {

if(is.null(cores)) registerDoParallel()
else registerDoParallel(cores=cores)

g1<-any(data>1)
s0<-any(data<0)
if(g1&s0) {
   print("Any data greater 1 and smaller 0")
   break
}
if(g1&!s0) {
   print("Any data greater 1")
   break
}
if(s0&!g1) {
   print("Any data smaller 0")
   break
}
if(!is.matrix(data)) data <- as.matrix(data)
U <- data
S <- seq(1:(dim(U)[2]))
N <-  dim(U)[1]
k <-  1
vine <- list()
SSi <-  list()

help.env <- new.env()
if(is.null(lam.vec)&l.search) lam.vec<-c(10,seq(50,1550,by=250))
assign("lam.vec",lam.vec,help.env)
assign("fix.lambda",fix.lambda,help.env)
assign("l.search",l.search,help.env)
assign("test.ind",test.ind,help.env)
assign("S",S,help.env)
assign("max.iter",max.iter,help.env)
assign("p",length(S),help.env)
assign("U",U,help.env)
assign("lambda",lambda,help.env)
assign("base",base,help.env)
assign("m",m,help.env)
assign("order.stat","cAIC",help.env)
assign("q",q,help.env)
assign("selec",selec,help.env)
assign("id",id,help.env)
if(!is.null(RVM)) assign("RVM",reRVineMatrix(help.env),help.env)

if(base=="Bernstein") ddb <- K+1 #Anzahl Bernsteinpolynome
if(base=="B-spline") ddb <- K-1+2*(q-1)
if(base=="B-spline") K<-K+q-1
if(base=="Bernstein") int.bernstein.help(help.env)
assign("index.k",matrix(0:K),help.env)

order.vine(help.env,test.ind=test.ind1)

for ( i in 1:length(S))
   {
     vine.knot <-  list(j1=S[i],j2=NULL,D=NULL,v=NULL, U = U[,S[i]],follow=c())
     SSi <- append(SSi, list(vine.knot))
   }
vine <- append(vine, list(SSi))
rm(SSi)
level <-  1
log.like <- AIC <- cAIC <- 0
pairs.fit <- get("pairs.fit",help.env)
pairs.new <- get("pairs.new",help.env)
mcoptions<-list(preschedule=FALSE)
vine <- append(vine, list(foreach(i=2:length(S),.options.multicore=mcoptions) %dopar% {
   index.j1 <- pairs.fit[i-1,1]
   index.j2 <- pairs.fit[i-1,2]
   if(selec=="cAIC") {
     for(kk in 1:dim(pairs.new)[1]) {
        if(all(pairs.new[kk,]==pairs.fit[i-1,])) model.l <- get("fit.level2",help.env)[[kk]]
     }
   }
   if(selec=="ken.tau") model.l <- get("fit.level2",help.env)[[i-1]]
   dim2<-length(model.l)
   vine.knot <-  list(j1=index.j1,j2=index.j2,D=NULL,v=model.l$v,U=U[,c(index.j1,index.j2)],log.like = model.l$log.like, AIC=model.l$AIC,cAIC=model.l$cAIC,lambda=model.l$lambda,f.hat=model.l$f.hat,indep=model.l$indep,knots1=model.l$knots1,knots2=model.l$knots2,int.base1=model.l$int.base1,int.base2=model.l$int.base2,Index.basis.D=model.l$Index.basis.D,ddb=model.l$ddb,K=model.l$K)
}))
foreach(i=2:length(S),.options.multicore=mcoptions) %dopar% {
   index.j1 <- pairs.fit[i-1,1]
   index.j2 <- pairs.fit[i-1,2]
   vine[[level]][[index.j1]]$follow <- c(vine[[level]][[index.j1]]$follow,i-1)
   vine[[level]][[index.j2]]$follow <- c(vine[[level]][[index.j2]]$follow,i-1)
   vine[[level+1]][[i-1]]$before <- c(vine[[level]][[i-1]]$before,pairs.fit[i-1,1],pairs.fit[i-1,2])
}

rm(list=c("fit.level2"),envir=help.env)
level <-  2
assign("level",level,help.env)
for(i in 1:length(vine[[level]])) {
  log.like <- log.like+vine[[level]][[i]]$log.like
  AIC <- AIC+vine[[level]][[i]]$AIC
  cAIC <- cAIC+vine[[level]][[i]]$cAIC
}
assign("rest.indep",FALSE,help.env)
  for(level in 3:length(S)) {
    #print(paste("level=",level,sep=""))
    assign("level",level,help.env)
    Tree.l1 <- vine[[level-1]]
    if(level<length(S)) {
      order.vine.level(tree=Tree.l1,help.env)
      pairs.fit <- get("pairs.fit",help.env)
      pairs.new <- get("pairs.new",help.env)
      vine <- append(vine, vinepart <- list(foreach(i=1:dim(pairs.fit)[1],.options.multicore=mcoptions) %do%  {
         help.j1 <- c(Tree.l1[[pairs.fit[i,1]]]$j1,Tree.l1[[pairs.fit[i,1]]]$j2,Tree.l1[[pairs.fit[i,1]]]$D)
         help.j2 <- c(Tree.l1[[pairs.fit[i,2]]]$j1,Tree.l1[[pairs.fit[i,2]]]$j2,Tree.l1[[pairs.fit[i,2]]]$D)
         index.j1 <- help.j1[!help.j1%in%help.j2]#Vorgaenger
         index.j2 <- help.j2[!help.j2%in%help.j1]#
         D.help <- c(Tree.l1[[pairs.fit[i,1]]]$D,Tree.l1[[pairs.fit[i,2]]]$D)
         D <- sort(unique(c(help.j1[help.j1%in%help.j2],D.help)))
         if(selec=="cAIC") {
           for(kk in 1:dim(pairs.new)[1]) {
              if(all(pairs.new[kk,]==pairs.fit[i,])) model.l <- get(paste("fit.level",level,sep=""),help.env)[[kk]]
           }
         }      
         if(selec=="ken.tau") model.l <- get(paste("fit.level",level,sep=""),help.env)[[i]] 
         list(j1=index.j1,j2=index.j2,D=D,v=get("ck.val",model.l),U=get("Y",model.l),log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.hat=get("f.hat.val",model.l),knots1=get("knots1",model.l),knots2=get("knots2",model.l),int.base1=get("int.base1",model.l),int.base2=get("int.base2",model.l),Index.basis.D=get("Index.basis.D",model.l),ddb=get("ddb",model.l),K=get("K",model.l))
      }))
      foreach(i=1:dim(pairs.fit)[1],.options.multicore=mcoptions,.multicombine=TRUE) %do% {
         index.j1 <- pairs.fit[i,1]
         index.j2 <- pairs.fit[i,2]
         vine[[level-1]][[index.j1]]$follow <- c(vine[[level-1]][[index.j1]]$follow,i)
         vine[[level-1]][[index.j2]]$follow <- c(vine[[level-1]][[index.j2]]$follow,i)
         vine[[level]][[i]]$before <- c(vine[[level]][[i]]$before,pairs.fit[i,1],pairs.fit[i,2])
      }
      add.VineMatrix(vinepart[[1]],help.env,level)
      for(i in 1:length(vine[[level]])) {
        log.like <- log.like+vine[[level]][[i]]$log.like
        AIC <- AIC+vine[[level]][[i]]$AIC
        cAIC <- cAIC+vine[[level]][[i]]$cAIC
      }
      rm(list=c(paste("fit.level",level,sep="")),pos=help.env)
      if(all(foreach(i=1:length(vine[[level]]),.combine=c) %do% vine[[level]][[i]]$indep)) assign("rest.indep",TRUE,help.env)
      #print(paste("level ",get("rest.indep",help.env),sep=""))
   }
   if(level==length(S)) {
      len.vine.before <- length(vine[[level-1]])
      pairs.fit <- matrix(c(1,2),1,2)
      help.j1 <- c(Tree.l1[[pairs.fit[1,1]]]$j1,Tree.l1[[pairs.fit[1,1]]]$j2,Tree.l1[[pairs.fit[1,1]]]$D)
      help.j2 <- c(Tree.l1[[pairs.fit[1,2]]]$j1,Tree.l1[[pairs.fit[1,2]]]$j2,Tree.l1[[pairs.fit[1,2]]]$D)
      j1 <- help.j1[!help.j1%in%help.j2]
      j2 <- help.j2[!help.j2%in%help.j1]
      D.help <- c(Tree.l1[[pairs.fit[1,1]]]$D,Tree.l1[[pairs.fit[1,2]]]$D)
      D <- sort(unique(c(help.j1[help.j1%in%help.j2],D.help)))
      index <- list(c(j1,D),c(j2,D)) # list of involved indices in current knot, i.e {j1,D} and {j2,D}
      UU <- c()
      for(j in 1:2) # Loop over both index sets, i.e. {j1,D} and {j2,D}
        {
          indexi <- sort(index[[j]])
          index.ancestor <-c()
          for (ml in 1:length(Tree.l1))
            {
              index.ancestor <- c(index.ancestor, all(indexi==sort(c(Tree.l1[[ml]]$j1,Tree.l1[[ml]]$j2,Tree.l1[[ml]]$D))))
            }
          ancestor.knot <- Tree.l1[index.ancestor][[1]] # Ancestor knot of j-th Element in Index set, i.e. knot with indices {j1,D} for j=1 and {j2,D} for j=2, respectively.
          
          if(j==1) {
            diff.help <- c("u1","u2")[!(c(Tree.l1[[pairs.fit[1,1]]]$j1,Tree.l1[[pairs.fit[1,1]]]$j2)%in%j1)]
          }
          if(j==2) {
            diff.help <- c("u1","u2")[!(c(Tree.l1[[pairs.fit[1,2]]]$j1,Tree.l1[[pairs.fit[1,2]]]$j2)%in%j2)]
          }
          UU <-cbind(UU,cond.cop(ancestor.knot$U,coef=ancestor.knot$v,K=ancestor.knot$K,diff=diff.help,ancestor.knot$Index.basis.D,base=base,q=q,env=help.env,kn1=ancestor.knot$knots1,kn2=ancestor.knot$knots2,int.base1=ancestor.knot$int.base1,int.base2=ancestor.knot$int.base2,ddb=ancestor.knot$ddb))
      }
      if(any(UU>1)) UU[which(UU>1)] <-1
      if(any(UU<0)) UU[which(UU<0)] <-0

          if(!l.search) model.l <- paircopula(data=UU,K=K,lambda=lambda,pen=pen,base=base,m=m,q=q,fix.lambda=fix.lambda,max.iter=get("max.iter",help.env))
          if(l.search) model.l<-lam.search(data=UU,K=get("K",help.env),lam=get("lam.vec",help.env),m=get("m",help.env),max.iter=get("max.iter",help.env),q=get("q",help.env),base=get("base",help.env),fix.lambda=fix.lambda,pen=pen,id=id)
          vine[[level-1]][[1]]$follow<-1
          vine[[level-1]][[2]]$follow<-1
          vine <- append(vine, vinepart <- list(list(list(j1=j1,j2=j2,D=D,v=get("ck.val",model.l),U=UU,log.like = get("log.like",model.l), AIC=get("AIC",model.l),cAIC=get("cAIC",model.l),lambda=get("lambda",model.l),f.hat=get("f.hat.val",model.l),knots1=get("knots1",model.l),knots2=get("knots2",model.l),int.base1=get("int.base1",model.l),int.base2=get("int.base2",model.l),Index.basis.D=get("Index.basis.D",model.l),ddb=get("ddb",model.l),K=get("K",model.l)))))#,p.value=test.res))))
        #}
      #}
      vine[[level]][[1]]$before<-c(1,2)
      add.VineMatrix(vinepart[[1]],help.env,level)
      log.like <- log.like+vine[[level]][[1]]$log.like
      AIC <- AIC+vine[[level]][[1]]$AIC
      cAIC <- cAIC+vine[[level]][[1]]$cAIC
   }
  }
  RVineMatrix(help.env)
  data.est <- rep(1,dim(data)[1])
  for(i in 2:length(S))  {
    for(j in 1:length(vine[[i]]))  {
      data.est <- data.est*vine[[i]][[j]]$f.hat
    }
  }
  class(vine) <- "penVine"
return(list(vine=vine,log.like=log.like,AIC=AIC,cAIC=cAIC,K=K,S=S,N=N,base=base,q=q,RVineMatrix=get("VineMatrix",help.env),env=help.env,f.est=data.est,selec=selec))
}
