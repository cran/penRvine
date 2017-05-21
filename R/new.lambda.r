new.lambda <- function(penden.env) {
  lambda <- get("lambda",penden.env)
  if(get("lambda",penden.env)!=0) {
    eps <- 0.0025*get("lambda",penden.env)
    #vorher mit 0.0025
  }
  else eps <-1
  #print(paste("eps",eps,sep=""))
  epsdf <- 0
  p <- get("p",penden.env)
  calc <- TRUE
  u <- t(get("ck.val.temp",penden.env))%*%get("DDD.sum",penden.env)%*%get("ck.val.temp",penden.env)
  pen.mat <- get("DDD.sum",penden.env)
  help2 <- get("eigen.pen.mat",penden.env)
  index <- get("index.eigen.pen.mat",penden.env)
  Utilde <- get("Utilde.eigen.pen.mat",penden.env)

  t.Utilde <- t(Utilde)
  diag.help2 <- diag(help2$values[index])
  hh <-1
  pos.def<-c(TRUE)
  while(calc) {
    #print(paste("lam=",lambda,".i=",hh,sep=""))
    if(hh==51) {
      assign("df.val",df.val,penden.env)
      #print("151")
      #print(lambda)
      break
      assign("lambda.out",TRUE,penden.env)
    }
    Derv2(penden.env,temp=TRUE,lambda=lambda[hh],lam.fit=TRUE)
    df.val <- sum(diag(x=solve(a=t.Utilde%*%(-get("Derv2.cal.temp",penden.env))%*%Utilde+lambda[hh]*diag.help2,tol=1e-50)%*%(t.Utilde%*%(-get("Derv2.cal.temp",penden.env))%*%Utilde)))
    aa<-try(lam.tem<-c(df.val/u))
    if(class(aa)=="try-error") return(lambda[hh])
    bb<-try(help1 <- abs(lam.tem - lambda[hh]))
    if(class(bb)=="try-error") browser()
    #help3 <- abs(df.val/u - lambda[1])
    #if((df.val/u)<0) {
    #  print("BLA")
    #  assign("df.val",df.val,penden.env)
    #  return(lambda[hh-1])
    #}
    if(help1<eps) {
      #print(paste("help1=",help1,sep=""))
      #lam.tem<-c(df.val/u)
      #Derv1(penden.env,temp=TRUE,lambda=lam.tem,lam.fit=TRUE)
      #if(any(-get("Derv1.pen.temp",penden.env)>0)) {
      #   print("AAA")
      #  c.vec<-c(0,0.0000001,0.000001,0.00001,0.0001,seq(0.001,0.99,by=0.005))
      #  res1<-foreach(ii=1:length(c.vec),.combine=c) %do% {
      #     Derv1(penden.env,temp=TRUE,lambda=lam.tem*c.vec[ii],lam.fit=TRUE)
      #     all(-get("Derv1.pen.temp",penden.env)<0)
      #  }
      #  if(length(c.vec[res1])>0) val1<-max(c.vec[res1])
      #  if(length(c.vec[res1])==0) {
      #    #print("C")
      #    return(lambda[hh])
      # }
      #}
      #else return(lam.tem)
      #Derv2(penden.env,temp=TRUE,lambda=lam.tem*val1,lam.fit=TRUE)
      #if(!any(eigen(-get("Derv2.pen.temp",penden.env))$value<=0)) pos.def<-c(pos.def,TRUE)
      #if(any(eigen(-get("Derv2.pen.temp",penden.env))$value<=0)) browser()
      assign("df.val",df.val,penden.env)
      #print("A")
      #print(lam.tem*val1)
      return(lam.tem)
    }
    #if((help3<eps)&hh>10) {
    #  assign("df.val",df.val,penden.env)
    #  #assign("lambda.out",TRUE,penden.env)
    #  lambda[hh+1]<-df.val/u
    #  #print("out")
    #  #print(lambda)
    #  return(lambda[3])
    #}
    #lam.tem<-c(df.val/u)
    #val1<-1
    #Derv1(penden.env,temp=TRUE,lambda=lam.tem,lam.fit=TRUE)
    #if(any(-get("Derv1.pen.temp",penden.env)>0)) {
    #  print("AAAB")
    #  browser()
    #  c.vec<-c(0,0.0000001,0.000001,0.00001,0.0001,seq(0.001,0.99,by=0.005))
    #  res1<-foreach(ii=1:length(c.vec),.combine=c) %do% {
    #     Derv1(penden.env,temp=TRUE,lambda=lam.tem*c.vec[ii],lam.fit=TRUE)
    #     all(-get("Derv1.pen.temp",penden.env)<0)
    #  } 
    # if(length(c.vec[res1])>0) val1<-max(c.vec[res1])
    #  if(length(c.vec[res1])==0) {
    #    print("non-regularA")
    #    return(lambda[hh])
    #  }
    #}
    #Derv2(penden.env,temp=TRUE,lambda=lam.tem*val1,lam.fit=TRUE)
    #if(!any(eigen(-get("Derv2.pen.temp",penden.env))$value<=0)) pos.def<-c(pos.def,TRUE)
    #if(any(eigen(-get("Derv2.pen.temp",penden.env))$value<=0)) browser()
    if(abs((lam.tem)-lambda[hh])>eps) {
      lambda[hh+1] <- lam.tem
      hh <- hh+1
    }
    #if(abs((lam.tem*val1)-lambda[hh])<eps) {  
    #  print("non-regularB")
    #  return(lam.tem*val1)
    #}
  }
  #print("regular")
  return(lambda[hh])
}
