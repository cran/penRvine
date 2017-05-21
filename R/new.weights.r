new.weights <- function(penden.env,lambda.temp=NULL,start=FALSE) {
  dd <- get("dd",penden.env)
  p <- get("p",penden.env)
  DD <- get("DD",penden.env)
  calc<-TRUE
  l.A <- length(get("T.marg",penden.env)[,1,1])
  l.A2 <- length(get("T.marg",penden.env)[,1,2])
  if(l.A>1) vec <- seq(1,(l.A-1)) else vec<-1
  if(l.A2>1) vec2 <- seq(1,(l.A2-1)) else vec2<-1
  #eps <- 1e-14
  eps<-get("eps.we",penden.env)
  ind.val<-get("ind.val",penden.env)
  assign("help.lambda2",NULL,penden.env)
  if (get("base", penden.env) == "Bernstein") {
    if(is.null(ind.val)) {
      assign("AA.help", t(get("T.marg", penden.env)[vec, ,1]), penden.env)
      assign("AA.help", cbind(get("AA.help", penden.env), t(get("T.marg", penden.env)[vec2, , 2])), penden.env)
    }
    if(!is.null(ind.val)) {
      assign("AA.help", t(get("T.marg", penden.env)[vec, ,1][,-ind.val]), penden.env)
      assign("AA.help", cbind(get("AA.help", penden.env), t(get("T.marg", penden.env)[vec2, , 2][,-ind.val])), penden.env)
    }
  }
  if (get("base", penden.env) == "B-spline")  {
     if(is.null(ind.val)) {
         assign("AA.help", t(get("A.Restrict", penden.env)[vec, , 1]), penden.env)
         assign("AA.help", cbind(get("AA.help", penden.env), t(get("A.Restrict", penden.env)[vec2, , 2])), penden.env)
     }
     if(!is.null(ind.val)) {
         assign("AA.help", t(get("A.Restrict", penden.env)[vec, , 1])[-ind.val,], penden.env)
         assign("AA.help", cbind(get("AA.help", penden.env), t(get("A.Restrict", penden.env)[vec2, , 2])[-ind.val,]), penden.env)
     }
  }
  meq <- 1+length(vec)+length(vec2)

  if(get("base",penden.env)=="B-spline") {
    if(is.null(ind.val)) bvec <- c(rep(0,meq),-get("ck.val",penden.env)+eps)
    if(!is.null(ind.val)) bvec <- c(rep(0,meq),-get("ck.val",penden.env)[-ind.val]+eps)
    if(is.null(ind.val)) assign("Amat",cbind(matrix(1,DD,1),get("AA.help",penden.env),diag(1,DD)),penden.env)
    if(!is.null(ind.val)) assign("Amat",cbind(matrix(1,DD-length(ind.val),1),get("AA.help",penden.env),diag(1,DD-length(ind.val))),penden.env)
  } 
  if(get("base", penden.env) == "Bernstein") {
     if(is.null(ind.val)) {
       bvec <- c(0, rep(1/(get("K", penden.env) + 1), length(vec2)+ length(vec)) - t(get("AA.help", penden.env)) %*% get("ck.val", penden.env),-get("ck.val", penden.env)+eps)
       assign("Amat", cbind(matrix(1, DD, 1), get("AA.help", penden.env), diag(1, DD)), penden.env)
     }
     if(!is.null(ind.val)) {
       bvec <- c(0, rep(1/(get("K", penden.env) + 1), length(vec2)+ length(vec)) - t(get("AA.help", penden.env)) %*% get("ck.val", penden.env)[-ind.val],-get("ck.val", penden.env)[-ind.val]+eps)
       assign("Amat", cbind(matrix(1, DD-length(ind.val), 1), get("AA.help", penden.env), diag(1, DD-length(ind.val))), penden.env)
     }

  }
  Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
  Derv2(penden.env,temp=TRUE,lambda=lambda.temp)
  sc <- base::norm(-get("Derv2.pen.temp",penden.env),"2")
  Derv2.pen.help <- -get("Derv2.pen.temp",penden.env)/sc
  Derv1.help<- get("Derv1.pen.temp",penden.env)/sc

  aa <- try(obj <- solve.QP(Dmat=Derv2.pen.help,dvec=Derv1.help,Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE),silent=FALSE)

  if((class(aa)=="try-error")&start) {
   Dv1<- get("Derv1.pen.temp",penden.env)
   ind0<-which(Dv1<0&Dv1>-1e-8)  
   Dv1[ind0]<-0
   Derv1.help<- Dv1/sc
   aa <- try(obj <- solve.QP(Dmat=Derv2.pen.help,dvec=Derv1.help,Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE),silent=FALSE)
  }

  if(class(aa)=="try-error") {
      if(any(-get("Derv1.pen.temp",penden.env)>0)) {
        print("AB")
        c.vec<-seq(0.01,0.99,by=0.01)
        ii<-1
        res1<-foreach(ii=length(c.vec):1,.combine=rbind) %do% {
           lam.tem<-lambda.temp*c.vec[ii]
           Derv1(penden.env,temp=TRUE,lambda=lam.tem)
           Derv2(penden.env,temp=TRUE,lambda=lam.tem)
           sc <- base::norm(-get("Derv2.pen.temp",penden.env),"2")
           Derv2.pen.help <- -get("Derv2.pen.temp",penden.env)/sc
           Derv1.help<- get("Derv1.pen.temp",penden.env)/sc
           aa <- try(obj <- solve.QP(Dmat=Derv2.pen.help,dvec=Derv1.help,Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE),silent=FALSE)
           if(class(aa)=="try-error") vali<-c(c.vec[ii],FALSE)
           else vali<-c(c.vec[ii],TRUE)
        }
        lambda.temp<-max(res1[res1[,2],1])*lambda.temp
        Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
        Derv2(penden.env,temp=TRUE,lambda=lambda.temp)
        sc <- base::norm(-get("Derv2.pen.temp",penden.env),"2")
        Derv2.pen.help <- -get("Derv2.pen.temp",penden.env)/sc
        Derv1.help<- get("Derv1.pen.temp",penden.env)/sc
        assign("help.lambda2",lambda.temp,penden.env)
        aa <- try(obj <- solve.QP(Dmat=Derv2.pen.help,dvec=Derv1.help,Amat=get("Amat",penden.env),bvec=bvec,meq=meq,factorized=FALSE),silent=FALSE)
     }
     else return("fehler")
  }
  assign("delta",obj$solution,penden.env)

  if(get("i",penden.env)>1& any(is.na(obj$solution))) {
     dimA<-dim(get("Amat",penden.env))
     aa <- try(obj <- solve.QP(Dmat=Derv2.pen.help,dvec=Derv1.help,Amat=get("Amat",penden.env)[,-dimA[2]],bvec=bvec[-dimA[2]],meq=meq,factorized=FALSE),silent=FALSE)
     if(any(is.na(obj$solution))) {
       print("is.na in obj$solution")
       assign("wrong.lambda",TRUE,penden.env)
       return("fehler")
     }
  }

  if(is.null(get("ind.val",penden.env))) assign("ck.val.temp",get("ck.val",penden.env)+obj$solution,penden.env)
  if(!is.null(get("ind.val",penden.env))) {
    ck.val<-get("ck.val",penden.env)
    ind.seq<-seq(1,length(ck.val))
    ck.val[ind.seq[-ind.val]]<-ck.val[-ind.val]+obj$solution
    ck.val[ck.val<=0]<-0
    ind<-which(ck.val<=0)
    #print(ind)
    assign("ind.val",ind,envir=penden.env) 
    assign("ck.val.temp",ck.val,penden.env)
  }

  if(any(get("ck.val.temp",penden.env)<=0)&is.null(ind.val)) {
     ck.val.temp<-get("ck.val.temp",penden.env)
     ind<-which(ck.val.temp<=0)
     #print(ind)
     ck.val.temp[ind]<-0
     assign("ck.val.temp",ck.val.temp,penden.env)
     assign("ind.val",ind,envir=penden.env) 
  }

  Derv1(penden.env,temp=TRUE,lambda=lambda.temp)
  Derv2(penden.env,temp=TRUE,lambda=lambda.temp)
  f.hat.val(penden.env,temp=TRUE)
  pen.log.like(penden.env,temp=TRUE)
  assign("help.lambda",lambda.temp,penden.env)
  if(get("no",penden.env)) return("fehler")
  else return("kein fehler")
}
