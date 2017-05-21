Derv2 <- function(penden.env,temp=FALSE,lambda=NULL,lam.fit=FALSE) {
  ind.val<-get("ind.val",penden.env)
  if(temp) {
    ck<-get("ck.val.temp",penden.env)
    Fy<-get("f.hat.val.temp",penden.env)
  }
  if(!temp) {
    ck<-get("ck.val",penden.env) 
    Fy<-get("f.hat.val",penden.env)
  }
  if(any(Fy<1e-12)) {
     ind<-which(Fy<1e-12)
     if(is.null(ind.val)) {
       if(!temp) {
         assign("Derv2.pen",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         #assign("Derv2.pen", -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*%ck)%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2])))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2])))),penden.env)
         #assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck)%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2]))),penden.env)
       }
       if(temp) {
         if(is.null(lambda)) lambda <- get("lambda",penden.env)
         assign("Derv2.pen.temp", -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         #assign("Derv2.pen.temp", -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*%ck)%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2])))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         assign("Derv2.cal.temp",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2])))),penden.env)
         #assign("Derv2.cal.temp", -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*%ck)%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,])[2]))),penden.env)
       }
     }
     if(!is.null(ind.val)) {
       if(!temp) {
         assign("Derv2.pen",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         #assign("Derv2.pen",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/((get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val])%*% matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2])))-get("lambda",penden.env)*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2])))),penden.env)
         #assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/((get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val])%*% matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2]))),penden.env)
       }
       if(temp) {
         if(is.null(lambda)) lambda <- get("lambda",penden.env)
         if(!lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2])))),penden.env)
         #if(!lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/((get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val])%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2]))),penden.env)
         if(!lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2]))))-lambda*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         #if(!lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]/((get("tilde.PSI.d.D",penden.env)[-ind,-ind.val]%*% ck[-ind.val])%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[-ind,-ind.val])[2])))-lambda*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         if(lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2])))),penden.env)
         #if(lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck)%*% matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))),penden.env)
         if(lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/(kronecker(get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))))-lambda*(get("DDD.sum",penden.env)),penden.env)
         #if(lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[-ind,]/((get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck)%*%matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2])))-lambda*(get("DDD.sum",penden.env)),penden.env)
       }
     }
  }
  else {
     if(is.null(ind.val)) {
       if(!temp) {
         assign("Derv2.pen", -crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2])))),penden.env)
       }
       if(temp) {
         if(is.null(lambda)) lambda <- get("lambda",penden.env)
         assign("Derv2.pen.temp", -crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)),penden.env)
         assign("Derv2.cal.temp",-crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2])))),penden.env)
       }
     }
     if(!is.null(ind.val)) {
       if(!temp) {
         assign("Derv2.pen",-crossprod(get("tilde.PSI.d.D",penden.env)[,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[,-ind.val])[2]))))-get("lambda",penden.env)*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         assign("Derv2.cal",-crossprod(get("tilde.PSI.d.D",penden.env)[,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[,-ind.val])[2])))),penden.env)
       }
       if(temp) {
         if(is.null(lambda)) lambda <- get("lambda",penden.env)
         if(!lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[,-ind.val])[2])))),penden.env)
         if(!lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)[,-ind.val]/(kronecker(get("tilde.PSI.d.D",penden.env)[,-ind.val]%*% ck[-ind.val], matrix(1,1,dim(get("tilde.PSI.d.D",penden.env)[,-ind.val])[2]))))-lambda*(get("DDD.sum",penden.env)[-ind.val,-ind.val]),penden.env)
         if(lam.fit) assign("Derv2.cal.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2])))),penden.env)
         if(lam.fit) assign("Derv2.pen.temp",cp<- -crossprod(get("tilde.PSI.d.D",penden.env)/(kronecker(get("tilde.PSI.d.D",penden.env)%*% ck, matrix(1,1,dim(get("tilde.PSI.d.D",penden.env))[2]))))-lambda*(get("DDD.sum",penden.env)),penden.env)
       }
     }
  }
}
