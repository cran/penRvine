Derv1 <- function(penden.env,temp=FALSE,lambda=NULL,lam.fit=FALSE) {
  ind.val<-get("ind.val",penden.env)
  if(temp) {
    ck<-get("ck.val.temp",penden.env)
    Fy<-get("f.hat.val.temp",penden.env)
  }
  if(!temp) {
    ck<-get("ck.val",penden.env) 
    Fy<-get("f.hat.val",penden.env)
  }
  assign("tilde.PSI.d.D.t",get("tilde.PSI.d.D",penden.env),penden.env)  
  if(any(Fy<1e-12)) {
    ind<-which(Fy<1e-12)
    Fy<-get("tilde.PSI.d.D",penden.env)[-ind,]%*% ck
    assign("tilde.PSI.d.D.t",get("tilde.PSI.d.D",penden.env)[-ind,],penden.env)  
  }
  if(!temp) {
    if(is.null(ind.val)) assign("Derv1.pen",matrix(colSums(get("tilde.PSI.d.D.t",penden.env)/kronecker(Fy, matrix(1,1,dim(get("tilde.PSI.d.D.t",penden.env))[2]))),get("DD",penden.env),1)-get("lambda",penden.env)*get("DDD.sum",penden.env)%*%ck,penden.env)
    if(!is.null(ind.val)) assign("Derv1.pen",matrix(colSums(get("tilde.PSI.d.D.t",penden.env)[,-ind.val]/kronecker(Fy, matrix(1,1,dim(get("tilde.PSI.d.D.t",penden.env)[,-ind.val])[2]))),get("DD",penden.env)-length(ind.val),1)-get("lambda",penden.env)*get("DDD.sum",penden.env)[-ind.val,-ind.val]%*%ck[-ind.val],penden.env)
  }
  if(temp) {
    if(is.null(lambda)) lambda <- get("lambda",penden.env)
    if(is.null(ind.val)&!lam.fit) assign("Derv1.pen.temp",matrix(colSums(get("tilde.PSI.d.D.t",penden.env)/kronecker(Fy, matrix(1,1,dim(get("tilde.PSI.d.D.t",penden.env))[2]))),get("DD",penden.env),1)-lambda*get("DDD.sum",penden.env)%*%ck,penden.env)
    if(!is.null(ind.val)&!lam.fit) assign("Derv1.pen.temp",matrix(colSums(get("tilde.PSI.d.D.t",penden.env)[,-ind.val]/kronecker(Fy, matrix(1,1,dim(get("tilde.PSI.d.D.t",penden.env)[,-ind.val])[2]))),get("DD",penden.env)-length(ind.val),1)-lambda*get("DDD.sum",penden.env)[-ind.val,-ind.val]%*%ck[-ind.val],penden.env)
   if(lam.fit) assign("Derv1.pen.temp",matrix(colSums(get("tilde.PSI.d.D.t",penden.env)/kronecker(Fy, matrix(1,1,dim(get("tilde.PSI.d.D.t",penden.env))[2]))),get("DD",penden.env),1)-lambda*get("DDD.sum",penden.env)%*%ck,penden.env)
   } 
}
 
