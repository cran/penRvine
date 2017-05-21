marg.likelihood <- function(penden.env,pen.likelihood,temp=FALSE) {
  help <- get("eigen.pen.mat",penden.env)
  index <- get("index.eigen.pen.mat",penden.env)
  diag.help2 <- Matrix(diag(help$values[index]),sparse=FALSE,doDiag=TRUE)
  if(is.null(get("ind.val",penden.env))) {
    Utilde <- Matrix(get("Utilde.eigen.pen.mat",penden.env))
    evalues <- help$values[index]
    t.Utilde <- t(Utilde)
  }
  if(!is.null(get("ind.val",penden.env))) {
    ind.val<-get("ind.val",penden.env)
    pen.mat <- get("DDD.sum",penden.env)[-ind.val,-ind.val]
    help2 <- eigen(pen.mat)
    index <- which(help2$values>1e-16)
    Utilde <- help2$vectors[,index]
    evalues <- help2$values[index]
    t.Utilde <- t(Utilde)
  }
  k1 <- 0.5*sum(log(get("lambda",penden.env)*evalues))
  k2 <- pen.likelihood
  if(!temp) eneu <- eigen(t.Utilde%*%-get("Derv2.pen",penden.env)%*%Utilde)$values
  else eneu <- eigen(t.Utilde%*%-get("Derv2.pen.temp",penden.env)%*%Utilde)$values
  if(any(eneu<=0)) print("eneu <=0")
  k3 <- -0.5*sum(log(eneu))
  if(!temp) assign("marg.log.like",k1+k2+k3,penden.env)
  else assign("marg.log.like.temp",k1+k2+k3,penden.env)
}
