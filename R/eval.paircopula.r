eval.paircopula <- function(val,K,int=FALSE,Index.basis.D,ck.val,base,q,knots.place,ddb,int.base1,int.base2,kn1,kn2) {
  if(!is.matrix(val)) {
    if(is.data.frame(val)) val <- as.matrix(val) else stop("val has to be a data.frame or a matrix")
  }
  tilde.Psi.d <-  array(NA, dim=c((length(val)/2),ddb,2))
  val <- matrix(val,(length(val)/2),2)
  penden.env <- new.env()
  assign("K",K,penden.env)
  assign("int.base1",int.base1,penden.env)
  assign("int.base2",int.base2,penden.env)
  assign("ddb",ddb,penden.env)

  index.b <- matrix(0:K)
  if(base=="Bernstein") {
    int.bernstein.help(penden.env)
    for (j in 1:2)
      {
        if(int) tilde.Psi.d[,,j] <-  int.bernstein(penden.env,Y=val[,j])
        else tilde.Psi.d[,,j] <- apply(index.b,1,bernstein,x=val[,j],n=K)
      }
  }
  if(base=="B-spline") {
    assign("n",dim(val)[1],penden.env)
    assign("q",q,penden.env)
    assign("x.help",seq(0, 1, length = 501),penden.env)
    for (j in 1:2)
      {
        if(int) tilde.Psi.d[,,j] <- int.bspline2(penden.env,Y=val[,j],k=j)
        else {
           if(j==1) kn<-kn1
           if(j==2) kn<-kn2
           tilde.Psi.d[,,j] <- my.bspline(y=val[,j],K=K,q=q,kn=kn)$base.den
        }
      }
  }
  val<-(tilde.Psi.d[,Index.basis.D[,1],1]*tilde.Psi.d[,Index.basis.D[,2],2])%*%ck.val
  val[val<0]<-0
  return(val)
}
