cond.cop <- function(data,coef,K,diff="u2",Index.basis.D,base,q=2,env,kn1,kn2,int.base1,int.base2,ddb) {
  p <- 2
  u1<-data[,1]
  u2<-data[,2]
  if(base=="B-spline"){
    assign("int.base1",int.base1,env)
    assign("int.base2",int.base2,env)
    assign("x.help",seq(0,1,length=501),env)
  }

  tilde.Psi.d <-  array(NA, dim=c(dim(data)[1],ddb,p))
  for (j in 1:p)
    {
      obj <- paste("u",j,sep="")
      if(base=="Bernstein") {
        if(obj==diff) tilde.Psi.d[,,j] <- apply(get("index.k",env),1,bernstein,eval(parse(text=(paste("u",j,sep="")))),n=get("K",env))
        if(obj!=diff) tilde.Psi.d[,,j] <- int.bernstein(env,eval(parse(text=(paste("u",j,sep="")))),k=j)
      }
      if(base=="B-spline") {
        if(obj==diff) {
          if(diff=="u2") kn<-kn2
          if(diff=="u1") kn<-kn1
          tilde.Psi.d[,,j] <- my.bspline(y=eval(parse(text=(paste("u",j,sep="")))),K=get("K",env),q=get("q",env),kn=kn)$base.den
        }
        if(obj!=diff) {
          tilde.Psi.d[,,j] <- int.bspline2(env,Y=eval(parse(text=(paste("u",j,sep="")))),k=j)
          ind<-which(tilde.Psi.d[,,j]>1)
          tilde.Psi.d[,,j][ind]<-1
        }
      }
    }
  val<-(tilde.Psi.d[,Index.basis.D[,1],1]*tilde.Psi.d[,Index.basis.D[,2],2])%*%coef
  val[which(val<0)]<-0
  return(round(val,6))
}

