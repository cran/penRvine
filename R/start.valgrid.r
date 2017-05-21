start.valgrid <- function(penden.env) {

  ddb<-get("ddb",penden.env)
  q <- get("q",penden.env)
  p <- get("p",penden.env)
  X.knots <- matrix(NA,get("DD",penden.env),p)
  tilde.Psi.d.knots.start.r <-  array(NA, dim=c(dim(X.knots)[1],get("ddb",penden.env),p))
  if(get("base",penden.env)=="B-spline") {
    for(kk in 1:2) {
      if(get("q",penden.env)==2|get("q",penden.env)==1) {
         knots.help1 <- get("knots1",penden.env)
         knots.help2 <- get("knots2",penden.env)
         assign("knots.help1",knots.help1,penden.env)
         assign("knots.help2",knots.help2,penden.env)

      }
      if(get("q",penden.env)>=3) {
        x.help<-c()
        knots<-get(paste("knots",kk,sep=""),penden.env)
        for(j in 1:(length(get(paste("knots",kk,sep=""),penden.env))-1)) {
            x.help <- c(x.help,(knots[j+1]-knots[j])/2+knots[j])
        }
        if(get("q",penden.env)==3) assign(paste("knots.help",kk,sep=""),c(0,x.help,1),penden.env)

        if(get("q",penden.env)==4) {
        x.help2<-c()
        x.help<-c(0,x.help,1)
        for(j in 1:(length(x.help)-1)) {
            x.help2 <- c(x.help2,(x.help[j+1]-x.help[j])/2+x.help[j])
        }
        assign(paste("knots.help",kk,sep=""),c(0,x.help2,1),penden.env)
        }
      }
    }   
    for(j in 1:p)  X.knots[,j]  <- get(paste("knots.help",j,sep=""),penden.env)[get("Index.basis.D",penden.env)[,j]]
    if(is.null(get("length.cond",penden.env))) length.cond <- length(get("knots.help1",penden.env))
  }
  if(get("base",penden.env)=="Bernstein") {
     assign("knots.help1",get("knots1",penden.env),penden.env)
     assign("knots.help2",get("knots2",penden.env),penden.env)
     for(j in 1:p)  X.knots[,j]  <- get(paste("knots.help",j,sep=""),penden.env)[get("Index.basis.D",penden.env)[,j]]
     if(is.null(get("length.cond",penden.env))) length.cond <- length(get("knots.help1",penden.env))
  }

  env.extend <- list()
  for(j in 1:p) {
    name <- c(paste("y",j,sep=""))
    env.extend[[noquote(name)]] <- get(paste("knots.help",j,sep=""),penden.env)
  }
  assign("X.knots.g.all",expand.grid(env.extend),penden.env)
  tilde.Psi.d.knots.start.g.all <- array(NA, dim=c(dim(get("X.knots.g.all",penden.env))[1],get("ddb",penden.env),p))
  
  assign("X.knots",X.knots,penden.env)
  for (j in 1:p)
    {
      if(get("base",penden.env)=="Bernstein") {
        index.b <- matrix(0:get("dd",penden.env))
        tilde.Psi.d.knots.start.r[,,j] <-  apply(index.b,1,bernstein,x=X.knots[,j],n=get("dd",penden.env))
        tilde.Psi.d.knots.start.g.all[,,j] <-  apply(index.b,1,bernstein,x=get("X.knots.g.all",penden.env)[,j],n=get("dd",penden.env))
      }
      if(get("base",penden.env)=="B-spline") {
        tilde.Psi.d.knots.start.r[,,j] <-  my.bspline(y=X.knots[,j],K=get("K",penden.env),q=get("q",penden.env),kn=get(paste("knots",j,sep=""),penden.env))$base.den
        tilde.Psi.d.knots.start.g.all[,,j] <-  my.bspline(y=get("X.knots.g.all",penden.env)[,j],K=get("K",penden.env),q=get("q",penden.env),kn=get(paste("knots",j,sep=""),penden.env))$base.den
      }
    }  
   assign("tilde.PSI.d.D.knots.start.r",tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,1],1],penden.env)
  assign("tilde.PSI.d.D.knots.start.g.all",tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,1],1],penden.env)
                                       
  for (j in 2:p) {
    assign("tilde.PSI.d.D.knots.start.r",get("tilde.PSI.d.D.knots.start.r",penden.env) * tilde.Psi.d.knots.start.r[,get("Index.basis.D",penden.env)[,j],j],penden.env)
    assign("tilde.PSI.d.D.knots.start.g.all",get("tilde.PSI.d.D.knots.start.g.all",penden.env) * tilde.Psi.d.knots.start.g.all[,get("Index.basis.D",penden.env)[,j],j],penden.env)
   }
  assign("ck.val",matrix(solve(get("tilde.PSI.d.D.knots.start.r",penden.env),rep(1,get("DD",penden.env)),tol=1e-50)),penden.env)
  assign("ck.val.start",get("ck.val",penden.env),penden.env)
  assign("ck.val.temp",get("ck.val",penden.env),penden.env)
}
