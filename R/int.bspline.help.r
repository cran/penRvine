int.bspline.help <- function(penden.env) {
  K <- get("K",penden.env)
  q<-get("q",penden.env)
  x.help <- seq(0,1,length=501)
  assign("int.base1",distr.func.help(base=get("base",penden.env),knots=get("knots1",penden.env),penden.env,q=get("q",penden.env),y=seq(0,1,length=501)),penden.env)
  assign("x.help",x.help,penden.env)
  assign("int.base2",get("int.base1",penden.env),penden.env)
}
