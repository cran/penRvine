int.bernstein <- function(penden.env,Y,k) {
  int.base <- matrix(0,length(Y),(get("K",penden.env)+1))
  index.k <- matrix(1:length(Y))
  for(j in 1:(get("K",penden.env)+1)) int.base[,j] <- spline(get(x="x.help",penden.env),y=get(paste("int.base",k,sep=""),penden.env)[,j],xout=Y)$y
  return(int.base)
}
