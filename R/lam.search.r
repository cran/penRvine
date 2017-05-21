lam.search<-function(data,K,lam,m,max.iter,q,base,fix.lambda,pen,id) {
    jj<-1
    res<-foreach(jj=1:length(lam),.combine=rbind) %do% {
       cop<- try(paircopula(data=data,K=K,max.iter=max.iter,lambda=lam[jj],m=m,base=base,q=q,pen=pen,fix.lambda=fix.lambda))
       save(file=paste("cop.",jj,".id=",id,".Rdata",sep=""),list=c("cop"))
       if(class(cop)!="try-error") {
           val<- c(lam[jj],get("log.like",cop),get("cAIC",cop),get("i",cop))
       } 
       else val <- c(lam[jj],0,0,0)
       val
    }
    ind<-which.min(res[,3])
    set<-seq(1,length(lam))[-ind]
    load(paste("cop.",ind,".id=",id,".Rdata",sep=""))
    for(jk in 1:length(lam)) system(paste("rm cop.",jk,".id=",id,".Rdata",sep=""))
    return(cop)
}
