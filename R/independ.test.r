independ.test<-function(data) {
   #cond.level<-0.075
   #ord <- order(UU[,1])
   N<-dim(data)[1] 
   #if(N%%2!=0) ord <- ord[-round(N/2,0)]
   #st1<-UU[ord[(1:(length(ord)/2))],]
   #st2<-UU[ord[((length(ord)/2+1):length(ord))],]
   #library(TwoCop)
   #set.seed(27)
   #if(N>=200) sz<-100 else sz<-round((N/2)/100,0)*100
   #teste<-foreach(j=1:39,.combine=rbind,.multicombine=TRUE) %do% {
   #      prob1<-sample(x=c(1:(N/2)),size=sz,replace=FALSE)
   #      prob2<-sample(x=c(1:(N/2)),size=sz,replace=FALSE)
   #         #prob1<-sample(x=c(1:N),size=sz,replace=FALSE)
   #         #prob2<-sample(x=c(1:N),size=sz,replace=FALSE)
   #      TwoCop(x=st1[prob1,],y=st2[prob2,],alpha=0.95,paired=TRUE)$pvalue
   #}
   #print(c("test.results",sort(teste),apply(teste,2,median)))
   #test.res <- apply(teste,2,median)
   #if(test.res<0.05) return(FALSE) else return(TRUE)


###pre-test#####
    ken.tau<-cor(data[,1],data[,2],method="kendall")
    T <- ( (9*N*(N-1)) / (2*(2*N+5)) )^0.5 * abs(ken.tau)
    pval<-(2*(1-pnorm(T)))
    if(pval<=0.05) return(pval)
    else return(indepTest(data,indepTestSim(n=dim(data)[1],p=dim(data)[2],m=2,N=100))$pvalues)
###########

   #d <- indepTestSim(n=dim(data)[1],p=dim(data)[2],m=2,N=100)
   #test <- indepTest(data,indepTestSim(n=dim(data)[1],p=dim(data)[2],m=2,N=100))
   #return(indepTest(data,indepTestSim(n=dim(data)[1],p=dim(data)[2],m=2,N=100))$pvalues)

}
