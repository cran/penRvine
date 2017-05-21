plot.paircopula <- function(x,val=NULL,marg=TRUE,plot=TRUE,int=FALSE,main.txt=NULL,
                    sub.txt=NULL,contour=FALSE,cuts=20,cex=1,cex.axes=1,
                    xlab=NULL,ylab=NULL,zlab=NULL,xlim=NULL,ylim=NULL,zlim=NULL,margin.normal=FALSE,...) {
  if(!class(x)=="paircopula") stop("obj has to be of class paircopula")
  env <- list()
  p <- get("p",x)
  d <- get("d",x)
  ddb <- get("ddb",x)
  Index.basis.D <- get("Index.basis.D",x)
  ck <- get("ck.val",x)
  D <- get("K",x)
  alpha <- 0
  if(margin.normal) assign("margin.normal",TRUE,x)
  base <- get("base",x)
  index.b <- matrix(0:get("dd",x))

  if(is.null(val))
    {
      if(p!=2) stop("geht nicht!")
      else {
        x.grid <- seq(0,1,length=51)
        grid <- expand.grid(y1=x.grid, y2=x.grid)  
        if(margin.normal) {
               x.grid.unif <- seq(-5,5,length=51)
               x.grid <- pnorm(x.grid.unif)
               grid <- expand.grid(y1=x.grid, y2=x.grid) 
        }
        tilde.Psi.d <-  array(NA, dim=c(dim(grid)[1],get("ddb",x),p))
        for (j in 1:p)
          {
            if(base=="Bernstein") {
              if(int) tilde.Psi.d[,,j] <-  int.bernstein(x,Y=grid[,j],k=j)
              else tilde.Psi.d[,,j] <-  apply(index.b,1,bernstein,x=grid[,j],n=get("dd",x))
            }
            if(base=="B-spline") {
              if(int) tilde.Psi.d[,,j] <- int.bspline2(x,Y=grid[,j],k=j)
              else tilde.Psi.d[,,j] <- my.bspline(y=grid[,j],K=get("K",x),q=get("q",x),margin.normal=margin.normal,kn=get(paste("knots",j,sep=""),x))$base.den
            }  
          }
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
        
        for (j in 2:p)
          {        
            tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,j],j]
          }
        grid[["plot"]] <- tilde.PSI.d.D%*%ck
        grid[["plot"]][which(grid[["plot"]]<0)]<-0
        if(is.null(zlim)) zlim <- c(0,round(max(grid[["plot"]]),2))
        if(is.null(ylim)) ylim <- c(0,1)
        if(is.null(xlim)) xlim <- c(0,1)
         
        lam1 <- get("lambda",x)[1]               
        if(is.null(main.txt)) {
          main.txt <- substitute("K="*a*","*c*"="*d,list(a=D,c=parse(text="lambda")[[1]],d=lam1))
          main.txt <- as.expression(main.txt)
        }
        if(margin.normal) {
           plot1<-grid$plot
           grid <- expand.grid(y1=x.grid.unif, y2=x.grid.unif)
           grid[["plot"]]<-plot1*dnorm(grid$y1)*dnorm(grid$y2)
           if(is.null(xlim)) xlim<-c(-3,3)
           if(is.null(ylim)) ylim<-c(-3,3)
           if(is.null(zlim)) zlim <- c(0,round(max(grid[["plot"]]),2))
        }
        if(!margin.normal) {
          if(is.null(zlim)) zlim <- c(0,round(max(grid[["plot"]]),2))
          if(is.null(ylim)) ylim <- c(0,1)
          if(is.null(xlim)) xlim <- c(0,1)
        }
        k <- dim(x$liste)[1]
        log.like <- round(get("log.like",x),3)
        pen.log.like <- round(get("pen.log.like",x),3)
        if(is.null(sub.txt)) sub.txt <- paste("log like=",log.like,", pen. log like= ",pen.log.like,", cAIC=",round(get("cAIC",x),3),sep="")
        if(is.null(zlab)) {
             if(int) z.txt <- "distribution" 
             if(!int) z.txt <- "density"
        }

        hh <- c("y1","y2")
        values <- as.formula(paste("plot~",paste(hh,collapse="*"),sep=""))
        if(!contour) obj1 <- wireframe(values,data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(cex=cex.axes,
                                     label=xlab),ylab=list(cex=cex.axes,label=ylab),scales=list(arrows=FALSE,col="black",font=3,
                                     x=list(cex=cex),y=list(cex=cex),z=list(cex=cex)),main=main.txt,xlim=xlim,ylim=ylim,
                                     shade=TRUE,par.settings = list(axis.line = list(col = "transparent")),par.box = c(col = "black"),zlim=zlim)
        else obj1 <- contourplot(values, data=grid,outer=TRUE,sub=sub.txt,zlab=list(label=z.txt,cex=cex.axes),xlab=list(label=xlab,cex=cex.axes),
                                 ylab=list(label=ylab,cex=cex.axes),scales=list(arrows=FALSE,col="black",font=3,cex=cex,x=list(lim=xlim),
                                 y=list(lim=ylim)),zlim=zlim,main=main.txt,shade=TRUE,cuts=cuts)
        if(marg) {
          T.marg <- get("T.marg",x)
          if(base=="B-spline") K<-get("ddb",x)
          if(base=="Bernstein") K<-get("K",x)+1
          xx <- rep(seq(0,1,length=K),p)
          density <- c()
          fac <- c()
          
          for(j in 1:p)
            {
              if(base=="B-spline") density <- c(density,round(get(paste("tilde.Psi.knots.d",j,sep=""),x)%*%(T.marg[,,j]%*%ck),5))
              if(base=="Bernstein") density <- c(density,round(get("tilde.Psi.knots.d",x)%*%(T.marg[,,j]%*%ck),5))
              fac <- c(fac,rep(j,get("ddb",x)))
            }
          datafr <- data.frame(xx,density,fac)
          graph.sets <-list(superpose.line=list(col=c(1:p),superpose.symbol = list(col = c(1:p))))
          
          obj2 <- xyplot(density~xx|fac,type="l",auto.key=list(space="right",title="marginal densities",sort=FALSE),par.settings=graph.sets)
          if(plot) {
            print(obj1,position=c(0,0.35,1,1),more=TRUE)
            print(obj2,position=c(0,0,1,0.35))
          }
          else return(list(density=obj1,marg.density=obj2))
        }
        else if(plot) {
          print(obj1)
          check <- any(grid$plot<0)
          print(check)
          if(check) print(grid[grid$plot<0,])
        }
        else return(list(density=obj1,grid=grid))
      }
    }
  else
    {
      if(!is.matrix(val)) {
        if(is.data.frame(val)) val <- as.matrix(val) else stop("val has to be a data.frame or a matrix")
      }
        tilde.Psi.d <-  array(NA, dim=c((length(val)/p),get("ddb",x),p))
        val <- matrix(val,(length(val)/p),p)
        if(int) tilde.Psi.d[,,1] <- int.bspline2(x,Y=val[,1],k=1)
        if(!int) tilde.Psi.d[,,1] <- my.bspline(y=val[,1],K=get("K",x),q=get("q",x),kn=get("knots1",x))$base.den
        tilde.PSI.d.D <- tilde.Psi.d[,Index.basis.D[,1],1]
        if(int) tilde.Psi.d[,,2] <- int.bspline2(x,Y=val[,2],k=2)
        if(!int) tilde.Psi.d[,,2] <- my.bspline(y=val[,2],K=get("K",x),q=get("q",x),kn=get("knots2",x))$base.den
        tilde.PSI.d.D <- tilde.PSI.d.D * tilde.Psi.d[,Index.basis.D[,2],2]
        val2<-tilde.PSI.d.D%*%ck
        ind<-which(val2<1e-12)
        datafr <- data.frame(val,val2)
        colnames(datafr)[p+1] <- "fit"
        return(datafr)
      }
}
