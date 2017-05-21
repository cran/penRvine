\name{penRvine-package}
\alias{penRvine-package}
\docType{package}
\title{Flexible R-vines Estimation Using Bivariate Penalized Splines
}
\description{
The package 'penRvine' offers routines for estimating densities and
distribution of R-vines. For details see the description of the function vine().
}
\details{
\tabular{ll}{
Package: \tab penRvine\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2017-05-19\cr
License: GPL (>= 2)
LazyLoad: \tab yes\cr
}
The packages contributes the function 'vine()' for estimating densities
and distributions of vines using penalized splines techniques. 
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\examples{
  #Simulating from a three-dimensional frank copula with
  #kendell's tau equal to 0.25, sample size N.set=100.
  #Please enlarge N.set for further studies.
  #require(copula)
  #N.set<-100
  #cop <- archmCopula(family = "frank", dim = 3, param =2.39)
  #parMarg<-list(list(min=0,max=1),list(min=0,max=1),list(min=0,max=1))
  #distr.cop <- mvdc(cop, margins=rep("unif",3), paramMargins = parMarg,marginsIdentical=TRUE)
  #c.X <- rMvdc(mvdc=distr.cop, n=N.set)
  #Y <- punif(c.X)
  #vine.copula<-vine(Y,K=8,base="B-spline",q=2,m=2,pen=1,cores=1,lambda=12500)
}
