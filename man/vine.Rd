\name{vine}
\alias{vine}
\title{Flexible Pair-Copula Estimation in vines with Penalized Splines}
\description{Calculating vines with penalized B-splines or penalized
  Bernstein polynomials}
\usage{vine(data,K=8,lambda=100,pen=1,base="B-spline",m=2,cores=NULL,
q=2,test.ind=FALSE,test.ind1=FALSE,selec="cAIC",max.iter=51,RVM=NULL,
lam.vec=NULL,l.search=FALSE,fix.lambda=FALSE,id=NULL)}
     
\arguments{
  \item{data}{'data' contains the data. 'data' has to be a matrix or a
    data.frame with two columns.}
  \item{K}{K is the degree of the Bernstein polynomials. In the case of
    linear B-spline basis functions, K+1 nodes are used for the basis functions.}
  \item{lambda}{Starting value for lambda, default is lambda=100.}
  \item{pen}{'pen' indicates the used penalty. 'pen=1' for the
    difference penalty of m-th order. 'pen=2' is only implemented for
    Bernstein polynomials, it is the penalty based
    on the integrated squared second order derivatives of the Bernstein polynomials.}
  \item{base}{Type of basis function, default is
    "B-spline". An alternative is base="Bernstein".}
  \item{m}{Indicating the order of differences to be penalised. Default is "m=2".}
  \item{cores}{Default=NULL, the number of cpu cores used for parallel computing can be specified.}
  \item{q}{Order of B-splines. Default is q=2, NULL if Bernstein polynomials are used.}
  \item{test.ind}{Default=FALSE, if TRUE each pair-copula starting in level 3 of the R-vine is tested for independence.}
  \item{test.ind1}{Default=FALSE, if TRUE each pair-copula in level 2 of the R-vine is tested for independence.}
  \item{selec}{Default="cAIC", determines the selection criteria of the vine structure. selec="ken.tau" chooses Kendells tau for selection.}
  \item{max.iter}{maximum number of iteration, the default is max.iter=51.}
  \item{RVM}{Default=NULL, RVM is RVine-Matrix determining the structure of the vine.}
  \item{l.search}{Default=FALSE, indicating if a search about several starting values for lambda should be performed. If search is selected, the starting value 'lambda' does not work anymore.}
  \item{lam.vec}{Vector of candidate values for penalty parameter lambda}.
  \item{fix.lambda}{Default=FALSE, indicating if lambda is fixed or not.}
  \item{id}{Indification number}
}
\details{The calculation of the vine is done stepwise. From the second level, each paircopula is calculated
  (parallel or not) until the highest level. The specifications in 'vine' are done for every paircopula in the
  vine. There is no option to change parameters for some pair-copulas.
}
\value{Returning a list containing
  \item{vine}{an object of class 'penVine'}
  \item{log.like}{the estimated log-likelihood}
  \item{AIC}{AIC value}
  \item{cAIC}{corrected AIC value}
  \item{K}{Number of K}
  \item{order}{the used order of the first level}
  \item{S}{Sequence seq(1:(dim(data)[2]))}
  \item{N}{Number of observations, that is dim(data)[1]}
  \item{base}{Used basis function}
  }

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
