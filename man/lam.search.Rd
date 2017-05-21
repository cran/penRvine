\name{lam.search}
\alias{lam.search}

\title{Search optimal starting vlaue for lambda}

\description{Starts the estimation of the paircopula with different starting values of the penalty parameter 'lambda'.
  
}

\usage{
lam.search(data,K,lam,m,max.iter,q,base,fix.lambda,pen,id)
}

\arguments{
  \item{data}{'data' contains the data. 'data' has to be a matrix or a
    data.frame with two columns.}
  \item{K}{K is the degree of the Bernstein polynomials. In the case of
    linear B-spline basis functions, K+1 nodes are used for the basis functions.}
  \item{lam}{Vector of potential starting values for lambda}
  \item{m}{Indicating the order of differences to be penalised. Default is "m=2".}
  \item{max.iter}{maximum number of iteration, the default is max.iter=51.}
  \item{q}{Order of B-spline basis, i.e. default q=2 means linear B-spline basis.}
  \item{base}{Type of basis function.}
  \item{fix.lambda}{Indicator if fix lambda is used.}
  \item{pen}{Indicating type of penalisation, corresponding to basis type 'basis'.}
  \item{id}{Indification number}
 }

 \details{Fitting of pair-copulas using each potential starting value for the penalty parameter 'lambda'. Returning the pair-copula with lowest cAIC in the fitted list.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{math}

\references{
  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
