\name{paircopula}
\alias{paircopula}
\title{Flexible Pair-Copula Estimation in R-vines using Bivariate Penalized Splines}
\description{Calculating paircopula with penalized B-splines or penalized
  Bernstein polynomials}
\usage{paircopula(data,K=8,base="Bernstein",max.iter=51,lambda=100,
       data.frame=parent.frame(),m=2,fix.lambda=FALSE,pen=1,q=2,
      length.cond=NULL,eps.we=1e-8)}
\arguments{
  \item{data}{'data' contains the data. 'data' has to be a matrix or a
    data.frame with two columns.}
  \item{K}{K is the degree of the Bernstein polynomials. In the case of
    linear B-spline basis functions, K+1 nodes are used for the basis functions.}
  \item{base}{Type of basis function, default is
    "Bernstein". An alternative is base="B-spline".}
  \item{max.iter}{maximum number of iteration, the default is max.iter=51.}
  \item{lambda}{Starting value for lambda, default is lambda=100.}
  \item{data.frame}{reference to the data. Default reference is the
    parent.frame().}
  \item{m}{Indicating the order of differences to be penalised. Default is "m=2".}
  \item{fix.lambda}{Determining if lambda is fixed or if the
    iteration for an optimal lambda is used, default 'fix.lambda=FALSE'.}
  \item{pen}{'pen' indicates the used penalty. 'pen=1' for the
    difference penalty of m-th order. 'pen=2' is only implemented for
    Bernstein polynomials, it is the penalty based
    on the integrated squared second order derivatives of the Bernstein
    polynomials. Due to numerical difficulties handling the integral of Bernstein polynomials (that is the beta function), this approach works only for K<=15.}
  \item{q}{Order of B-spline basis, i.e. default q=2 means linear B-spline basis.}
  \item{length.cond}{Number of marginal equidistant evaluation points for positivity of copula density, default is number of knots}
  \item{eps.we}{Default eps.we=1e-8}
}
\details{
  Each paircopula is calculated using Bernstein polynomials or B-spline
  densities as basis functions. Optimal coefficients and optimal penalty
  parameter lambda are selected iteratively using quadratic programming.
}
\value{Returning an object of class 'paircopula', consisting of the
  environment 'penden.env', which includes all values.}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
