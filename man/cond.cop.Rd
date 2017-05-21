\name{cond.cop}
\alias{cond.cop}
\title{Flexible Pair-Copula Estimation in R-vines with Penalized
  Splines}
\description{Calculation of the conditional paircopulas.}
\usage{cond.cop(data,coef,K,diff="u2",Index.basis.D,base,q=2,env,kn1,kn2,int.base1,int.base2,ddb)}
\arguments{
  \item{data}{Considered bivariate data, later used as "u1" and "u2".}
  \item{coef}{Considered coefficients of the splines.}
  \item{K}{Number of marginal knots.}
  \item{diff}{Default="u2", alternatively diff="u1". Determining in which
    direction the pair-copula is differentiated.}
  \item{Index.basis.D}{Vector of indices built in the program before.}
  \item{base}{"B-spline" or "Bernstein".}
  \item{q}{Order of the B-spline, default=2}
  \item{env}{Environment with needed data.}
  \item{kn1}{Marginal knots for the first covariate.}
  \item{kn2}{Marginal knots for the second covariate.}
  \item{int.base1}{Integrated marginal density basis of the first covariate.}
  \item{int.base2}{Integrated marginal density basis of the second covariate.}
  \item{ddb}{Number of coefficients}.
}

\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
