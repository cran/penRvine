\name{pen.matrix}
\alias{pen.matrix}

\title{Calculating the penalty matrix P}

\description{
Calculating the penalty matrix P depends on the used basis
function and the selected kind of penalty.
}
\usage{
pen.matrix(penden.env)
}

\arguments{
  \item{penden.env}{Containing all information, environment of
    paircopula().}
}

\details{If paircopula is selected with 'pen=1', the differences of
  order 'm' are penalized. Only for Bernstein polynomials: If paircopula is selected with 'pen=2', the
  integrated squared second order derivatives are used as penalty.
}

\value{
  \item{DDD.sum}{Penalty matrix P}
  Matrix is saved in the environment.
}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\keyword{nonparametric}

