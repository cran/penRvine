\name{marg.likelihood}
\alias{marg.likelihood}

\title{Calculating the marginal likelihood}

\description{
Calculating the marginal likelihood of paircopula().
}
\usage{
marg.likelihood(penden.env,pen.likelihood,temp=FALSE)
}

\arguments{
  \item{penden.env}{Containing all information, environment of
    paircopula().}
  \item{pen.likelihood}{Actual penalized likelihood for calculation,
    temporary or not.}
  \item{temp}{Default=FALSE, indicating if temporary values throughout
  iteration are calculated.}
}

\value{
  \item{marg.log.like}{Marginal log-likelihood, saved in the environment}
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
\keyword{nonparametric}

