\name{Derv2}
\alias{Derv2}
\title{Calculating the second order derivative of the paircopula likelihood function w.r.t. parameter b}
\description{
 Calculating the second order derivative of the paircopula likelihood function
 w.r.t. parameter v. 
}
\usage{
Derv2(penden.env,temp=FALSE,lambda=NULL,lam.fit=FALSE)
}

\arguments{
	\item{penden.env}{Containing all information, environment of
	  paircopula().}
	\item{temp}{Default=FALSE,if TRUE temporary calculations of
	  optimal parameters are done.}
	\item{lambda}{Default=NULL, i.e. the saved smoothing parameter lambda in the
	  environment is used. Alternatively, temporary values of lambda
	  are used for optimization of lambda.}
        \item{lam.fit}{Default=FALSE, indicating if the first derivative is calculated to determine the next optimal penalty parameter lambda.}
      }
\details{
  We approximate the second order derivative in this approach with the negative fisher information. 
}

\value{
  \item{Derv2.pen}{second order derivative w.r.t. v with penalty}
  \item{Derv2.cal}{second order derivative w.r.t. v without
    penalty. Needed for calculating of e.g. AIC.}
  Derv2.cal and Derv2.pen are saved in the environment.
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
\keyword{math}
