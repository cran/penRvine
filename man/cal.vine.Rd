\name{cal.vine}
\alias{cal.vine}
\title{Flexible Pair-Copula Estimation in R-vines with Penalized
  Splines}
\description{Calculating the density of the estimated R-vine at the
  point(s) val.}
\usage{cal.vine(obj,val,cores)}
\arguments{
  \item{obj}{Vine object of class 'penRvine'.}
  \item{val}{Values in which the current R-vine should be evaluated.}
  \item{cores}{Default=NULL, the number of cpu cores used for parallel computing can be specified.}

}
\details{The current R-vine is evaluated in val and the corresponding
  density values are returned.}
\value{The returing values are the density of the current R-vine at the
  point(s) 'val'.}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{Flexible Pair-Copula Estimation in D-vines using Bivariate Penalized
  Splines, Kauermann, G. and Schellhase, C. (2014), Statistics and Computing 24(6): 1081-1100).

  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
