\name{order.vine,test.ind}
\alias{order.vine}
\alias{order.vine.level}
\title{Ordering the first level of the R-vine.}

\description{
Ordering the levels of the R-vine, depending on the pairwise cAIC- or kendells tau-values.
}

\usage{
order.vine(help.env,test.ind=FALSE)
order.vine.level(tree,help.env)
}

\arguments{
  \item{help.env}{Containing all information, environment of
    vine()}
  \item{test.ind}{Deafult is FALSE, due to (n over 2) possible pair-copulas in the second level of a R-vine with n variables. If TRUE, testing of independence for each pair-copula in the second level of the R-vine is performed.}
  \item{tree}{Forthgoing tree, which used for estimating the next tree, only needed for level 3 and higher levels.} 
}

\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}
\references{
  Nonparametric estimation of simplified vines: comparison of methods, Nagler N., Schellhase, C. and Czado, C. (2017) Dependence Modeling.
 }
\keyword{nonparametric}
