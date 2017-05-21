\name{independ.test}
\alias{independ.test}

\title{Testing for independence between two margins of pair-copulas}

\description{
Two margins of a pair-copulas are tested for independence in front of the estimation of the pair-copula, following Genest, C. and Remillard, B. (2004).
}
\usage{
independ.test(data)
}

\arguments{
  \item{data}{Bivariate data for testing of independence.}
}

\details{
  The test procedure is from the 'copula'-package.
}

\value{
  The p-value of the test statistic is returned.
}
\author{Christian Schellhase <cschellhase@wiwi.uni-bielefeld.de>}

\references{
     Genest, C. and Remillard, B. (2004) Tests of independence and
     randomness based on the empirical copula process.  _Test_ *13*,
     335-369.
}
