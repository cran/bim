\name{fisch}
\alias{fisch}
\alias{fisch.bim}
\non_function{}
\title{Eight QTL Stephens and Fisch simulated data}
\description{
Contains genotypes and phenotypes for data simulated using model in
Stephens and Fisch (1998) but with 90 percent heritability.
}
\usage{
data(fisch)
data(fisch.bim)
}
\format{
  \code{fisch} is f2 (see \code{\link[qtl]{read.cross}} for format).
  \code{fisch.bim} produced by \code{link{run.bmapqtl}}.
}
\seealso{\code{\link[qtl]{read.cross}}, \code{\link{plot.bim}}, \code{\link{run.bmapqtl}}}
\author{ Brian S. Yandell, \url{mailto:yandell@stat.wisc.edu}}
\source{
Patrick J. Gaffney (\url{mailto:paga@lubrizol.com}), Lubrizol Corp.
}
\references{
\url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl}
}
\examples{
data(fisch)
data(fisch.bim)
}
\keyword{datasets}
