\name{vernpois.bim}
\alias{vernpois.bim}
\non_function{}
\title{MCMC samples from posterior for 8-week vernalization}
\description{
  Samples use Poisson prior on number of QTL with prior mean 3.
  See \code{read.bim} for details.
}
\usage{data(vernpois.bim)}
\references{ \url{www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\seealso{\code{\link{read.bim}},\code{\link{plot.bim}}}

\examples{
data(vernpois.bim)
plot( vernpois.bim )
}
\keyword{datasets}
