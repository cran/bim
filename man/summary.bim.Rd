\name{summary.bim}
\alias{summary.bim}
\alias{summary.bim.model}
\alias{summary.bim.qtl}
\title{Summary of Bayesian interval mapping samples}
\description{
Quick summary of Bayesian interval mapping samples.
}
\usage{
summary.bim( object, \dots )
summary.bim.model( object, \dots )
summary.bim.qtl( object, \dots )
}
\arguments{
  \item{object}{ object of class \code{bim} (from \code{\link{run.bmapqtl}}),
    \code{bim.model} (from \code{\link{bim.model}}), or \code{bim.qtl} (from
    \code{\link{bim.qtl}} or \code{\link{bim.effects}})}
  \item{\dots}{additional arguments to \code{\link[base]{summary}}}
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{run.bmapqtl}},\code{\link{bim.model}},\code{\link{bim.qtl}},
  \code{\link{bim.effects}}}

\examples{
data( verngeo.bim )
\testonly{library(qtl)}
summary( verngeo.bim )
verngeo.model <- bim.model( verngeo.bim )
summary( verngeo.model )
## estimate QTL density and locate peak(s)
verngeo.qtl <- bim.qtl( verngeo.bim )
## augment bim.qtl with estimates of effects
verngeo.qtl <- bim.effects( verngeo.bim, qtl = verngeo.qtl )
## same idea but with just one call
verngeo.qtl <- bim.effects( verngeo.bim )
}
\keyword{utilities}
