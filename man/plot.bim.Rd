\name{plot.bim}
\alias{plot.bim}
\title{Diagnostics plots for Bayesian interval mapping}
\description{
  Diagnostic plots highlight putative QTL loci and effects as well as
  providing graphical model assessment tools.
}
\usage{
plot.bim(x, cross, nqtl=1, pattern=NULL, exact=FALSE, ... )
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{...}{graphical parameters can be given as arguments to \code{\link[base]{plot}}}
}
\details{
  Import object \code{cross} using \code{read.cross} and simulate
  Bayesian interval mapping sample object \code{x} using
  \code{run.bmapqtl}, respectively. This routine creates five plots via
  calls to several other plot routines for bim data.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[base]{plot}},\code{\link[qtl]{read.cross}},\code{\link{run.bmapqtl}},
  \code{\link{plot.bim.mcmc}}, \code{\link{plot.bim.loci}},
  \code{\link{plot.bim.model}}, \code{\link{plot.bim.effects}},  \code{\link{plot.bim.diag}} }

\examples{
data( vern )
data( verngeo.bim )
plot.bim( verngeo.bim, vern )
}
\keyword{hplot}
\keyword{models}
