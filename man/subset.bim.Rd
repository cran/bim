\name{subset.bim}
\alias{subset.bim}
\title{Subsetting Bayesian interval mapping data}
\description{
Subset Bayesian interval mapping iterations on number of QTL and/or
chromosome pattern of QTL, using exact match or inclusive subsetting.
}
\usage{
subset.bim(x, cross, nqtl=1, pattern=NULL, exact=FALSE, chr, \dots )
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or
    chromosome names)}
  \item{\dots}{additional arguments to \code{\link[base]{subset}}}
}
\details{
Subset to include only iterations with at least \code{nqtl} number of
QTL and at least the \code{pattern} across chromosomes. \code{pattern}
  is a vector of chromosome identices, with repeats for
  multiple linked QTL on a chromosome. If \code{exact=FALSE}, then all
  iterations with at least the given \code{pattern} and \code{nqtl} are
  included. \code{nqtl} will be reset to \code{length(pattern)} if it is
  smaller than this value. Note that \code{pattern} should be number
  codes corresponding to those used in the \code{x} object. At
  present, chromosome names are not allowed. Further subsets to only
  include QTL from these iterations that are on chromosomes \code{chr}.
}
\references{ \url{www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{read.bim}}}

\examples{
\dontrun{
bim223 <- subset.bim( bim, pattern = c(2,2,3) )}
}
\keyword{utilities}
