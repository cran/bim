\name{plot.bim.loci}
\alias{plot.bim.loci}
\title{Jittered plot of Bayesian QTL loci samples by chromosome}
\description{
Each point is one locus from the Bayesian QTL estimates, plotted
vertically by chromosome, jittered to give a sense of density.
}
\usage{
plot.bim.loci(x, cross, nqtl=1, pattern=NULL, exact=FALSE,
  chr, labels=TRUE, amount=.35, cex, \dots )
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} (extracted by
  \code{bim.cross} if not provided) }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or chromosome names)}
  \item{labels}{include marker labels if TRUE}
  \item{amount}{amount of jitter (between 0 and .45)}
  \item{cex}{ character expansion (may be invisible if too
    small--default set by \code{bim.cex})}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\details{
  Focuses attention on chromosome lengths and concentration of QTL loci
  estimates. Horizontal blue lines at markers if \code{cross}
  included. Adjust \code{amount} and \code{cex} to modify look of 
  plot. Most useful when looking at multiple chromosomes.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link[base]{jitter}},\code{\link{subset.bim}}}

\examples{
data( vern )
data( verngeo.bim )
\testonly{library(qtl)}
plot.bim.loci( verngeo.bim, vern, 2 )
}
\keyword{models}
