\name{plot.bim.model}
\alias{plot.bim.model}
\title{Graphical model assessment for Bayesian interval mapping}
\description{
  Model-averaged posteriors and posterior/prior ratios for graphical Bayes
  factor assessment. First row concerns number of QTL, second row
  evaluates pattern of QTL across chromosomes.
}
\usage{
plot.bim.model( x, cross = NULL, nqtl = 1, pattern=NULL, exact=FALSE,
  cutoff = 1, assess = <<see below>>, \dots )
}
\arguments{
  \item{x}{ object of class \code{bim} or class \code{bim.model} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{cutoff}{percent cutoff for inclusion in model selection}
  \item{assess}{object of class \code{bim.model} from \code{\link{bim.model}}}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\details{
  \code{plot.bim.model} uses results (\code{assess}) of \code{bim.model}
  and arranges plots on a single page.  If \code{x} is of class
  \code{bim.model}, then \code{assess} is set to its value and the other
  arguments are ignored. Left plot is of posterior against  model
  identifier, while right plot assesses Bayes factors.  Since Bayes
  factors are ratios of \code{posterior/prior} ratios, a semi-log plot
  of \code{posterior/prior} against model identifier (\emph{m} = number
  of QTL or \emph{M} = model pattern) provides a graphical model
  assessment tool with a BF threshold yardstick.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{bim.model}}}

\examples{
data( vern )
data( verngeo.bim )
\testonly{library(qtl)}
plot.bim.model( verngeo.bim, vern, 2 )
}
\keyword{models}
