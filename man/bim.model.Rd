\name{bim.model}
\alias{bim.model}
\alias{bim.nqtl}
\alias{bim.pattern}
\title{Bayesian model selection for number and pattern of QTL across genome}
\description{
Posterior number and pattern of QTL, along with posterior/prior Bayes
factor ratios.
}
\usage{
bim.model( bim, cross, nqtl = 1, pattern=NULL, exact=FALSE,
  cutoff = 1 )
bim.nqtl( bim )
bim.pattern( bim, cross, nqtl = 1, pattern=NULL, exact=FALSE,
  cutoff = 1 )
}
\arguments{
  \item{bim}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} (extracted by
  \code{bim.cross} if not provided)}
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{cutoff}{percent cutoff for inclusion in model selection}
}
\details{
  \code{bim.model} creates results from both \code{bim.nqtl} and
  \code{bim.pattern}.
  
  \code{bim.nqtl} estimates posterior frequency of number of QTLs as the
  margine over all other model parameters.  However, note that
  posterior may be influenced by prior, while Bayes factor is
  empirically less sensitive for QTL model selection.  Bayes factors are
  ratios of \code{bf=posterior/prior} ratios. 

  \code{bim.pattern} shows at most 15 model patterns with at least
  \code{cutoff} \% posterior are returned.
  Patterns are comma-separate list of chromosomes, with asterisk 
  \code{*} for multiple QTL per chromosome.
  \code{bim} is first subsetted using \code{\link{subset.bim}}.
}
\value{
  List with items \code{nqtl} and \code{pattern}, each containing:
  \item{posterior}{posterior for number of QTL}
  \item{prior}{prior for number of QTL}
  \item{bf}{rank-ordered posterior/prior ratios rescaled so bf[1] = 1}
  \item{bfse}{approximate standard error for bf computed using binomial
    variance}
  In addition, there is an object \code{param} with values for
  \code{nqtl}, \code{pattern}, \code{exact} and \code{cutoff}.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.bim}},\code{\link{plot.bim.model}},\code{\link{subset.bim}}}

\examples{
data( verngeo.bim )
\testonly{library(qtl)}
bim.model( verngeo.bim )
}
\keyword{models}
