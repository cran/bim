\name{bim.fdr}
\alias{bim.fdr}
\alias{plot.bim.fdr}
\title{Bayesian False Discovery Rate for QTL mapping}
\description{
Bayesian false discover rate (FDR) applied to QTL mapping. Use
plot.bim.fdr to determine FDR for HPD level used in plot.bim.effects.
}
\usage{
bim.fdr(x, cross, nqtl = 1, pattern = NULL, exact = FALSE, chr,
  \dots, levels=<<see below>>,df=3,  qtl = <<see below>>)
plot.bim.fdr(x, cross, \dots, fdr = <<see below>>,
  critical.value=<<see below>>, hpd = NULL )
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or
    chromosome names)} 
  \item{levels}{probability levels for HPD regions (default every 1 percent)}
  \item{df}{degrees of freedom for smoothing spline estimator of prior
    on size}
  \item{qtl}{QTL estimate information from \code{\link{bim.qtl}}}
  \item{fdr}{object returned by call to \code{\link{bim.fdr}}}
  \item{critical.value}{pFDR guidelines (ignored if hpd provided;
    default every .05 to .25)} 
  \item{hpd}{HPD guidelines}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\value{
  \code{bim.fdr} returns an enhanced version of the
  \code{\link{bim.qtl}} object, with additional elements:
  \item{levels}{HPD \code{levels}}
  \item{size}{relative size of HPD region across \code{levels}}
  \item{fdr}{false discover rate across \code{levels}}
  \item{hyp}{estiamtes for probabilities of hypotheses H0: no QTL at
    locus, M0: no QTL anywhere, M1: QTL somewhere}
  \item{prob}{cumulative prior probability}
  \item{spline}{spline estimate of cumulative prior probability}
The \code{size} and \code{fdr} objects are data frames, one element per
chromosome plus the last being \code{all} chromosomes. \code{spline} is
an object of class \code{smooth.spline}.

\code{plot.bim.fdr} returns
  \item{hyp}{estiamtes for probabilities of hypotheses H0: no QTL at
    locus, M0: no QTL anywhere, M1: QTL somewhere}
  \item{fdr}{FDR for corresponding HPD levels (chosen by \code{hpd} or
    \code{critical[value} settings)}
  }
\details{
FDR follows ideas of Storey (\url{http://www.stat.berkeley.edu/~storey})
adapted to Bayesian setting.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.bim.effects}},\code{\link{bim.qtl}}}

\examples{
data( vern )
data( verngeo.bim )
\testonly{library(qtl)}
plot.bim.fdr( verngeo.bim, vern, pattern=c(1,1) )
}
\keyword{models}
\keyword{hplot}
