\name{bim.}
\alias{bim.qtl}
\alias{plot.bim.qtl}
\title{Bayesian QTL estimation and mapping of loci}
\description{
  Smooth density, point estimates and HPD regions for QTL loci. Loci are
  estimated at peaks of density.
}
\usage{
bim.qtl(x, cross, nqtl=1, pattern=NULL, exact=FALSE,
  chr, bw=2, levels = seq(.5,.95,by=.05) )
plot.bim.qtl(x, cross, nqtl = 1, pattern = NULL,
  exact = FALSE, chr, bw = 2, qtl=<<see below>>,
  level = 0.80, col = "black", add = FALSE, \dots) 
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or chromosome names)}
  \item{bw}{bandwidth for loci density}
  \item{levels}{saved levels for HPD region}
  \item{qtl}{QTL estimate information (see value)}
  \item{level}{level for HPD region from \code{seq(.5,.95,by=.05)}}
  \item{col}{line color}
  \item{add}{create new plot if FALSE}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\value{
  \code{bim.qtl} returns (and \code{plot.bim.qtl} silently returns) a
  list containing: 
  \item{loci}{data frame with \code{chr} = chromosome name, \code{x} =
    locus on chromosome, \code{y} = height of density}
  \item{dens}{list of density objects by chromosome}
  \item{hpd} {density critical values for high probability density (HPD) region}
}
\details{
Model averaging conditional on at least \code{nqtl} QTL and at least
chromosome \code{pattern} in model.  \code{bim.qtl} does estimation of
density and loci as well as high probability density (HPD) region.
\code{plot.bim.qtl} calls \code{bim.qtl} and plots the density curve.
See \code{\link{plot.bim.effects}} for enhanced plots with HPD and effects
estimates.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.bim.effects}}}

\examples{
data( vern )
data( verngeo.bim )
\testonly{library(qtl)}
plot.bim.qtl( verngeo.bim, vern, 2 )
}
\keyword{models}
\keyword{hplot}
