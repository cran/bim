\name{bim.effects}
\alias{bim.effects}
\alias{plot.bim.effects}
\title{Bayesian QTL map of loci and effects}
\description{
  Histogram of QTL loci plus scatter plots and smooth estimates of QTL
  effects (additive and dominance--Mather-Jinks model).
}
\usage{
plot.bim.effects(x, cross, nqtl = 1, pattern = NULL, exact = FALSE,
  chr, bw = 2, qtl = <<see below>>, cex, level = .80,
  project=<<see below>>, main=<<see below>>, mfcol=<<see below>>, \dots )
bim.effects(x, cross, nqtl = 1, pattern = NULL, exact = FALSE, chr,
  bw = 2, qtl = <<see below>>)

}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{cross}{ corresponding object of class \code{cross} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{chr}{subset of chromosomes to plot (numerical indices or
    chromosome names)} 
  \item{bw}{bandwidth for loci density}
  \item{qtl}{QTL estimate information (see value)}
  \item{cex}{ character expansion for scatter plots }
  \item{prefix}{ prefix for chromosome names, typically \code{c} }
  \item{level}{level for HPD region from \code{seq(.5,.95,by=.05)}}
  \item{project}{project identifier, default is name of \code{x} argument}
  \item{main}{title, with default \code{expression(project summaries 
      with m >= nqtl)} and substitutions for \code{project} and 
    \code{nqtl}} 
  \item{mfcol}{parameter to \code{par} set to 2 or 3 depending on if
    dominance is present}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\value{
  \code{bim.effects} returns an enhanced version of the
  \code{\link{bim.qtl}} object, with additional elements
  \item{add}{additive smoothing spline by chromosome}
  \item{dom}{dominance smoothing spline by chromosome (if present)}
  \item{est}{estimates of QTL}
  The first two elements are lists of
  \code{\link[modreg]{smooth.spline}} objects, one for each chromosome.
  The third element is a data frame with the following elements, with
  one entry per QTL. The final entry is for the mean.
  \item{chrom}{chromosome name}
  \item{loci}{loci position along chromosome}
  \item{add}{additive effect}
  \item{add.sd}{SD for additive effect}
  \item{dom}{dominance effect (if present)}
  \item{dom.sd}{SD for dominance effect}
  }
\details{
Model averaging conditional on at least \code{nqtl} and at least
chromosome \code{pattern} QTL in model.
\code{bim.effects} uses \code{\link{bim.qtl}} to estimate the QTL 
First row has histogram of loci overlaid with smoothed density estimate
(blue line) and HPD regions (red circles along origin).
Second row is additive effect; third if present has dominance effect.
Smoothing spline fit along scatter plot to estimate mean effects (solid
blue line), plus or minus two SDs (dashed blue line; estimated as well
by smoothing spline). Simple linear regression used if 50 or fewer
samples for a chromosome.
Purple triangles indicate marker locations if supplied. Red circles and
vertical dashed lines at estimated loci and effects. Vertical black
lines delineate chromosomes.
}
\references{ \url{http://www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.bim}},\code{\link{subset.bim}},\code{\link{bim.qtl}}}

\examples{
data( vern )
data( verngeo.bim )
\dontrun{
vern.qtl <- bim.effects( verngeo.bim, vern, 2 )
plot.bim.effects( verngeo.bim, vern, 2, qtl = vern.qtl )
}
}
\keyword{models}
\keyword{hplot}
