\name{plot.bim.mcmc}
\alias{plot.bim.mcmc}
\title{Bayesian MCMC sequence plots for burnin and iterations.}
\description{
  Plot MCMC burnin and iteration sequences for
  number of QTL, LOD score, environmental SD, and heritability.
  
}
\usage{
plot.bim.mcmc(x, element = c("burnin","iter"),
  xlab = c("burnin sequence","mcmc sequence"),
  items = <<see below>>, ylabs = items,
  types = <<see below>>, \dots )

}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{element}{ elements of \code{x} to plot}
  \item{xlab}{x-axis labels for elements}
  \item{items}{diagnostic item names to be plotted as time series
    (default is all names in \code{x$iter}}
  \item{ylabs}{labels for items}
  \item{types}{line types (default is both for first item--usually
    number of QTL--and lines for the rest)}
  \item{\dots}{graphical parameters can be given as arguments to
    \code{\link[base]{plot}}} 
}
\details{
  Plots show grey line for saved sequences, usually nearly vertical due
  to mixing, plus a smoothed curve (blue).
Burnin plots should show initial jump in number of QTL during pre-burnin
sequence, ideally settling into a well-behaved Markov chain by the end
of the burn-in. That is, look for visual evidence of stationarity as
flat blue line and even spread of grey. Trend or cyclic pattern
indicates Markov chain is not mixing well and needs to be run much
longer. No formal assessment is provided.
}
\references{ \url{www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{run.bmapqtl}}}

\examples{
data( verngeo.bim )
plot.bim.mcmc( verngeo.bim )
}
\keyword{models}
