\name{plot.bim.diag}
\alias{plot.bim.diag}
\title{Marginal and model-conditional summaries of Bayesian interval
  mapping diagnostics}
\description{
A density histogram is drawn for model-averaged summary diagnostics such
as LOD, variance, or heritability.
}
\usage{
plot.bim.diag(x, nqtl=1, pattern=NULL, exact=FALSE,
  items=c("herit", "LOD"), mains=items,
  mfrow=c(length(items), 2), \dots )
}
\arguments{
  \item{x}{ object of class \code{bim} }
  \item{nqtl}{subset on number of QTL}
  \item{pattern}{subset on chromosome pattern of QTL}
  \item{exact}{subset on exact pattern or number of QTL if true}
  \item{items}{diagnostics to be summarized; must be column of data}
  \item{mains}{titles for items}
  \item{mfrow}{plot arrangement parameter for \code{par()}}
  \item{\dots}{graphical parameters can be given as arguments to \code{\link[base]{plot}}}
}
\details{
  Model-averaged density is smooth kernel estimate similar to ordinary
  histogram. A \code{\link[base]{boxplot}} (without outliers) is
  overlaid for comparison   with conditional boxplots. Conditional
  boxplots by number of QTL may   show   indication of model bias for
  small number of QTL. This and   \code{\link{bim.nqtl}} can help
  suggest the minimal model. Diagnostic   items that make sense to plot
  are \code{"LOD"}, \code{"envvar"}   (environmental variance),
  \code{"herit"} (heritability), \code{"mean"} (grand mean),
  \code{"addvar"} (variance of \code{add}), \code{"domvar"} (variance of
  \code{add}). Marginal and conditional medians are printed.
}
\references{ \url{www.stat.wisc.edu/~yandell/qtl/software/Bmapqtl} }
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\seealso{\code{\link{plot.bim}}, \code{\link{density}}, \code{\link[base]{boxplot}}}

\examples{
data( verngeo.bim )
plot.bim.diag( verngeo.bim, 2, items = c("LOD","herit") )
}
\keyword{models}
