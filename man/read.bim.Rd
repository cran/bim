\name{read.bim}
\alias{read.bim}
\title{ Read samples from WinQTL output}
\description{
MCMC samples from WinQTL are organized into an object of class
\code{bim} for diagnostic plots and analysis. This can be used to import
MCMC samples from WinQTL or from an optional file saved by run.bmapqtl.
}
\usage{
read.bim( dir, bimfile, nvalfile = "nval.dat", na.strings="." )
}
\arguments{
  \item{dir}{ directory containing files }
  \item{bimfile}{ name of MCMC sample file (currently with \code{.z} suffix)}
  \item{nvalfile}{parameter file for \code{Bmapqtl} program (currently hard-wired)}
  \item{na.strings}{ missing value code }
}
\value{
  An object of class `bim', which is a list with components
  \item{bmapqtl}{parameter settings for \code{Bmapqtl} run}
  \item{burnin}{burnin sequence data frame, containing
    \code{niter} = iteration number,
    \code{nqtl} = number of QTL,
    \code{LOD} = LOD score based on simulated QTL genotypes,
    \code{mean} = grand mean of phenotype,
    \code{envvar} = environmental variance,
    \code{addvar} = additive variance,
    \code{domvar} = dominance variance (if included), and
    \code{herit} = heritability.}
  \item{iter}{MCMC iteration summary data frame, with same elements as
    \code{burnin}.}
  \item{loci}{data frame with MCMC estimates across loci by iteration,
    containing
    \code{niter} = iteration number,
    \code{nqtl} = number of QTL,
    \code{chrom} = chromosome number,
    \code{locus} = locus position in centi-Morgans,
    \code{add} = additive effect,
    \code{dom} = dominance effect (if included).}
}
}
\author{ Brian S. Yandell, yandell@stat.wisc.edu }

\details{Bmapqtl used to operate as a separate external command, but was
recently incorporated as a C dynamic library. However, \code{read.bim}
is still useful for output from WinQTL. For information on WinQTL visit
\url{http://statgen.ncsu.edu/qtlcart}.}

\seealso{\code{\link{plot.bim}},\code{\link{run.bmapqtl}},\code{\link{read.bmapqtl}}}

\examples{
\dontrun{verngeo.bim <- read.bim( ".", "verngeo.z" )}
}
\keyword{file}
