\name{run.bmapqtl}
\alias{run.bmapqtl}
\title{Run Bmapqtl reversible jump MCMC}
\description{Run Bmapqtl reversible jump MCMC
}

\usage{
run.bmapqtl(cross, pheno=1, chrom=0, result.file="")
}

\arguments{
  \item{cross}{An object of class "cross"}
  \item{pheno}{phenotypic trait number or name}
  \item{chrom}{chromosome number (0 = all). Note that at this stage,
    it can either be all chromosomes (0) or a single chromosome.}
  \item{result.file}{Result file name. If this is specified, the result
    will be written to a tab-delimited text file. If it is an empty string
    (""), result will not be output as a file.}
}

\value{
  The result is an object of class "bim", which is a list with the
  following components:
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

\author{Hao Wu, hao@jax.org}

\seealso{\code{\link{plot.bim}},
  \code{\link{read.bim}}
}

\examples{
\dontrun{
data(vern)
bim = run.bmapqtl(vern)
plot(bim)
}
}

\keyword{models}

