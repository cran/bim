\name{bmapqtl.options}
\alias{bmapqtl.options}
\title{Options Settings for BmapQTL}
\description{This function allows the users to set and examine the
  options and priors for reversible jump MCMC QTL model search. The
  options are stored in global hidden variable .bmapqtl.options.}

\usage{
bmapqtl.options(\dots,reset=FALSE)
}

\arguments{
  \item{reset}{reset options to default settings}
  \item{\dots}{any options can be defined, using `name = value'. The
    available options are listed below.
    \describe{
      \item{prior.nqtl}{one of \code{"poisson"}, \code{"geometric"} or
	\code{"uniform"} (partial match OK).  Default is
	\code{"geometric"}.}
      \item{mean.nqtl}{prior for number of QTL (default = 3).}
      \item{niter}{number of iterations (default = 400000).}
      \item{by}{record every \code{by} iterations (default = 400).}
      \item{burnin}{proportion of iterations for burnin (default = .05).}
      \item{preburn}{proportion of burnin used for preburn (default = .05).}
      \item{nqtl}{initial number of QTL (default = 0).}
      \item{init}{2-vector with initial mean and variance (see details).}
      \item{prior.mean}{2-vector with prior for mean (see details).}
      \item{prior.var}{2-vector with prior for variance (see details).}
      \item{prior.add}{2-vector with prior for additive effects (see details).}
      \item{prior.dom}{2-vector with prior for dominance effects (see
	details).}
      \item{seed}{Random number seed. If 0(default), seed will be generated
	randomly.}
    }
  }
}

\details{
  Initial mean and variance \code{init} are used directly if the
  variance is positive. For values \code{c(a,-b)} use mean 0 and
  variance \eqn{.5s^2}.

  Prior arguments have coded meaning as follows, with \code{a} and
  \code{b} non-negative. For the grand mean, \code{c(a,b)} signifies
  \code{Normal(a,b)}, while \code{c(a,-b)} stands for
  \code{Normal}(\eqn{\bar{Y},bS^2}), with \eqn{\bar{Y}} and
  \eqn{S^2} the phenotype mean and variance, respectively.  Additive and
  dominance priors coded  as \code{c(a,b)} have \code{Normal(a,b)}; as
  \code{c(a,-b)} have \code{Normal}(\eqn{0,aS^2}); as \code{c(-a,-b)} have
  \code{Normal}(\eqn{0,\beta S^2}) with \eqn{\beta\sim}\code{Beta(2,10)}
  if \code{a=b=0} or \eqn{\beta\sim}\code{Beta(a,b)} otherwise. Variance
  has prior inverse Gamma: \code{IG(a,b)} for \code{c(a,b)}, or
  \code{IG(a,\eqn{s^2})} for \code{c(a,-b)}.
}

\value{
  A list (in any case) with the previous values of the options
  changed, or all options when no arguments were given.
}

\author{Hao Wu, hao@jax.org}

\seealso{\code{\link[base]{options}}}

\examples{
bmapqtl.options()
bmapqtl.options(niter=100000, by=1000)
}
    
\keyword{utilities}
