\name{read.bmapqtl}
\alias{read.bmapqtl}
\alias{write.bmapqtl}
\title{Read and write options for WinQTL}

\description{
  write.bmapqtl takes the options and priors in current workspace from
  \code{\link{bmapqtl.options}} and write it to a (user specified) text 
  file. read.bmapqtl opens the file and read the options and priors in.
}

\usage{
write.bmapqtl( dir = ".", nvalfile = "nval.dat" )
read.bmapqtl( dir = ".", nvalfile = "nval.dat" )
}

\arguments{
  \item{dir}{ directory containing files }
  \item{nvalfile}{file name for options and priors.}
}


\value{
  There is no return value from \code{write.bmapqtl}.

  \code{read.bmapqtl} returns a list for options and priors. To replace
  the options and priors in current workspace, use
  \code{.bmapqtl.options <- read.bmapqtl()}.

}

\author{ Brian S. Yandell, yandell@stat.wisc.edu,
  Hao Wu, hao@jax.org}
\details{Bmapqtl used to operate as a separate external command using
  external file \code{nval.dat}, but was
  recently incorporated as a C dynamic library. However,
  \code{read.bmapqtl} and \code{write.bmapqtl} are still useful for
  WinQTL. For information on WinQTL visit
\url{http://statgen.ncsu.edu/qtlcart}.}
\seealso{\code{\link{read.bim}},
  \code{\link{run.bmapqtl}},
  \code{\link{bmapqtl.options}}
}

\examples{
\dontrun{
write.bmapqtl( ".", "nval.dat")
read.bmapqtl(".", "nval.dat" )
}
}

\keyword{file}
