#####################################################################
##
## $Id: bmapqtl.R,v 1.0 2002/07/04 yandell@stat.wisc.edu Exp $
##
##     Copyright (C) 2002 Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
##############################################################################


read.bmapqtl <- function( dir = ".",
                         nvalfile = "nval.dat" )
{
  cat( "Bmapqtl parameter file is", nvalfile, "in directory", dir,  "\n" )
  if (!missing(dir)) {
    n <- nchar(dir)
    if (substring(dir, n, n) == "/") 
      dir <- substring(dir, 0, n - 1)
    nvalfile <- file.path(dir, nvalfile )
  }
  nvals <- scan( nvalfile, 0, nmax = 6, quiet = TRUE )
  nval1 <- nvals[1]
  prior <- c("poisson","geometric","uniform")
  if( nval1 < 1 | nval1 > 3 )
    stop( paste( "can only handle these priors:", paste( prior, collapse = ", " )))
  nval1 <- prior[nval1]
  name.vals <- c( nval1,"burnin","preburn","niter","by","nqtl")
  nqtl <- nvals[6]
  
  nvals <- scan( nvalfile, 0, nmax = 6 + 2 * nqtl + 11, quiet = TRUE )
  
  if( nqtl > 0 ) {
    name.vals <- c( name.vals,
                   paste( "chrom", seq( nqtl ), sep = "" ),
                   paste( "locus", seq( nqtl ), sep = "" ))
  }
  bmapqtl <- list( prior.nqtl = nval1, mean.nqtl = nvals[17+2*nqtl],
                  niter = nvals[4], by = nvals[5],
                  burnin = nvals[2], preburn = nvals[3],
                  nqtl = nvals[6] )
  if( nqtl > 0 ) {
    bmapqtl$chrom <- nvals[6+seq(nqtl)]
    bmapqtl$locus <- nvals[6+nqtl+seq(nqtl)]
  }
  bmapqtl$init <- nvals[6+2*nqtl+1:2]
  bmapqtl$prior.mean <- nvals[8+2*nqtl+1:2]
  bmapqtl$prior.var <- nvals[10+2*nqtl+1:2]
  bmapqtl$prior.add <- nvals[12+2*nqtl+1:2]
  bmapqtl$prior.dom <- nvals[14+2*nqtl+1:2]
  bmapqtl$runfile <- nvalfile
  for( i in names(bmapqtl))
    .bmapqtl.options[[i]] <<- bmapqtl[[i]]
  bmapqtl
}


##############################################################################
write.bmapqtl <- function( dir = ".", nvalfile = "nval.dat" )
{
  # get variables from options
  prior.nqtl <- .bmapqtl.options$prior.nqtl
  mean.nqtl <- .bmapqtl.options$mean.nqtl
  mean.nqtl <- .bmapqtl.options$mean.nqtl
  niter <- .bmapqtl.options$niter
  by <- .bmapqtl.options$by
  burnin <- .bmapqtl.options$burnin
  preburn <- .bmapqtl.options$preburn
  nqtl <- .bmapqtl.options$nqtl
  init <- .bmapqtl.options$init
  prior.mean <- .bmapqtl.options$prior.mean
  prior.var <- .bmapqtl.options$prior.var
  prior.add <- .bmapqtl.options$prior.add
  prior.dom <- .bmapqtl.options$prior.dom
  # make chrom and locus from nqtl
  chrom <- rep( 1, nqtl )
  locus <- rep( 1, nqtl )
  
  cat( "Creating Bmapqtl parameter file", nvalfile, "in directory", dir,  "\n" )
  if (!missing(dir)) {
    n <- nchar(dir)
    if (substring(dir, n, n) == "/") 
      dir <- substring(dir, 0, n - 1)
    nvalfile <- file.path(dir, nvalfile )
  }
  if( file.exists( nvalfile )) {
    warning( paste( "previous file", nvalfile, "moved to *.mov" ))
    file.rename( nvalfile, paste( nvalfile, "mov", sep = "." ))
  }
  ## choice of prior for number of QTL
  priors <- c("poisson","geometric","uniform","exponential")
  prior.nqtl <- pmatch( tolower( prior.nqtl ), priors )
  if( is.na( prior.nqtl ))
    stop( paste( "prior must be one of", paste( priors[1:3], collapse = ", " )))
  if( prior.nqtl == 4 )
    prior.nqtl <- 2
  ## burnin and preburn
  write( prior.nqtl, nvalfile, append=FALSE)
  write( paste( burnin, preburn ), nvalfile, append=TRUE)
  write( paste( as.integer( c( niter, by )), collapse = " " ), nvalfile, append=TRUE)
  write( nqtl, nvalfile, append=TRUE)
  if( nqtl > 0 ) {
    write( paste( chrom, collapse = " " ), nvalfile, append=TRUE)
    write( paste( locus, collapse = " " ), nvalfile, append=TRUE)
  }
  write( paste( init, collapse = " " ), nvalfile, append=TRUE)
  write( paste( prior.mean, collapse = " " ), nvalfile, append=TRUE)
  write( paste( prior.var, collapse = " " ), nvalfile, append=TRUE)
  write( paste( prior.add, collapse = " " ), nvalfile, append=TRUE)
  write( paste( prior.dom, collapse = " " ), nvalfile, append=TRUE)
  write( mean.nqtl, nvalfile, append = TRUE )

  read.bmapqtl( nvalfile = nvalfile )
}
