#####################################################################
##
## bmapqtl.options.R, 08/14/2003, hao@jax.org
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

## default values for options, these are global variables
bmapqtl.options.init <- function() {
.bmapqtl.options <<- NULL
.bmapqtl.options$prior.nqtl <<- "geometric"
.bmapqtl.options$mean.nqtl <<- 3 # prior for number of QTL
.bmapqtl.options$niter <<- 400000
.bmapqtl.options$by <<- 400 # number of iterations, recorded by
.bmapqtl.options$burnin <<- 0.05
.bmapqtl.options$preburn <<- 0.05 # burn-in and pre-burn-in
.bmapqtl.options$nqtl <<- 0 # initial number of QTL
.bmapqtl.options$init <<- c(0.5, -1) # normal(0,.5*s^2)
.bmapqtl.options$prior.mean <<-  c(1,-1) # normal(0,s^2)
.bmapqtl.options$prior.var <<- c(3,-1) # IG(3,s^2)
.bmapqtl.options$prior.add <<- c(0,0) # Beta(2,10)
.bmapqtl.options$prior.dom <<- c(0,0) # Beta(2,10)
.bmapqtl.options$seed <<- 0 # random seed
assign(".bmapqtl.options",.bmapqtl.options,1)
}
bmapqtl.options <- function(...,reset=FALSE)
{
  # take the arguments
  args <- list(...)
  nargu <- length(args)

  # return variable
  if( reset |!exists(".bmapqtl.options")) {
      bmapqtl.options.init()
    result <- .bmapqtl.options
  }
  else
    result <- list(NULL)
  # assign values
  if(nargu&!reset) {
    for (i in 1:nargu) {
      argname <- names(args)[i] # argument name
      argvalue <- args[[i]] # argument value
      if(is.null(argname)) { # trying to get an option
        result[[i]] <- .bmapqtl.options[[argvalue]]
        names(result)[i] <- argvalue
      }
      else {
        # trying to assign an option
        # error checking stuff here
        switch( argname,
               "mean.nqtl" = {
                 if(any(argvalue < 0))
                   stop("Prior for number of QTL need to be greater than or equal to zero")
               },
               "niter" = {
                 if(any(argvalue <= 0))
                   stop("Number of iterations need to be greater than zero")
               },
               "seed" = {
                 if(any(argvalue < 0))
                   stop("Random number seed need to be greater than or equal to zero")
               },
               "prior.nqtl" = {
                 priors = c("geometric","poisson","uniform")
                 argvalue <- priors[ pmatch( tolower( argvalue ), priors, nomatch = 1 ) ]
               }
               )
        # assign values
        .bmapqtl.options[[argname]] <<- argvalue
        result[[i]] <- .bmapqtl.options[[argname]]
        names(result)[i] <- argname
      }
    }
    # if nqtl is bigger than zero, need to add chrom and locus fields
    #if( .bmapqtl.options$nqtl > 0 ) {
    #  .bmapqtl.options$chrom <<- rep(1, .bmapqtl.options$nqtl)
    #  .bmapqtl.options$locus <<- rep(1, .bmapqtl.options$nqtl)
    #}
  }
  cat( "simulate", as.integer( .bmapqtl.options$niter ), "MCMC steps, recording by",
      as.integer( .bmapqtl.options$by ), "with", .bmapqtl.options$burnin,
      "burnin and", .bmapqtl.options$preburn, "pre-burnin\n" )
  cat( paste( "prior for number of QTL: ", .bmapqtl.options$prior.nqtl, "(",
             .bmapqtl.options$mean.nqtl, ")\n", sep = "" ))
  cat( "initial number of QTL:", .bmapqtl.options$nqtl, "\n" )
  cat( "hyperparameters for priors:\n" )
  print(t(data.frame(.bmapqtl.options[8:12])))
  cat( "random seed:", .bmapqtl.options$seed, "\n" )
  invisible(result)
}

