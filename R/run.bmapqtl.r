#####################################################################
##
## $Id: run.bmapqtl.r,v 1.1 2004/04/30 14:04:19 jgentry Exp $
## run.bmapqtl.R, 08/14/2003 hao@jax.org
##
## Part of the R/bim package
##
## run.bmapqtl calls MCMC simulation
##
##     Copyright (C) 2003 Hao Wu, The Jackson Lab, & Brian S. Yandell, UW-Madison
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

run.bmapqtl <- function(cross, pheno=1, chrom=0, result.file="")

{
  if( !exists( ".bmapqtl.options" ))
    bmapqtl.options()
  ## error checking stuff
  if (class(cross)[2] != "cross")
    stop("The first input variable is not an object of class cross.")
  if(is.character(pheno))
    pheno <- pmatch(pheno,names(cross$pheno),nomatch=0)
  if(pheno <= 0)
    stop("Phenotype column number need to be positive")
  if(chrom < 0)
    stop("Chromosome number cannot be negative")

  # prepare data
  # number of chromosomes
  if(chrom == 0)
    chrom <- 1:nchr(cross)
  
  # map data
  chr.array <- NULL
  nmar.array <- NULL
  pos.array <- NULL
  dist.array <- NULL
  # loop thru chromosomes
  for(i in chrom) {
    # marker position
    pos.i <- cross$geno[[i]]$map
    pos.array <- c(pos.array, pos.i)
    # chromosome
    chr.i <- rep(i, length(pos.i))
    chr.array <- c(chr.array, chr.i)
    # marker order
    nmar.i <- 1:length(pos.i)
    nmar.array <- c(nmar.array, nmar.i)
    # distance
    dist.i <- c( (pos.i[-1]-pos.i[-max(nmar.i)]), 0)
    dist.array <- c(dist.array, dist.i)
  }

  ## find the index for individuals with missing phenotypes
  ## and exclude them later
  idx.missing <- which(is.na(cross$pheno[,pheno]))
  
  ## there are more cross types in QTLCart than in R/qtl
  ## need to extend that in R/qtl (add one more class?)
  cross.str <- class(cross)[1]

  ## make genotype matrix
  genodata <- NULL
  for(i in chrom) {
    if(length(idx.missing) != 0)
      tmp <- as.matrix(cross$geno[[i]]$data[-idx.missing,])
    else
      tmp <- as.matrix(cross$geno[[i]]$data)
    
    ## change the coding to:
    ## AA -> 1, AB -> 0, BB -> -1, missing -> -3
    ## not AA -> -2, not BB -> 2
    ## but R/qtl's riself and risib changed BB to AB, so change back
    
    #tmp[is.na(tmp)] <- -3 # missing
    #tmp[tmp==2] <- 0 # AB
    #tmp[tmp==3] <- -1 #BB
    #tmp[tmp==4] <- 2 # not BB
    #tmp[tmp==5] <- -2 # not AA
    #tmp[is.na(tmp)] <- -3 # missing

    tmp[is.na(tmp)] <- -3 # missing
### tmp[tmp==1] <- 1 # AA
    tmp[tmp==2] <- ifelse(cross.str=="riself"|cross.str=="risib",1,0) # 0 AB
    tmp[tmp==3] <- -1 #BB
    tmp[tmp==4] <- 2 # not BB
    tmp[tmp==5] <- -2 # not AA

    ## bind to genodata
    genodata <- cbind(genodata, tmp)
  }
  
  ## make a seed if the input one is zero
  if(.bmapqtl.options$seed == 0) {
    if(exists(".Random.seed",1))
      rm(.Random.seed,pos=1)
    runif(1)
    seed <- .Random.seed[2]
  }
  else 
    seed <- .bmapqtl.options$seed

  ## make cross type parameters
  ## translation of R/qtl cross types
  cross.str <- switch( cross.str,
                      bc = "B1",
                      f2 = "RF2",
                      riself = "RI1",
                      risib = "RI2",
                      cross.str)

  switch( cross.str,
         "B1"= {
           crosstype.1 <- 1
           crosstype.2 <- 1
         },
         "B2"= {
           crosstype.1 <- 2
           crosstype.2 <- 1
         },
         "RF2"= {
           crosstype.1 <- 4
           crosstype.2 <- 2
         },
         { switch(substr(cross.str,1,2),
                  "SF"= {
                    crosstype.1 <- 3
                    crosstype.2 <- as.integer(substr(cross.str,3,10000))
                  },
                  "RI"= {
                    crosstype.1 <- 5
                    crosstype.2 <- as.integer(substr(cross.str,3,10000))
                  },
                  stop(paste("cross type", cross.str,"not recognized"))
                  )
         })
  
### should check that cross type is consistent with data
### otherwise call below will give many "denom = 0.0 in cond_prob..." messages
  
  ## call engine function
  nind <- nind(cross)-length(idx.missing)
  nchr <- length(chrom)
  nmark <- nmar(cross)[chrom]
  if(length(idx.missing) != 0)
    y <- cross$pheno[-idx.missing,pheno]
  else
    y <- cross$pheno[,pheno]

  # size for return variable
  nret <- ceiling((.bmapqtl.options$niter * (1+.bmapqtl.options$burnin)) /
    .bmapqtl.options$by * 30)

  cat("Bayesian interval mapping MCMC run in progress.",
      "\nCount of 1000 iterations shown separated by dots (negative for burnin):\n")
  z <- .C("R_mcmc",
          # input variables
          as.integer(nind), # number of individuals
          as.integer(nchr), # number of chromosomes
          # the following four items are for genetic map
          as.integer(chr.array),
          as.integer(nmar.array),
          as.double(pos.array),
          as.double(dist.array),
          # genotype data matrix
          as.integer(genodata),
          # number of markers per chromosome
          as.integer(nmark),
          # selected phenotype values
          as.double(y),
          # parameters for cross type
          as.integer(crosstype.1),
          as.integer(crosstype.2),
          # the following are parameters
          as.double(.bmapqtl.options$burnin), # burnin
          as.double(.bmapqtl.options$preburn), # pre burnin
          as.integer(.bmapqtl.options$niter), # iterations
          as.integer(.bmapqtl.options$by), # increment
          as.integer(seed), # random seed
          # the following are priors
          as.double(.bmapqtl.options$prior.add[1]), # mean add
          as.double(.bmapqtl.options$prior.add[2]), # var add
          as.double(.bmapqtl.options$prior.dom[1]), # mean dom
          as.double(.bmapqtl.options$prior.dom[2]), # var dom
          as.double(.bmapqtl.options$prior.mean[1]), # mean of mean
          as.double(.bmapqtl.options$prior.mean[2]), # var of mean
          as.double(.bmapqtl.options$prior.var[1]), # mean of var
          as.double(.bmapqtl.options$prior.var[2]), # var of var
          as.character(.bmapqtl.options$prior.nqtl), # prior type of QTL
          # parameter for QTL, mean for poisson, range for uniform
          as.double(.bmapqtl.options$mean.nqtl),
          # initial values for chain
          as.double(.bmapqtl.options$init[1]), # init mu
          as.double(.bmapqtl.options$init[2]), # init s2
   
          # return variables
          result=as.double(rep(0, 13*nret)),
          PACKAGE="bim"
          )

  tbl <- as.data.frame(matrix(z$result, ncol=13, byrow=T))
  # get rid of zeros
  idx <- apply(tbl, 1, function(x) all(x==0))
  tmp <- min(which(idx))
  tbl <- tbl[1:(tmp-1),]
  colnames(tbl) <- c("niter", "nqtl", "iqtl", "chrom",
                     "LOD", "mu", "sigmasq", "addvar",
                     "domvar", "add", "dom", "locus", "esth")

  # write result to a file if specified
  if(nchar(result.file) != 0)
    write.table(tbl, file=result.file, quote=F, sep=" ",
                na=".", row.name=F)

  # make output object
  burnin <- tbl[tbl$niter < 0 & tbl$iqtl == 1,
                c("niter", "nqtl","LOD", "mu", "sigmasq", "addvar", "domvar", "esth")]
  iter <- tbl[tbl$niter >= 0 & tbl$iqtl == 1, c("niter", "nqtl",
        "LOD", "mu", "sigmasq", "addvar", "domvar", "esth")]
  names(burnin) <- names(iter) <- c("niter", "nqtl", "LOD",
        "mean", "envvar", "addvar", "domvar", "herit")
  no.dom <- all(is.na(iter$domvar))
  if (no.dom)
    iter$domvar <- NULL
  cols <- c("niter", "nqtl", "chrom", "locus", "add")
  if (!no.dom)
    cols <- c(cols, "dom")
  loci <- tbl[tbl$niter >= 0, cols]
  sim <- list(burnin = burnin, iter = iter,
              loci = loci, bmapqtl = .bmapqtl.options)
  class(sim) <- "bim"
  sim

}
