#####################################################################
##
## $Id: bim.R,v 1.0 2002/07/13 yandell@stat.wisc.edu Exp $
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
read.bim <- function( dir, bimfile, nvalfile = "nval.dat", na.strings = "." )
{
  if (!missing(dir)) {
    n <- nchar(dir)
    if (substring(dir, n, n) == "/") 
      dir <- substring(dir, 0, n - 1)
    bimfile <- file.path(dir, bimfile)
    nvalfile <- file.path(dir, nvalfile)
  }
  cat( "reading", bimfile, "created by Bmapqtl\n" )
  
  bmapqtl <- read.bmapqtl( nvalfile = nvalfile )
  cat( "MCMC prior for number of QTL assumed to be", bmapqtl$prior.nqtl,
      "with mean", bmapqtl$mean.nqtl, "\n" )

  tbl <- read.table( bimfile, header = TRUE, na.strings = na.strings )
  burnin <- tbl[ tbl$niter < 0 & tbl$iqtl == 1,
    c("niter","nqtl","LOD","mu","sigmasq","addvar","domvar","esth") ]
  iter <- tbl[ tbl$niter >= 0 & tbl$iqtl == 1,
    c("niter","nqtl","LOD","mu","sigmasq","addvar","domvar","esth") ]
  names( burnin ) <- names( iter ) <- 
    c("niter","nqtl","LOD","mean","envvar","addvar","domvar","herit")
  no.dom <- all( is.na( iter$domvar ))
  if( no.dom )
    iter$domvar <- NULL
  cols <- c("niter","nqtl","chrom","locus","add")
  if( !no.dom)
    cols <- c(cols,"dom")
  loci <- tbl[ tbl$niter >= 0, cols ]
  sim <- list( bmapqtl = bmapqtl, burnin = burnin, iter = iter, loci = loci )
  class( sim ) <- "bim"
  sim
}
##############################################################################
summary.bim <- function( object, ... )
{
  cat ("Bayesian interval mapping object", substitute( object ), "\n" )
  cat( "had", as.integer( object$bmapqtl$inter ), "iterations recorded at each",
      as.integer( object$bmapqtl$by ), "steps\n" )
  burnin <- round( object$bmapqtl$burnin * object$bmapqtl$niter )
  cat( "with", as.integer( round( object$bmapqtl$preburn * burnin )),
      "pre-burn-in and", as.integer( burnin ), "burn-in steps.\n" )
  cat( paste( "Prior for number of QTL was ", object$bmapqtl$prior.nqtl,
             "(", object$bmapqtl$mean.nqtl, ").\n\n", sep = "" ))
  cat( "Percentages for number of QTL detected:" )
  print( round( 100 * table( object$iter$nqtl ) / nrow( object$iter )))
  cat( "\nDiagnostic summaries:\n" )
  print( summary( object$iter[,1] ))
  invisible( )
}
##############################################################################
bim.legacy <- function( bim )
{
  if( !is.null( bim$bmapqtl ))
    return( bim )
  bim$bmapqtl <- read.bmapqtl()
  bim$bmapqtl$prior.nqtl <- switch( bim$prior,
                                   exp =, geo =, geometric = "geometric",
                                   pois =, poisson = "poisson",
                                   unif =, uniform = "uniform" )
  bim$bmapqtl$mean.nqtl <- bim$mean.prior
  bim$prior <- NULL
  bim$mean.prior <- NULL
  bim
}
##############################################################################
bim.prior <- function( bim, range )
{
  ## legacy code
  if( is.null( bim$bmapqtl )) {
    warning( "legacy data: please run > bim <- bim.legacy( bim )" )
    bim <- bim.legacy( bim )
  }
  prior <- bim$bmapqtl$prior.nqtl
  mean <- bim$bmapqtl$mean.nqtl
  pr <- switch( prior,
         geometric = ( 1 / mean ) ^ range / ( 1 - ( 1 / mean )),
         poisson   = dpois( range, mean ),
         uniform   = rep( 1 / ( 1 + mean ), length( range )))
  if( is.null( pr ))
    stop( "only geometric, poisson, and uniform priors recognized\n" )
  names( pr ) <- range
  pr
}
##############################################################################
subset.bim <- function( x, cross = bim.cross( x ),
                       nqtl = 1, pattern = NULL, exact = FALSE, chr, ... )
{
  nqt <- nqtl[1]
  if( !is.null( pattern )) {
    nqtl <- max( nqtl, length( pattern ))
    if( is.character( pattern ))
      stop( "pattern must be numeric, not character" )
  }
  if( exact ) {
    x$iter <- x$iter[ x$iter$nqtl == nqtl, ]
    if( !nrow( x$iter ))
      stop( paste( "empty object: no iterations with number of QTL =", nqtl ))
    x$loci <- x$loci[ x$loci$nqtl == nqtl, ]
  }
  else {
    x$iter <- x$iter[ x$iter$nqtl >= nqtl, ]
    if( !nrow( x$iter ))
      stop( paste( "empty object: no iterations with number of QTL >=", nqtl ))
    x$loci <- x$loci[ x$loci$nqtl >= nqtl, ]
  }
  if( !is.null( pattern )) {
    mypat <- table( pattern )
    if( exact ) {
      mypat <- c( mypat, extra = 0 )
      patfn <- function( x, mypat, yourpat ) {
        tbl <- table( x )
        tmp <- match( names( tbl ), names( mypat ), nomatch = length( mypat ))
        yourpat[tmp] <- tbl[ tmp > 0 ]
        all( yourpat == mypat )
      }
    }
    else {
      patfn <- function( x, mypat, yourpat ) {
        tbl <- table( x )
        tmp <-  match( names( tbl ), names( mypat ), nomatch = 0 )
        yourpat[tmp] <- tbl[ tmp > 0 ]
        all( yourpat >= mypat )
      }
    }
    blank <- rep( 0, length( mypat ))
    names( blank ) <- names( mypat )
    iters <- unlist( tapply( x$loci$chrom, x$loci$niter, patfn,
                            mypat, blank, simplify = FALSE ))
    x$iter <- x$iter[iters,]
    if( !nrow( x$iter ))
      stop( paste( "empty object: no patterns like", pattern ))
    iters <- x$iter$niter
    x$loci <- x$loci[ !is.na( match( x$loci$niter, iters )), ]
  }
  if (!missing(chr)) {
    n.chr <- nchr( cross )
    if (is.logical(chr)) {
      if (length(chr) != n.chr) 
        stop(paste("If logical, chr argument must have length", 
                   n.chr))
      chr <- (1:n.chr)[chr]
    }
    if (is.numeric(chr)) {
      if (all(chr < 1)) 
        chr <- (1:n.chr)[chr]
      if (any(chr < 1 | chr > n.chr)) 
        stop("Chromosome numbers out of range.")
    }
    else {
      if (any(is.na(match(chr, names(x$geno))))) 
        stop("Not all chromosome names found.")
      chr <- match(chr, names(cross$geno))
    }
    chr <- sort( unique( chr ))
    kept <- match( seq( length( cross$geno )), chr, nomatch = 0)
    x$loci$chrom <- kept[x$loci$chrom]
    x$loci <- x$loci[x$loci$chrom > 0, ]
  }
  x
}
##############################################################################
bim.match <- function( bim, pattern )
{
  patfn <- function(x) {
    tmp <- paste( ifelse( x > 1,
                         paste( x, "*", sep=""),
                         "" ),
                 names( x ), collapse = ",", sep = "" )
    if( length( x ) > 1 )
      tmp <- paste( sum( x ), tmp, sep = ":" )
    tmp
  }
  mypat <- patfn( table( pattern ))
  counts <- tapply( bim$loci$chrom, bim$loci$niter, function( x ) table( x ),
    simplify = FALSE )
  patterns <- unlist( lapply( counts, patfn ))
  as.numeric( names( counts )[ !is.na( match( patterns, mypat )) ] )
}
##############################################################################
bim.cross <- function( bim )
{
  ## fake a cross based on bim object
  cross <- list( geno = tapply( bim$loci$locus, bim$loci$chrom,
                   function( x ) list( data = matrix(NA,0,0),
                                      map = c(0, max( x ))),
                   simplify = FALSE ),
                pheno = matrix( NA, 0, 0 ))
  for(i in seq(length(cross$geno)))
    class(cross$geno[[i]]) <- "A"
  class( cross ) <- c("f2","cross")
  names( cross$geno ) <- as.numeric( seq( along = cross$geno ))
  cross  
}
##############################################################################
bim.cex <- function( bim )
{
  2 ^ ( 2 - min( 4, max( 2, ( log10( nrow( bim$loci ))))))
}
##############################################################################
##############################################################################
plot.bim <- function( x, cross = bim.cross( x ),
                     nqtl = 1, pattern = NULL, exact = FALSE,
                     ... )
{
  cat( "time series of burnin and mcmc runs\n" )
  plot.bim.mcmc( x, ... )
  cat( "jittered plot of quantitative trait loci by chromosome...\n" )
  plot.bim.loci( x, cross, nqtl, pattern, exact, ... )
  cat( "model selection plots: number of QTL and chromosome pattern...\n" )
  model <- bim.model( x, cross, nqtl, pattern, exact, ... )
  summary( model )
  plot( model )
  cat( "quantitative trait loci (histogram) and effects (scatter plot)...\n" )
  qtl <- plot.bim.effects( x, cross, nqtl, pattern, exact, ... )
  summary( qtl )
  cat( "summary diagnostics as histograms and boxplots by number of QTL\n" )
  plot.bim.diag( x, nqtl, ... )
  invisible()
}
##############################################################################
plot.bim.mcmc <- function( x, element = c("burnin","iter"),
                          xlab = c("burnin sequence","mcmc sequence"),
                          items = names( x$iter )[-1],
                          ylabs = items,
                          types = c("b",rep("l", length( items ) - 1 )),
                          ... )
{
  require( modreg )
  tmpar <- par( mfcol = c( length( items ), length( element )),
               mar = c(3.1,3.1,0.1,0.1) )
  on.exit( par( tmpar ))

  for( s in seq( element )) for( i in seq( items )) {
    plot( x[[ element[s] ]]$niter, x[[ element[s] ]][[ items[i] ]],
         type = types[i], col = "grey20", xlab = "", ylab = "", ... )
    tmp <- bim.smooth( x[[ element[s] ]]$niter,
                      x[[ element[s] ]][[ items[i] ]] )
    lines( tmp$x, tmp$y, lwd = 3, col = "blue" )
    mtext( xlab[s], 1, 2, cex = 1 )
    mtext( ylabs[i], 2, 2, cex = 1 )
  }
  invisible()
}
##############################################################################
bim.smooth <- function( x, y )
{
  require( modreg )
  
  ux <- unique( x )
  if( length( ux ) < 50 ) {
#    smo <- list( x = sort( ux ),
#                y = rep( mean( y ), length( ux )))
#    smo$sd <- rep( mad( y ), length( smo$x ))
    lmy <- lm( y ~ x )
    smo <- list( x = sort( ux ),
                y = predict( lmy, data.frame( x = sort( ux ))),
                sd = rep( sqrt( sum( resid( lmy ) ^ 2 ) / lmy$df.resid ), length( ux )))
  }
  else {
    smo <- smooth.spline( x, y )
    smo$sd <- sqrt( pmax( 0, smooth.spline( x, ( y - predict( smo, x )$y ) ^ 2 )$y ))
  }
  smo  
}
##############################################################################
plot.bim.loci <- function( x, cross = bim.cross( x ),
                          nqtl = 1, pattern = NULL, exact = FALSE,
                          chr, labels = TRUE, amount = .35,
                          cex = bim.cex( x ), ... )
{
  require( qtl )
  amount <- max( 0, min( 0.45, amount ))
  x <- subset( x, cross, nqtl, pattern, exact, chr )
  
  cross <- subset( cross, chr )
  map <- pull.map( cross )
  nmap <- names( map )

  loci <- x$loci[ , c("chrom","locus") ]
  if( 0 == nrow( loci )) {
    warning( "no mcmc samples on chosen chromosomes: ",
      paste( chr, collapse = "," ))
    return( )
  }
  plot( range( loci$chrom ) + c(-.5,.5), range( loci$locus ),
       type = "n", xaxt = "n", xlab = "", ylab = "", ... )
  points( jitter( loci$chrom, , amount ), loci$locus, cex = cex )
  uchrom <- unique( loci$chrom )
  axis( 1, uchrom, nmap[uchrom] )
  mtext( "chromosome", 1, 2 )
  mtext( "MCMC sampled loci", 2, 2 )
  cxy <- par( "cxy" )[2] / 4

  if( !is.null( cross )) {
    for( i in uchrom ) {
      tmp <- map[[ nmap[i] ]]
      text( rep( i-.5, length( tmp )), cxy + tmp, names( tmp ),
        adj = 0, cex = 0.5, col = "blue" )
      for( j in tmp )
        lines( i + c(-.5,.5), rep( j, 2 ), col = "blue",
              lwd = 2 )
    }
  }
}
##############################################################################
bim.nqtl <- function( bim )
{
  ## posterior number of QTL
  posterior <- table( bim$iter$nqtl )
  ntrial <- sum( posterior )
  posterior <- posterior / ntrial
  ## prior number of QTL
  prior <- bim.prior( bim, as.numeric( names( posterior )))

  ## posterior/prior ratios for Bayes factor
  bf <- posterior / prior
  if( bf[1] > 0 & prior[1] > 0 )
    bf <- bf / bf[1]
  
  ## bfse = approximate Monte Carlo standard error for bf (actually binomial error )
  ## note that this is rescaled since bf[1] forced to be 1
  list( nqtl = 
       list( posterior = posterior, prior = prior, bf = bf,
            bfse = sqrt(( 1 - posterior ) / ( posterior * ntrial )) * bf ))
}
##############################################################################
bim.pattern <- function( bim, cross = bim.cross( bim ),
                        nqtl = 1, pattern = NULL, exact = FALSE,
                        cutoff = 1 )
{
  require( qtl )
  bim <- subset( bim, cross, nqtl, pattern, exact )
  loci <- bim$loci
  counts <- tapply( loci$chrom, loci$niter, function( x ) table( x ),
    simplify = FALSE )
  pattern <- unlist( lapply( counts, function(x) {
    tmp <- paste( ifelse( x > 1,
                         paste( x, "*", sep=""),
                         "" ),
                 names( x ), collapse = ",", sep = "" )
    if( length( x ) > 1 )
      tmp <- paste( sum( x ), tmp, sep = ":" )
    tmp
  } ))
  posterior <- rev( sort( table( pattern )))
  posterior <- posterior / sum( posterior )
  tmp <- posterior >= cutoff / 100
  if( sum( tmp ))
  posterior <- posterior[tmp]
  else {
    cat( "warning: posterior cutoff of ", cutoff,
        "is too large and is being ignored\n",
        "posterior range is", range( posterior ), "\n" )
  }
  if( length( posterior ) > 15 )
    posterior <- posterior[1:15]
  ucount <- match( names( posterior ), pattern )

  ## prior for pattern
  rng <- max( bim$iter$nqtl )
  pr <- bim.prior( bim, 0:rng )
  bf <- posterior
  map <- pull.map( cross )
  chrlen <- unlist( lapply( map, max ))
  nchrom <- length( chrlen )
  chrlen <- chrlen / sum( chrlen )
  
  names( chrlen ) <- seq( nchrom )
  fact <- rep( 1, rng )
  for( i in 2:(rng+1) ) 
    fact[i] <- fact[i-1] * i
  for( i in seq( posterior )) {
    ct <- counts[[ ucount[i] ]]
    st <- sum( ct )
    bf[i] <- pr[st] * prod( chrlen[ names( ct ) ] ^ ct ) *
      fact[st] / prod( fact[ct] )
  }

  ntrial <- length( pattern )
  prior <- bf
  bf <- posterior / prior
  if( bf[1] > 0 & prior[1] > 0 )
    bf <- bf / bf[1]

  ## bfse = approximate Monte Carlo standard error for bf (actually binomial error )
  ## note that this is rescaled since bf[1] forced to be 1
  list( pattern = 
       list( posterior = posterior, prior = prior, bf = bf,
            bfse = sqrt(( 1 - posterior ) / ( posterior * ntrial )) * bf ))
}
##############################################################################
bim.model <- function( bim, cross = bim.cross( bim ),
                      nqtl = 1, pattern = NULL, exact = FALSE,
                      cutoff = 1 )
{
  bim <- subset.bim( bim, cross, nqtl, pattern, exact )
  assess <- list( )
  assess$nqtl <- bim.nqtl( bim )$nqtl
  assess$pattern <- bim.pattern( bim, cross, cutoff = cutoff )$pattern
  assess$param <- list( nqtl = nqtl, pattern = pattern, exact = exact,
                       cutoff = cutoff )
  class( assess ) <- "bim.model"
  assess
}
##############################################################################
summary.bim.model <- function( object, ... )
{
  if( !is.null( object$nqtl )) {
    cat( "posterior for number of QTL as %" )
    print( round( 100 * object$nqtl$posterior ))
    cat( "Bayes factor ratios for number of QTL")
    print( round( object$nqtl$bf, 1 ))
  }
  if( !is.null( object$pattern )) {
    cat( "model posterior above cutoff", object$param$cutoff, "as %\n" )
    print( round( 100 * object$pattern$posterior))
    cat( "Bayes factor ratios for chromosome pattern\n")
    print( round( object$pattern$bf, 1 ))
  }
  invisible()
}
##############################################################################
plot.bim.model <- function( x, cross = bim.cross( x ),
                           nqtl = 1, pattern = NULL, exact = FALSE,
                           cutoff = 1 ,
                           assess = bim.model( x, cross, nqtl, pattern, exact, cutoff ), ... )
{
  if( inherits( x, "bim.model" ))
    assess <- x
  is.nqtl <- !is.null( assess$nqtl )
  is.pattern <- !is.null( assess$pattern )
  tmpar <- par( mfrow = c( is.nqtl + is.pattern, 2 ))
  on.exit( par( tmpar ))
  if( is.nqtl )
    plot.bim.pattern( assess$nqtl, as.numeric( names( assess$nqtl$posterior )),
                     c("number of QTL","QTL posterior","QTL posterior"), NULL, ... )
  if( is.pattern )
    plot.bim.pattern( assess$pattern, ... )
  assess
}
##############################################################################
plot.bim.pattern <- function( x,
                             bars = seq( x$posterior ),
                             labels = c("model index","model posterior","pattern posterior"),
                             barlabels = names(x$posterior),
                             threshold = c( weak=3, moderate=10, strong=30 ),
                             units = 2, rescale = TRUE, ... )
{
  require( qtl )
  ## plot posterior
  bar <- barplot( x$posterior, col = "white", names = bars, ... )
  tmp <- if( rescale )
    x$prior * max( x$posterior ) / max( x$prior )
  else
    x$prior
  lines( bar, tmp, type = "b", col = "blue", lwd = 2 )
  if( !is.null( barlabels )) {
    cex <- 1
    usr <- par("usr")
    tmp <- as.numeric( 2 * ( x$posterior - usr[3] ) >= diff( usr[3:4] ))
    for( i in 0:1 ) {
      ii <- i == tmp
      if( any( ii ))
        text( bar[ii], x$posterior[ii], barlabels[ii], srt=90,
             adj= i, cex = cex )
    }
  }
  mtext( labels[1], 1, 2 )
  mtext( labels[2], 2, 2 )
  mtext( labels[3], 3, 0.5 )

  x$bf[ x$bf == 0 | x$prior == 0 ] <- NA

  plot( seq( bars ), x$bf, log = "y", xaxt = "n", xlim = c( 0.5, length( bars )),
    xlab = "", ylab = "", ... )
  mtext( labels[1], 1, 2 )
  mtext( "posterior / prior", 2, 2 )
  mtext( "Bayes factor ratios", 3, 0.5 )
  axis( 1, seq( bars ), bars )

  if( !is.null( barlabels )) {
    cxy <- par( "cxy" )[1] * cex / 2
    ## decypher number of QTL from pattern--kludge!
    nqtl <- strsplit( barlabels, "" )
    nqtl <- lapply( nqtl, function( x ) {
      colon <- seq( x )[ x == ":" | x == "*" ][1]
      x <- if( is.na( colon ))
        "1"
      else
        x[ seq( colon - 1 ) ]
      paste( x, collapse = "" )
    } )
    text( seq( barlabels ) - cxy, x$bf, nqtl, srt = 90, cex = cex )
  }

  usr <- 10^par( "usr" )[3:4]
  for( i in seq( bars ) ) {
    if( x$bfse[i] > 0 ) {
      bfbar <- x$bf[i] + c(-units,units) * x$bfse[i]
      bfbar[1] <- max( usr[1], bfbar[1] )
      bfbar[2] <- min( usr[2], bfbar[2] )
    }
    else
      bfbar <- usr
    lines( rep(i,2), bfbar )
  }
  ## put threshold yardstick on plot
  if( length( threshold )) {
    bars <- floor( mean( bars ) / 2 ) + 0.5
    maxusr <- usr[2]
    usr <- prod( usr ^ c(.95,.05) )
    lines( bars + c(-.25,.25), rep( usr, 2 ), lwd = 3, col = "blue" )
    texusr <- usr
    for( i in seq( length( threshold ))) {
      sigusr <- min( maxusr, usr * threshold[i] )
      if( texusr < maxusr )
        text( bars + 0.5, sqrt( texusr * sigusr ), names( threshold )[i],
             col = "blue", adj = 0 )
      arrows( bars, usr, bars, sigusr, 0.1, lwd = 3, col = "blue" )
      texusr <- sigusr
    }
  }
  invisible( x )
}
##############################################################################
### marginal histograms
##############################################################################
plot.bim.diag <- function( x,
                     nqtl = 1, pattern = NULL, exact = FALSE,
                     items= names( x$iter )[-(1:2)],
                     mains = items,
                     mfrow = c(nhist,2), ... )
{
  x <- subset( x,, nqtl, pattern, exact )
  nhist <- length( items )
  tmpar <- if( !is.null( mfrow ))
    par( mfrow = mfrow, mar=c(3.6,4.1,0.6,0.1) )
  else
    par( mar=c(3.6,4.1,0.6,0.1) )
  on.exit( par( tmpar ))
  mains[ match( "herit", mains, nomatch = 0 ) ] <- "heritability"
  for( i in seq( nhist )) {
    ## marginal histogram
    main <- as.expression( substitute( paste( "marginal ", main, ", ",
        italic(m) >=  nqtl),
        list( nqtl = nqtl, main = mains[i] )))
    tmp <- x$iter[ , items[i] ]
    tmp <- tmp[ !is.na( tmp ) ]
    med <- quantile( tmp, c(.25,.5,.75) )
    plot( density(tmp ), main = "", xlab = "", ylab = "", ... )
    mtext( "density", 2, 2.5 )
    mtext( main, 1, 2.5 )
    ## hand-made boxplot on its side
    b <- boxplot( tmp, plot = FALSE )
    tmp <- par("usr")[4] / 8
    up <- tmp * 1.5
    polygon( b$stats[c(2,4,4,2)], up+c(0,0,tmp,tmp))
    lines( rep(b$stats[3],2), up+c(0,tmp) )
    lines( b$stats[1:2], up+rep(tmp/2,2), lty = 2 )
    lines( b$stats[4:5], up+rep(tmp/2,2), lty = 2 )
    lines( rep(b$stats[1],2), up+tmp*c(1,3)/4 )
    lines( rep(b$stats[5],2), up+tmp*c(1,3)/4 )

    ## conditional boxplots
    tmp <- split( x$iter[[ items[i] ]], x$iter$nqtl )
    boxplot( tmp )
    mtext( paste( mains[i], "conditional on number of QTL" ), 1, 2.5 )
    mtext( mains[i], 2, 2.5 )

    cat( items[i], round( med[2], 3 ), "\n" )
    cat( "conditional", mains[i], "\n" )
    print( round( unlist( lapply( tmp, median, na.rm = TRUE )), 3 ))
  }
}
