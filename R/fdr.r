#####################################################################
##
## $Id: fdr.R,v 1.0 2003/09/14 yandell@stat.wisc.edu Exp $
##
##     Copyright (C) 2003 Brian S. Yandell
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
bim.fdr <- function( x, cross, nqtl = 1, pattern=NULL, exact=FALSE, chr, ...,
                    levels = seq( 0.01, 0.99, by = 0.01 ), df = 3,
                    qtl = bim.qtl( x, cross, nqtl, pattern, exact, ..., levels = levels ))
{
  levels <- names( qtl$hpd )

  ## subset bim and cross
  if (!is.null(pattern)) 
    nqtl <- max(nqtl, length(pattern))
  x <- subset(x, cross, nqtl, pattern, exact, chr )
  if( is.numeric( pattern ))
    pattern <- names( cross$geno )[pattern]
  cross <- subset( cross, chr )
  map <- pull.map( cross )
  pattern <- match( pattern, names( map ))
  pattern <- pattern[ !is.na( pattern ) ]

  ## collapse map for internal use
  map <- unlist( lapply( map, function( x ) diff( range( x ))))
  ## size under H0: no QTL at any locus
  size <- lapply( qtl$dens, function( ch, hpd, levels ) {
    width <- diff( ch$x[1:2] )
    size <- numeric( length( levels ))
    names( size ) <- levels
    for( i in levels )
      size[i] <- width * sum( ch$y > hpd[i] )
    size
  }, qtl$hpd, levels )
  size <- as.data.frame( size )
  nal <- names( size )
  size$all <- size[[1]] * 0
  for( i in nal) {
    tmp <- pmin( map[i], size[[i]] )
    size[[i]] <- pmin( 1, tmp / map[i] )
    size$all <- size$all + tmp
  }
  size$all <- size$all / sum( map )
  ## power under H1: QTL at locus = HPD levels
  levels <- as.numeric( levels )
  
  ## prior probability of no QTL at locus
  require(modreg)
  prob <- ( 1 - levels ) / ( 1 - size$all )
  prob[ prob == Inf ] <- NA
  tmp <- !is.na( prob )
  spline <- smooth.spline( size$all[tmp], prob[tmp], df = df )
  hyp <- c( H0 = min( spline$y ), M0 = mean( x$iter$nqtl == 0 ))
  hyp["M1"] <- 1 - hyp["M0"]
  print( hyp )
  
  ## positive false discovery rate (my Bayesian spin)
  pfdr <- size
  for( i in nal) {
    pfdr[[i]] <- hyp["H0"] * size[[i]] /
      ( hyp["M0"] * size[[i]] + hyp["M1"] * levels )
  }
  pfdr$all <- size$all / ( hyp["M0"] * size$all + hyp["M1"] * levels )
  invisible( list( levels = levels, size = size, fdr = pfdr,
                  hyp = hyp, prob = prob, spline = spline ))
}
################################################################################
plot.bim.fdr <- function( x, cross, ..., fdr = bim.fdr( x, cross, ... ),
                          critical.value = seq(0.05,0.25,by=.05), hpd = NULL )
{
  par( mfrow = c(1,2), mar = c(3.1,3.1,0.2,3.1))
  ## plot estimate prior of no QTL at any given locus
  ## as limit as size -> 0
  plot(fdr$size$all,fdr$prob, xlab = "", ylab = "", ylim=c(0,1))
  mtext( "relative size of HPD region", 1, 2 )
  mtext( "pr( H=0 | p>size )", 2, 2 )
  lines( fdr$spline$x, fdr$spline$y,col="blue",lwd = 3)
  ## plot pFDR and size vs. power
  plot( range( fdr$levels ), c(0,1), type = "n", xlab = "", ylab = "" )
  mtext( "pr( locus in HPD | m>0 )",1,2)
  mtext( "BH pFDR(-) and size(.)", 2, 2 )
  tmp <- pretty( c(0,fdr$hyp["H0"]) )
  axis( 4, tmp/fdr$hyp["H0"], as.character(tmp ))
  mtext( "Storey pFDR(-)", 4, 2 )
  lines( fdr$levels, fdr$size$all, lwd = 3, lty = 3 )
  lines( fdr$levels, fdr$fdr$all, lwd = 3 )
  if( is.null( hpd )) {
    ncrit <- length( critical.value )
    power <- numeric( ncrit )
    for( i in seq( ncrit )) {
      tmp <- critical.value[i] / fdr$hyp["H0"]
      power[i] <- max( fdr$levels[ fdr$fdr$all <= tmp ] )
      lines( c( rep( power[i], 2 ), 1 ), c( 0, rep( tmp, 2 )), col = "red", lwd = 3 )
    }
  }
  else {
    ncrit <- length( hpd )
    power <- critical.value <- numeric( ncrit )
    for( i in seq( ncrit )) {
      tmp <- fdr$levels <= hpd[i]
      power[i] <- max( fdr$levels[tmp] )
      critical.value[i] <- max( fdr$fdr$all[tmp] )
      lines( c( rep( power[i], 2 ), 1 ), c( 0, rep( critical.value[i], 2 )), col = "red", lwd = 3 )
    }
    critical.value <- critical.value * fdr$hyp["H0"]
  }
  names( power ) <- as.character( round( critical.value, 3 ))
  list(hyp = round( fdr$hyp, 3 ), fdr = power )
}
