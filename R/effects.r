#####################################################################
##
## $Id: effects.R,v 1.0 2002/07/12 yandell@stat.wisc.edu Exp $
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
bim.qtl <- function( x, cross = bim.cross( x ),
                    nqtl=1, pattern=NULL, exact = FALSE, chr,
                    bw = 2, levels = seq( 0.5, 0.95, by = 0.05 ))
{
  require(qtl)
  require(modreg)

  ## subset x and cross
  if (!is.null(pattern)) 
    nqtl <- max(nqtl, length(pattern))
  x <- subset(x, cross, nqtl, pattern, exact, chr )
  if( is.numeric( pattern ))
    pattern <- names( cross$geno )[pattern]
  cross <- subset( cross, chr )
  map <- pull.map( cross )
  pattern <- match( pattern, names( map ))
  pattern <- pattern[ !is.na( pattern ) ]

  ## smooth estimate of loci density by chromosome
  dens <- tapply( x$loci$locus, x$loci$chrom,
                 function( data, bw, len ) {
                   tmp <- density( data, bw = bw )
                   ## need to rescale height so total area is 1
                   tmp$y <- tmp$y * length( data ) / len
                   tmp
                 }, bw, length( x$loci$locus ),
                 simplify = FALSE )
  
  ## HPD region across genome: find density critical values
  y <- unlist( lapply( dens, function( chr ) chr$y ))
  o <- order( -y )
  p <- cumsum( unlist( lapply( dens, function( chr ) chr$y / ( chr$x[2] - chr$x[1] )))[o] )
  p <- p / max( p )

  hpd <- numeric( length( levels ))
  names( hpd ) <- levels
  for (i in seq(along = levels)) {
    tmp <- p <= levels[i]
    if( !any( tmp ))
      tmp <- 1
    hpd[i] <- min(y[o][tmp])
  }

  smo <- lapply( dens, function( tmp, nqtl ) {
    ## find where density peaks by downturn on both sides
    n <- length(tmp$y)
    index <- (tmp$y >= c(0, tmp$y[-n])) & (tmp$y >= c(tmp$y[-1], 0))
    index <- seq(n)[index]
    o <- order(-tmp$y[index])
    x <- tmp$x[index][o][seq(nqtl)]
    y <- tmp$y[index][o][seq(nqtl)]
    data.frame( x = x, y = y ) }, nqtl )

  findpeaks <- function( smo, nqtl )
  {
    peaks <- data.frame(chr = rep( names( smo ), rep(nqtl,length(smo))),
                        x = unlist( lapply( smo, function( chr ) chr$x )),
                        y = unlist( lapply( smo, function( chr ) chr$y )))
    peaks <- peaks[ !is.na( peaks$x ), ]
    peaks <- peaks[ order( -peaks$y )[ seq( nqtl ) ], ]
    peaks
  }
  if( !is.null( pattern )) {
    tbl <- table( pattern )
    ## estimate density across genome and find putative QTL peaks
    toploci <- data.frame( chr = ordered( character( nqtl ), names( smo )),
                          x = numeric( nqtl ), y = numeric( nqtl ))
    nloci <- 0
    for( i in names( tbl )) {
      toploci[ nloci + seq( tbl[i] ), c("x","y") ] <- smo[[i]][ seq( tbl[i] ), ]
      toploci[ nloci + seq( tbl[i] ), "chr" ] <- i
      smo[[i]][ seq( tbl[i] ), ] <- rep( NA, 2 )
      nloci <- nloci + tbl[i]
    }
    ## if nqtl larger than pattern, find next highest peaks
    if( nloci < nqtl )
      toploci[ nloci + seq( nqtl - nloci ), ] <-
        findpeaks( smo, nqtl )[ seq( nqtl - nloci ), ]
  }
  else
    toploci <- findpeaks( smo, nqtl )
  names( dens ) <- names( map )[ as.numeric( names( dens )) ]
  toploci$chr <- ordered( names( map )[ as.numeric( as.character( toploci$chr )) ],
                         names( map ))
  names( toploci ) <- c("chr","loci","dens")
  qtl <- list( loci = toploci, dens = dens, hpd = hpd )
  class( qtl ) <- "bim.qtl"
  qtl
}
##############################################################################
bim.effects <- function (x, cross = bim.cross( x ),
                         nqtl = 1, pattern = NULL, exact = FALSE, chr,
                         bw = 2, 
                         qtl = bim.qtl( x, cross, nqtl, pattern, exact,, bw ))
{
  require(qtl)
  require(modreg)
  estfn <- function(locus, add, chrom, loci, nchr ) {
    nloci <- names(loci)
    est <- sd <- double(length(nloci))
    names(est) <- nloci
    smo <- list()
    for (i in sort( unique(chrom))) {
      ii <- i == chrom
      esti <- nchr[i] == nloci
      smo[[ nchr[i] ]] <- bim.smooth(locus[ii], add[ii])
      ## now find estimate where this overlaps with loci!
      if (any(esti)) 
        for (j in seq(esti)[esti]) {
          x <- abs( smo[[ nchr[i] ]]$x - loci[j])
          est[j] <- mean( smo[[ nchr[i] ]]$y[x == min(x)])
          sd[j] <- sqrt( mean(smo[[ nchr[i] ]]$sd[x == min(x)]^2))
        }
    }
    list( smo = smo, est = est, sd = sd )
  }
  
  ## subset bim and cross
  if (!is.null(pattern)) 
    nqtl <- max(nqtl, length(pattern))
  x <- subset(x, cross, nqtl, pattern, exact, chr )
  if( !is.null( pattern ) & is.numeric( pattern ))
    pattern <- names( cross$geno )[pattern]
  cross <- subset( cross, chr )
  nchr <- names( cross$geno )
  if( !is.null( pattern )) {
    pattern <- match( pattern, nchr )
    pattern <- pattern[ !is.na( pattern ) ]
  }
  domhere <- !is.na(match("dom", names(x$loci)))

  ## grand mean
  mean <- mean( x$iter$mean )
  mean.sd <- sqrt( var( x$iter$mean ))
  loci <- qtl$loci$loci
  names( loci ) <- as.character( qtl$loci$chr )
  
  ## additive effect
  est <- estfn(x$loci$locus, x$loci$add, x$loci$chrom, loci, nchr )
  qtl$add <- est$smo

  tmp <- data.frame(chrom = c(names(loci),"mean"),
                    loci = c(qtl$loci$loci,NA),
                    add = c(est$est,mean),
                    add.sd = c(est$sd,mean.sd) )
  
  ## dominance effect if present
  if (domhere) {
    est <- estfn(x$loci$locus, x$loci$dom, x$loci$chrom, loci, nchr )
    tmp$dom <- c(est$est,NA)
    tmp$dom.sd <- c(est$sd,NA)
    qtl$dom <- est$smo
  }
  qtl$est <- tmp
  ## make sure object is of class bim.qtl
  class( qtl ) <- "bim.qtl"
  invisible( qtl )
}
##############################################################################
summary.bim.qtl <- function( object, ... )
{
  if( !is.null( object$loci )) {
    cat( "\nQTL loci and density peaks:\n" )
    print( object$loci )
  }
  if( !is.null( object$hpd )) {
    cat( "\nHPD region density cutoffs:\n" )
    print( object$hpd )
  }
  if( !is.null( object$est )) {
    cat( "\nQTL loci and effect estimates:\n" )
    print( object$est )
  }
  if( !is.null( object$dens )) {
    cat( "\nQTL density estimates by chromosome at",
        length( object$dens[[1]]$x ), "grid points with bw =",
        object$dens[[1]]$bw, "\n" )
  }
  if( !is.null( object$add )) {
    cat( "\nSmoothing spline parameters for additive effects:\n" )
    print( unlist(lapply( object$add, function(x) x$spar )))
  }
  if( !is.null( object$dom )) {
    cat( "\nSmoothing spline parameters for dominance effects:\n" )
    print( unlist(lapply( object$add, function(x) x$spar )))
  }
  invisible()
}
##############################################################################
plot.bim.effects <- function (x, cross = bim.cross( x ),
                              nqtl = 1, pattern = NULL, exact = FALSE, chr,
                              bw = 2,
                              qtl = bim.effects(x, cross,
                                nqtl, pattern, exact,, bw ),
                              cex = bim.cex(x), level = .80,
                              project = substitute(x), main = mains,
                              mfcol = c(2 + domhere, 1), ...) 
{
  project <- project
  require(qtl)
  require(modreg)
  mpos <- function(cross, cumchrlen, loci, est = rep(0, length(loci))) {
    ## place marker positions on map
    map <- pull.map( cross )
    usr <- par("usr")[3]
    for (i in seq(length(cumchrlen) - 1))
      points(cumchrlen[i] + map[[i]], rep(usr, length( map[[i]])), 
             pch = 2, col = "purple", lwd = 3)

    if (length(cumchrlen) > 2) 
      abline(v = cumchrlen - 2.5)
    points(loci, est, col = "red", lwd = 3, cex = 2)
    abline(v = loci, lty = 2, col = "red", lwd = 3)
    cchrlen <- cumchrlen
    cchrlen <- (cchrlen[-length(cchrlen)] + cchrlen[-1] - 
                5)/2
    mtext( names( cchrlen ), 1, 0.25, at = cchrlen,
          cex = 1 / ceiling( length( cchrlen ) / 15 ))
  }
  plotfn <- function(locus, add, chrom, main, ylab, smo, cumchrlen, cex, ...) {
    plot(locus, add, cex = cex, bty = "l", col = "grey40", 
         xlim = par( "usr" )[1:2], xaxs = "i", xlab = "", ylab = "", ...)
    mtext(ylab, 2, 2, cex = 1)
    mtext(main, 3, 1)

    for (i in names( smo )) {
      lines(cumchrlen[i] + smo[[i]]$x, smo[[i]]$y, lwd = 3, col = "blue")
      lines(cumchrlen[i] + smo[[i]]$x, smo[[i]]$y + 2 * smo[[i]]$sd, lwd = 3, col = "blue", lty = 2 )
      lines(cumchrlen[i] + smo[[i]]$x, smo[[i]]$y - 2 * smo[[i]]$sd, lwd = 3, col = "blue", lty = 2 )
    }
    abline(h = 0)
  }
  ## subset bim and cross
  x <- subset(x, cross, nqtl, pattern, exact, chr )
  if( !is.null( pattern ) & is.numeric( pattern ))
    pattern <- names( cross$geno )[pattern]
  cross <- subset( cross, chr )
  if( !is.null( pattern )) {
    fullpattern <- pattern
    pattern <- match( pattern, names( cross$geno ))
    pattern <- pattern[ !is.na( pattern ) ]
    nqtl <- length(pattern)
  }
  else
    fullpattern <- NULL
  domhere <- !is.na(match("dom", names(x$loci)))

  ## string chromosomes together with 5cM gap between
  map <- pull.map( cross )
  chrlen <- unlist(lapply(map, max))
  cumchrlen <- c(0, cumsum(5 + chrlen))
  names(cumchrlen) <- c(names(chrlen), "xxx")
  locus <- x$loci$locus + cumchrlen[x$loci$chrom]

  ## QTL loci
  loci <- qtl$loci$loci
  names( loci ) <- as.character( qtl$loci$chr )
  loci <- loci + cumchrlen[names(loci)]

  ## titles and plot setup
  if( is.null( pattern )) {
    mains <- paste( project, "summaries with number or QTL" )
    tmp <- nqtl
  }
  else {
    mains <- paste( project, "summaries with pattern" )
    tmp <- paste( as.character( fullpattern ), collapse = "," )
  }
  mains <- if( exact )
    as.expression( substitute( paste( main, dum == nqtl ),
        list( main = mains, dum = "", nqtl = tmp )))
    else
      as.expression( substitute( paste( main, dum >= nqtl ),
        list( main = mains, dum = "", nqtl = tmp )))
  tmp <- length(main)
  if (tmp < 3) 
    main[(1 + tmp):3] <- ""
  tmpar <- if (!is.null(mfcol)) 
    par(mfcol = mfcol, mar = c(3.1, 3.1, 3.1, 0.1))
  else par(mar = c(3.1, 3.1, 3.1, 0.1))
  on.exit(par(tmpar))

  ## conditional loci histogram and effect scatter plot
  aa <- hist(locus, breaks = seq(0, ceiling(max(cumchrlen)), 1),
             prob = TRUE, xlab = "", main = "", ylab = "", ...)
  mtext("loci histogram", 2, 2, cex = 1)
  mtext(main[1], 3, 1)
  ## HPD region (set level to 0.80 if no match)
  level <- match( round( level, 2 ), names( qtl$hpd ), nomatch = 7 )
  hpd <- qtl$hpd[level]
  for( i in names( chrlen )) {
    if( length( qtl$dens[[i]]$x )) {
    ## density lines
      lines( cumchrlen[i] + qtl$dens[[i]]$x, qtl$dens[[i]]$y, col = "blue",
            lwd = 3 )
      ## HPD regions
      in.hpd <- qtl$dens[[i]]$y > hpd
      points( cumchrlen[i] + qtl$dens[[i]]$x[in.hpd], rep( 0, sum( in.hpd )),
             col = "red" )
    }
  }
  mpos(cross, cumchrlen, loci)

  ## additive effects
  plotfn(locus, x$loci$add, x$loci$chrom, 
                main[2], "additive", qtl$add, cumchrlen, cex, ...)
  mpos(cross, cumchrlen, loci, qtl$est$add[ seq( nqtl ) ] )

  ## dominance effects if present
  if (domhere) {
    plotfn(locus, x$loci$dom, x$loci$chrom, 
                  main[3], "dominance", qtl$dom, cumchrlen, cex, ...)
    mpos(cross, cumchrlen, loci, qtl$est$dom[ seq( nqtl ) ] )
  }
  invisible( qtl )
}
##############################################################################
plot.bim.qtl <- function (x, cross = bim.cross( x ),
                          nqtl = 1, pattern = NULL, exact = FALSE, chr, bw = 2,
                          qtl = bim.qtl( x, cross, nqtl, pattern, exact,, bw ),
                          level = .8, col = "black", add = FALSE, ...) 
{
  require(qtl)
  require(modreg)
  mpos <- function(cross, cumchrlen ) {
    ## place marker positions on map
    map <- pull.map( cross )
    usr <- par("usr")[3]
    for (i in seq(length(cumchrlen) - 1))
      points(cumchrlen[i] + map[[i]], rep(usr, length( map[[i]])), 
             pch = 2, col = "purple", lwd = 3)

    if (length(cumchrlen) > 2) 
      abline(v = cumchrlen - 2.5)
    cchrlen <- cumchrlen
    cchrlen <- (cchrlen[-length(cchrlen)] + cchrlen[-1] - 
                5)/2
    mtext( names( cchrlen ), 1, 0.25, at = cchrlen,
          cex = 1 / ceiling( length( cchrlen ) / 15 ))
  }
  ## subset bim and cross
  x <- subset(x, cross, nqtl, pattern, exact, chr )
  cross <- subset( cross, chr )

  ## string chromosomes together with 5cM gap between
  map <- pull.map( cross )
  chrlen <- unlist(lapply(map, max))
  cumchrlen <- c(0, cumsum(5 + chrlen))
  names(cumchrlen) <- c(names(chrlen), "xxx")

  ## new plot?
  if( !add ) {
    tmpar <- par(mar = c(3.1, 3.1, 3.1, 0.1))
    rangey <- range( unlist( lapply( qtl$dens, function( chr ) chr$y )))
    plot( c(0,max( cumchrlen )), rangey,
         type = "n", xlab = "", ylab = "" )
    mtext("loci histogram", 2, 2, cex = 1)
    mtext("position in cM along genome", 1,2)
    mpos(cross, cumchrlen )
  }
  ## density lines
  for (i in names(chrlen))
    lines(qtl$dens[[i]]$x + cumchrlen[i], qtl$dens[[i]]$y, col = col, lwd = 3)
  ## HPD region (set level to 0.80 if no match)
  level <- match( round( level, 2 ), names( qtl$hpd ), nomatch = 7 )
  hpd <- qtl$hpd[level]
  abline( h = hpd, col = col, lty = 2 )

  invisible( qtl )
}
