#####################################################################
##
## $Id: zzz.R,v 1.1 2003/09/12 yandell@stat.wisc.edu Exp $
##
## Part of the R/bim package
##
## .First.lib is run when the package is loaded with library(bim)
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

.First.lib <- function(lib, pkg) {
  library.dynam("bim", pkg, lib)
  require(qtl)
  require(modreg)
  require(mva)
}

# end of zzz.R

