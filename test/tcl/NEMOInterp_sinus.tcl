#!/bin/sh
# the next line restarts using tclsh \
exec $SPI_PATH/tclsh "$0" "$@"
#============================================================================
# Environnement Canada
# Centre Meteorologique Canadien
# 2100 Trans-Canadienne
# Dorval, Quebec
#
# Projet     : Exemple de scripts.
# Fichier    : NEMOInterp.tcl
# Creation   : Mars 2016 - J.P. Gauthier - CMC/CMOE
# Description: Test NEMO ORCA grid interpolation
# Parametres :
#
# Retour:
#
# Remarques  :
#    . ssmuse-sh -d rpn/OCEAN/cstint-3.2.7
#    . ssmuse-sh -x eccc/cmd/cmds/apps/SPI/beta
#============================================================================

package require TclGeoEER
package require Logger

Log::Start [info script] 0.1

set orca   [lindex $argv 0]
set latlon [lindex $argv 1]

file delete -force ./out.csintrp.avg
puts [time {
exec cstintrp -ns XX -fs ${orca} -fstype custom -fr ${latlon} -nr " " -fdtype rpn -fd ./out/out.csintrp.avg -uwgdir /tmp/.cstintrp.avg 2>@1   
}]

file delete -force ./out.csintrp

puts [time {  
exec cstintrp -ns XX -fs ${orca} -fstype custom -intyp bilin -fr ${latlon} -nr " " -fdtype rpn -fd ./out/out.csintrp -uwgdir /tmp/.cstintrp 2>@1 
}]

file delete -force ./out.spi 

puts [time {
fstdfile open NEMO read ${orca}
fstdfile open GRID read ${latlon}
fstdfile open OUT write ./out/out.spi

fstdfield read LLGRID GRID -1 "" -1 -1 -1 "" " "
fstdfield stats LLGRID -nodata 0.0

foreach var { XX } {

   Log::Print INFO "Processing $var"
   
   foreach idx [fstdfield find NEMO -1 "" -1 -1 -1 "" "$var"] {
      fstdfield read FROM NEMO $idx     
      fstdfield gridinterp LLGRID FROM LINEAR index
      fstdfield write LLGRID OUT -32 True
   }
}

fstdfile close NEMO OUT GRID
}]
Log::End
