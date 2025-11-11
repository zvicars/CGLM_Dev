EXEDIR=$1
OFILE=$2
BOX_X=$3
BOX_Y=$4
BOX_Z=$5
$EXEDIR/genphi -i phi.input -o $OFILE.phi -box $BOX_X $BOX_Y $BOX_Z -int 7 -gs 0.184
$EXEDIR/phi2xyz -f $OFILE.phi -o $OFILE.xyz