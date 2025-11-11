#generates the best SWIPES system with given input parameters
UTILITIES_DIR=$1
X=$2
Y=$3
Z=$4
EPS=$5
NX=$(python3 -c "print(round($X/0.184))")
NY=$(python3 -c "print(round($Y/0.184))")
NZ=$(python3 -c "print(round($Z/0.184))")
WD=$PWD
$WD/generate_phiatoms.sh $EPS $X $Y $Z
$WD/generate_phifile.sh $UTILITIES_DIR out $NX $NY $NZ