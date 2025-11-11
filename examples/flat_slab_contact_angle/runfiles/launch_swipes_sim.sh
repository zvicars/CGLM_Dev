XSTAR=$1
KAPPA=$2
EPS=$3

X=4
Y=14
Z=30
EXEREL=../../../bin
EXEDIR=$(cd $EXEREL; pwd)
NX=$(python3 -c "print(round($X/0.184))")
NY=$(python3 -c "print(round($Y/0.184))")
NZ=$(python3 -c "print(round($Z/0.184))")

NX1=$(python3 -c "print(round($X/0.184)-1)")
NY1=$(python3 -c "print(round($Y/0.184)-1)")
NZ1=$(python3 -c "print(round($Z/0.184)-1)")
SLAB_WIDTH_Z=$(python3 -c "print(0.5*$Z - 7.0)")
SLAB_POS_Z=$(python3 -c "print($Z/2.0)")
NZ1BEGIN=$(python3 -c "print(round(($SLAB_POS_Z-$SLAB_WIDTH_Z + 1)/0.184))")
NZ1END=$(python3 -c "print($NZ1BEGIN+5)")
NZ2BEGIN=$(python3 -c "print(round(($SLAB_POS_Z+$SLAB_WIDTH_Z - 2)/0.184))")
NZ2END=$(python3 -c "print($NZ2BEGIN+5)")

mkdir temp
sed -e "s/#NX#/${NX}/g" \
-e "s/#NY#/${NY}/g" \
-e "s/#NZ#/${NZ}/g" \
-e "s/#NX1#/${NX1}/g" \
-e "s/#NY1#/${NY1}/g" \
-e "s/#NZ1#/${NZ1}/g" \
-e "s/#NZ1BEGIN#/${NZ1BEGIN}/g" \
-e "s/#NZ2BEGIN#/${NZ2BEGIN}/g" \
-e "s/#NZ1END#/${NZ1END}/g" \
-e "s/#NZ2END#/${NZ2END}/g" \
-e "s/#KAPPA#/${KAPPA}/g" \
-e "s/#XSTAR#/${XSTAR}/g"  run_template.input > temp/run.input

UTILITIES_DIR=$EXEDIR
TOP_DIR=$PWD
cd prep_swipes
./create_system.sh $UTILITIES_DIR $X $Y $Z $EPS 
mv out.phi ../temp/in.phi
mv out.xyz ../temp/in.xyz
cd $TOP_DIR

mkdir $EPS
cp ./temp/* ./${EPS}/
cp ./analysis.input ./${EPS}/analysis.input
rm -r ./temp
cd ./${EPS}

echo "Beginnining simulation with epsilon = ${EPS}"

$EXEDIR/CGLM run.input > log.txt 2> err.txt
$EXEDIR/analysis analysis.input > analysis_log.txt 2> analysis_err.txt
cp liq_average.stl $TOP_DIR/liquid.stl
cp solid_surfaces_average.stl $TOP_DIR/solid.stl
cd $TOP_DIR
echo "Finished"