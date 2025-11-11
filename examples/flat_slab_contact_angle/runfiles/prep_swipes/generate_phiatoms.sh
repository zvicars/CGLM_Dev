#create a plate with a specified width and height
EPS=$1
BOX_WIDTH_X=$2
BOX_WIDTH_Y=$3
BOX_WIDTH_Z=$4
SLAB_WIDTH_Z=$(python3 -c "print(0.5*$BOX_WIDTH_Z - 7.0)")
SLAB_POS_Z=$(python3 -c "print($BOX_WIDTH_Z/2.0)")
RIGHT_BOUNDARY_POS=$(python3 -c "print($BOX_WIDTH_Z-1.5)")
LEFT_BOUNDARY_POS=0.5

#substitute relevant input args into template file
sed -e "s/#EPS#/${EPS}/g" -e "s/#BWID#/${BOX_WIDTH}/g" \
    -e "s/#WID_X#/${BOX_WIDTH_X}/g"  \
    -e "s/#WID_Y#/${BOX_WIDTH_Y}/g"  \
    -e "s/#WID_Z#/${BOX_WIDTH_Z}/g"  \
    -e "s/#SWID_Z#/${SLAB_WIDTH_Z}/g"  \
    -e "s/#SLAB_Z#/${SLAB_POS_Z}/g"  \
    -e "s/#LB_POS#/${LEFT_BOUNDARY_POS}/g" \
    -e "s/#RB_POS#/${RIGHT_BOUNDARY_POS}/g" \
     phi_template.input > phi.input

