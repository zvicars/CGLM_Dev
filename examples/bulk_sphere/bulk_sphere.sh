INSTALL_REL=../../
INSTALL_FOLDER=$(cd $INSTALL_REL; pwd)
EXE_FOLDER=$INSTALL_FOLDER/bin
echo "Starting Simulation"
$EXE_FOLDER/CGLM run_bulk_sphere.txt > log.txt 2> err.txt
echo "Done"
echo "Creating XYZ file"
$EXE_FOLDER/b2xyz -b 0 -e 5000 -s 1 -f traj_out.traj -o traj_out.xyz 
echo "Done"
echo "Running Analysis"
$EXE_FOLDER/analysis analysis.input > analysis_log.txt 2> analysis_err.txt
echo "Done"