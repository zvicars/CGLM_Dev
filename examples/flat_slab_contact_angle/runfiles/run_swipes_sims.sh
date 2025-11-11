XSTAR=125000
KAPPA=0.00001
for i in 0.2 #0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 (truncated for simplicity)
do
./launch_swipes_sim.sh $XSTAR $KAPPA $i #& (backgrounding for running multiple sims simultaneously)
done