#!/bin/sh

# binom embar
./run_sbatch_binom_embar.sh
echo "Binom embar -- DONE"
# mc eur
./run_sbatch_mc_eur.sh
echo "MC eur -- DONE"
# mc amer
./run_sbatch_mc_amer.sh
echo "MC amer -- DONE"
# mc asia
./run_sbatch_mc_asia.sh
echo "MC asia -- DONE"
