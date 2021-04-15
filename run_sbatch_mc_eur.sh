#!/bin/sh

# BINOM EMBAR

# helper script, to add jobs to rocket cluster
# takes one parameter -- $1 is the ex number we want to run.

# # if filename is empty, exit
# if [[ -z "$1" ]]
# then
#   echo "No filename given, exiting script"
#   exit
# fi

results=results/results_mc_eur.csv

S0=100
r=0.02
sigma=0.75
T=1

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,T_overall,T_calculation,Result,Error" > ${results}

wait(){
  while [ -n "$(squeue -h -u utt)" ]
  do
    sleep 10
  done
}

common_cycle(){
  payoff_fun=$1
  E=$2
  echo "#pragma once
double comparison = $(./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
      echo "#!/bin/bash
    #SBATCH -J MC_EUR_SERIAL_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}
    #SBATCH -N 1
    #SBATCH --ntasks-per-node=1
    #SBATCH --mem=2048
    #SBATCH -t 00:10:00

    # The directory, where programs exist
    cd $HOME/parcompfin
    make mc_eur_bin
    ./bin/mc_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

    echo "Serial N=${N} -- DONE"
    for thread in 1 2 4 8 16 32 64 #128
    do

      echo "#!/bin/bash
    #SBATCH_-J_MC_EUR_OMP_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${thread}
    #SBATCH -N 1
    #SBATCH --ntasks-per-node=${thread}
    #SBATCH --mem=2048
    #SBATCH -t 00:10:00

    # The directory, where programs exist
    cd $HOME/parcompfin
    make mc_eur
    ./bin/mc_eur_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${thread} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

      echo "OMP N=${N}, thread=${thread} -- DONE"
    done
    for p in 1 2 4 8 16 32 64 #128
    do

      echo "#!/bin/bash
    #SBATCH -J MC_EUR_MPI_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}
    #SBATCH -N ${p}
    #SBATCH --ntasks-per-node=1
    #SBATCH --mem=2048
    #SBATCH -t 00:10:00

    # The directory, where programs exist
    cd $HOME/parcompfin
    make mc_eur
    mpirun -n $p --hostfile hostfile ./bin/mc_eur_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait
      
      echo "MPI N=${N}, processes=${p} -- DONE"
      if [ ${p} == 1 ];
      then 
          echo "#!/bin/bash
        #SBATCH -J MC_EUR_HYBRID_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}_${thread}
        #SBATCH -N 1
        #SBATCH --ntasks-per-node=1
        #SBATCH --mem=2048
        #SBATCH -t 00:10:00

        # The directory, where programs exist
        cd $HOME/parcompfin
        make mc_eur
        mpirun -n $p --hostfile hostfile ./bin/mc_eur_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} $thread >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait
      elif [ ${p} -le 8 ];
      then
        for thread in 2 4 8 #16 32 64 #128
        do
            echo "#!/bin/bash
          #SBATCH -J MC_EUR_HYBRID_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}_${thread}
          #SBATCH -N ${p}
          #SBATCH --ntasks-per-node=${thread}
          #SBATCH --mem=2048
          #SBATCH -t 00:10:00

          # The directory, where programs exist
          cd $HOME/parcompfin
          make mc_eur
          mpirun -n $p --hostfile hostfile ./bin/mc_eur_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} $thread >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

          echo "Hybrid N=${N}, processes=${p}, thread=${thread} -- DONE"
        done
      fi
    done
    echo "mc eur ${payoff_fun} w/ E=${E}: DONE -- $N"
  done
}


for E in 90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

for E in 110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E}
done

