#!/bin/sh

# mc_asia EMBAR

# helper script, to add jobs to rocket cluster
# takes one parameter -- $1 is the ex number we want to run.

# # if filename is empty, exit
# if [[ -z "$1" ]]
# then
#   echo "No filename given, exiting script"
#   exit
# fi

results=results/results_mc_asia.csv

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
  # echo "#pragma once
# double comparison = $(./bin/binom_vanilla_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  for N in  1000 5000 10000 50000 100000 #1000000 # paths
  do
    for M in 100 1000 10000 # steps in paths
    do

        echo "#!/bin/bash
      #SBATCH -J MC_ASIA_SERIAL_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}
      #SBATCH -N 1
      #SBATCH --ntasks-per-node=1
      #SBATCH --mem=2048
      #SBATCH -t 00:10:00

      # The directory, where programs exist
      cd $HOME/parcompfin
      make mc_asia_bin
      ./bin/mc_asia ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

      echo "Serial N=${N} M=${M} -- DONE"
      for thread in 1 2 4 8 16 32 64 #128
      do

        echo "#!/bin/bash
      #SBATCH -J MC_ASIA_OMP_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${thread}
      #SBATCH -N 1
      #SBATCH --ntasks-per-node=${thread}
      #SBATCH --mem=2048
      #SBATCH -t 00:10:00

      # The directory, where programs exist
      cd $HOME/parcompfin
      make mc_asia
      ./bin/mc_asia_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

        echo "OMP N=${N} M=${M}, thread=${thread} -- DONE"
      done
      for p in 1 2 4 8 16 32 64 #128
      do

        echo "#!/bin/bash
      #SBATCH -J MC_ASIA_MPI_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}
      #SBATCH -N ${p}
      #SBATCH --ntasks-per-node=1
      #SBATCH --mem=2048
      #SBATCH -t 00:10:00

      # The directory, where programs exist
      cd $HOME/parcompfin
      make mc_asia
      mpirun -n $p --hostfile hostfile ./bin/mc_asia_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait
        
        echo "MPI N=${N} M=${M}, processes=${p} -- DONE"
        if [ ${p} == 1 ];
        then
            echo "#!/bin/bash
          #SBATCH -J mc_asia_HYBRID_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}_${thread}
          #SBATCH -N 1
          #SBATCH --ntasks-per-node=1
          #SBATCH --mem=2048
          #SBATCH -t 00:10:00

          # The directory, where programs exist
          cd $HOME/parcompfin
          make mc_asia
          mpirun -n $p --hostfile hostfile ./bin/mc_asia_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} $thread >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait
        elif [ ${p} -le 8 ];
        then
          for thread in 2 4 8 #16 32 64 128
          do
              echo "#!/bin/bash
            #SBATCH -J mc_asia_HYBRID_${payoff_fun}_${S0}_${E}_${r}_${sigma}_${T}_${N}_${p}_${thread}
            #SBATCH -N ${p}
            #SBATCH --ntasks-per-node=${thread}
            #SBATCH --mem=2048
            #SBATCH -t 00:10:00

            # The directory, where programs exist
            cd $HOME/parcompfin
            make mc_asia
            mpirun -n $p --hostfile hostfile ./bin/mc_asia_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} $thread >> ${results}" > parcomfin.sbatch
        # sbatch parcomfin.sbatch && wait

            echo "Hybrid N=${N} M=${M}, processes=${p}, thread=${thread} -- DONE"
          done
        fi
      done
      echo "mc amer ${payoff_fun} w/ E=${E}: DONE -- $N"
    done
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

