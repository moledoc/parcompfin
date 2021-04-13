# !/bin/sh

results=results/results_mc_amer.csv

S0=100
r=0.02
sigma=0.75
T=1

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,T_overall,T_calculation,Result,Error" > ${results}

for payoff_fun in "call" "put"
do
  for E in 50 75 95
  do
    echo "#pragma once
double comparison = $(./bin/binom_vanilla_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
    make mc_amer_bin
    for N in 10000 #50000 #100000 #1000000 # paths
    do
      for M in 1000 10000 # steps in paths
      do
        ./bin/mc_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results} 
        echo "Serial N=${N} M=${M} -- DONE"
        for thread in 1 2 4 #8 #16 32 64 #128
        do
          ./bin/mc_amer_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}
          echo "OMP N=${N}, M=${M}, thread=${thread} -- DONE"
        done
        for p in 1 2 4 #8 #16 32 64 #128
        do
          mpirun -n $p --hostfile hostfile ./bin/mc_amer_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}
          echo "MPI N=${N}, M=${M}, processes=${p} -- DONE"
          for thread in 1 2 4 #8 #16 32 64 #128
          do
            mpirun -n $p --hostfile hostfile ./bin/mc_amer_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} $thread >> ${results}
            ecHo "hybrid N=${N}, M=${M}, processes=${p}, thread=${thread} -- DONE"
          done
        done
        echo "MC amer ${payoff_fun} w/ E=${E}: DONE -- M=$M"
      done
      echo "MC amer ${payoff_fun} w/ E=${E}: DONE -- N=$N"
    done
  done
done
