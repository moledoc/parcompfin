#!/bin/sh

results=results/results_mc_eur.csv

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
double comparison = $(./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
    make mc_eur_bin
    for N in 1 2 5 10 25 50 100 #1 10 100 #
    do
      ./bin/mc_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N}000000 >> ${results} 
      for thread in 1 2 4 8 16 32 64 #128
      do
        ./bin/mc_eur_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N}000000 $thread >> ${results}
      done
      for p in 1 2 4 8 #16 32 64 #128
      do
        mpirun -n $p --hostfile hostfile ./bin/mc_eur_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N}000000 >> ${results}
        for thread in 1 2 4 8 #16 32 64 #128
        do
          mpirun -n $p --hostfile hostfile ./bin/mc_eur_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N}000000 $thread >> ${results}
        done
      done
      echo "MC eur ${payoff_fun} w/ ${E} DONE -- $N"
    done
  done
done
