#!/bin/sh

results=results/results_binom_embar.csv

S0=100
r=0.02
sigma=0.75
T=1

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,T_overall,T_calculation,Result,Error" > ${results}

# payoff_fun="call"
# E=40
# echo "#pragma once
# double comparison = $(./bin/binom_embar ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h

for payoff_fun in "call" "put"
do
  for E in 50 75 95
  do
    echo "#pragma once
double comparison = $(./bin/binom_embar ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
    make binom_bin
    for N in 100 250 500 750 1000 1500
    do
      ./bin/binom_embar ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results} 
      for thread in 1 2 4 8 #16 32 64 #128
      do
        ./bin/binom_embar_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} $thread >> ${results}
      done
      for p in 1 2 4 8 #16 32 64 #128
      do
        mpirun -n $p --hostfile hostfile ./bin/binom_embar_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}
        for thread in 1 2 4 8 #16 32 64 #128
        do
          mpirun -n $p --hostfile hostfile ./bin/binom_embar_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} $thread >> ${results}
        done
      done
      echo "Binom embar ${payoff_fun} w/ E=${E}: DONE -- $N"
    done
  done
done
