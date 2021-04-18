#!/bin/sh

results=results/results_mc_eur.csv

S0=100
r=0.02
sigma=0.75
T=1

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,T_overall,T_calculation,Result,Error" > ${results}

common_cycle(){
  payoff_fun=$1
  E=$2
  echo "#pragma once
double comparison = $(./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  make mc_eur_bin
  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
    ./bin/mc_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results} 
    echo "Serial N=${N} -- DONE"
  done

  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
    for thread in 1 2 4 8 16 32 64 #128
    do
      ./bin/mc_eur_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${thread} >> ${results}
        echo "OMP N=${N}, thread=${thread} -- DONE"
    done
  done

  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
    for p in 1 2 4 8 16 32 64 #128
    do
      mpirun -np $p --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}
      echo "MPI N=${N}, processes=${p} -- DONE"
    done
  done

  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
    mpirun -np 1 --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} 1 >> ${results}
    echo "Hybrid N=${N},processes=1, thread=1 -- DONE"
  done

  for N in 1000000 2000000 5000000 10000000 25000000 50000000 100000000
  do
    for p in 2 4 8 #16 32 64 #128
    do
      for thread in 2 4 8 #16 32 64 #128
      do
        mpirun -np $p --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${thread} >> ${results}
        echo "Hybrid N=${N}, processes=${p}, thread=${thread} -- DONE"
      done
    done
  done
  echo "MC eur ${payoff_fun} w/ ${E}: DONE -- $N"
}

for E in 90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

for E in 110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E}
done

