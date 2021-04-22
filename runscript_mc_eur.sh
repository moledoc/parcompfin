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
  Ns=(1000000 2000000 5000000 7500000 10000000 25000000 50000000 75000000 100000000)
  thr=(1 2 4 8 16 32 64 128)
  proc=(1 2 4 8 16 32 64 128)
  hybr=(1 2 4 8 16 32 64 128)
  echo "#pragma once
double comparison = $(./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 100000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  make mc_eur_bin
  for N in ${Ns} 
  do
    ./bin/mc_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results} 
    echo "Serial N=${N} -- DONE"
  done

  for N in ${Ns} 
  do
    for thread in ${thr} 
    do
      ./bin/mc_eur_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${thread} >> ${results}
        echo "OMP N=${N}, thread=${thread} -- DONE"
    done
  done

  for N in ${Ns}
  do
    for p in ${proc}
    do
      mpirun -np ${p} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}
      echo "MPI N=${N}, processes=${p} -- DONE"
    done
  done

  for N in ${Ns}
  do
    for hybrid in ${hybr}
    do
      mpirun -np ${hybrid} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/binom_embar_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${hybrid} >> ${results}
      echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
    done
  done
  echo "MC eur ${payoff_fun} w/ ${E}: DONE -- $N"
}

for E in 110 #90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

for E in 90 #110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E}
done

