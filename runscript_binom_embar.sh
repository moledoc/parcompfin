#!/bin/sh

results=results/results_binom_embar.csv

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
  Ns=(100 250 750 1000 1500 5000 10000 50000 100000) # paths
  thr=(1 5 10 25 50 100)
  proc=(1 5 10 25 50 100)
  hybr=(1 10 25 50)
  echo "#pragma once
double comparison = $(./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 100000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  make binom_bin
  for N in ${Ns[@]}
  do
    ./bin/binom_embar ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results} 
    ./bin/binom_vanilla_eur ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results} 
    echo "Serial N=${N} -- DONE"
  done

  for N in ${Ns[@]}
  do
    for thread in ${thr[@]} 
    do
      ./bin/binom_embar_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${thread} >> ${results}
      echo "OMP N=${N}, thread=${thread} -- DONE"
    done
  done

  for N in ${Ns[@]}
  do
    for p in ${proc[@]}
    do
      mpirun -np ${p} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/binom_embar_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} >> ${results}
      echo "MPI N=${N}, processes=${p} -- DONE"
    done
  done

  for N in ${Ns[@]}
  do
    for hybrid in ${hybr[@]}
    do
      mpirun -np ${hybrid} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/binom_embar_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${hybrid} >> ${results}
      echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
    done
  done
  echo "Binom embar ${payoff_fun} w/ E=${E}: DONE -- $N"
}


for E in 110 #90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

for E in 90 #110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E}
done
