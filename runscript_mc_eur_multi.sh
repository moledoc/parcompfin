#!/bin/sh

results=results/results_mc_eur_multi.csv

# S0=100
# r=0.02
# sigma=0.75
# T=1

S0=100
r=0.1
sigma=0.2
T=1
rho=0.5
assets=(4 10)
compares=(11.90 11.60) # from article

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,T_overall,T_calculation,Result,Error" > ${results}

common_cycle(){
  payoff_fun=$1
  E=$2
  asset=$3
  compare=$4
  Ns=(1000000 2000000 5000000 10000000 25000000 50000000 75000000 100000000)
  thr=(1 5 10 25 32 50 64 100 125)
  proc=(1 5 10 25 32 50 64 100 125)
  hybr=(1 10 25 50)
  echo "#pragma once
double comparison = ${compare};" > include/comparison.h
  make mc_eur_multi_bin
  for N in ${Ns[@]} 
  do
    ./bin/mc_eur_multi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} >> ${results} 
    echo "Serial N=${N} -- DONE"
  done

  for N in ${Ns[@]} 
  do
    for thread in ${thr[@]} 
    do
      ./bin/mc_eur_multi_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} ${thread} >> ${results}
        echo "OMP N=${N}, thread=${thread} -- DONE"
    done
  done

  for N in ${Ns[@]}
  do
    for p in ${proc[@]}
    do
      mpirun -np ${p} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_multi_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} >> ${results}
      echo "MPI N=${N}, processes=${p} -- DONE"
    done
  done

  for N in ${Ns[@]}
  do
    for hybrid in ${hybr[@]}
    do
      mpirun -np ${hybrid} --hostfile hostfile --mca btl_base_warn_component_unused 0 ./bin/mc_eur_multi_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} ${hybrid} >> ${results}
      echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
    done
  done
  echo "MC eur ${payoff_fun} w/ ${E}: DONE -- $N"
}

for E in 100 #110 #90 95 100 105 110 120 130 140
do
  # for asset in ${assets[@]}
  for i in 0 1
  do
    common_cycle "call" ${E} ${assets[${i}]} ${compares[${i}]}
  done
done

# for E in 90 #110 105 100 95 90 80 70 60
# do
#   common_cycle "put" ${E}
# done

