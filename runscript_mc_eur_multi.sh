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

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,Nr_of_assets,T_overall,T_calculation,Result,Abs_Error,Error" > ${results}

common_cycle(){
  payoff_fun=$1
  E=$2
  asset=$3
  compare=$4
  Ns=(10000 80000 100000 160000 320000 640000 800000 1000000 1600000 3200000 6400000 8000000 10000000 16000000 32000000 64000000 80000000 100000000)
  thr=(1 8 16 32 64)
  proc=(1 8 16 32 64)
  hybr=(1 4 8)
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
      mpirun -np ${p} --hostfile hostfile ./bin/mc_eur_multi_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} >> ${results}
      echo "MPI N=${N}, processes=${p} -- DONE"
    done
  done

  for N in ${Ns[@]}
  do
    for hybrid in ${hybr[@]}
    do
      mpirun -np ${hybrid} --hostfile hostfile ./bin/mc_eur_multi_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${asset} ${rho} ${hybrid} >> ${results}
      echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
    done
  done
  echo "MC eur multi ${payoff_fun} w/ ${E}: DONE -- $N"
}

for E in 100 #110 #90 95 100 105 110 120 130 140
do
  # for asset in ${assets[@]}
  for i in 0 #1
  do
    common_cycle "call" ${E} ${assets[${i}]} ${compares[${i}]}
  done
done


