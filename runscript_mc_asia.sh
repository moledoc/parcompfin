# !/bin/sh

results=results/results_mc_asia.csv

S0=100
r=0.02
sigma=0.75
T=1

echo "#pragma once
double comparison = 0;" > include/comparison.h

echo "Method,Payoff,S0,E,r,sigma,T,N,M,Parallel,Nr_of_assets,T_overall,T_calculation,Result,Abs_Error,Error" > ${results}

common_cycle(){
  payoff_fun=$1
  E=$2
  M=$3
  compare=$4
  Ns=(10000 80000 100000 160000 320000 640000 800000 1000000 1600000 3200000 6400000 8000000 10000000) # paths
  # Ms=(200) # 1000) # steps in path
  thr=(1 8 16 32 64)
  proc=(1 8 16 32 64)
  hybr=(1 4 8)
  echo "#pragma once
double comparison = ${compare};" > include/comparison.h
  make mc_asia_bin
  for N in ${Ns[@]}
  do
    # for M in ${Ms[@]}
    # do
      ./bin/mc_asia ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results} 
      echo "serial n=${N}, m=${M} -- done"
    # done
  done

  for N in ${Ns[@]} 
  do
    # for M in ${Ms[@]} 
    # do
      for thread in ${thr[@]}
      do
        ./bin/mc_asia_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}
        echo "OMP N=${N}, M=${M}, thread=${thread} -- DONE"
      done
    # done
  done

  for N in ${Ns[@]}
  do
    # for M in ${Ms[@]}
    # do
      for p in ${proc[@]} 
      do
        mpirun -np ${p} --hostfile hostfile ./bin/mc_asia_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}
        echo "MPI N=${N}, M=${M}, processes=${p} -- DONE"
      done
    # done
  done

  for N in ${Ns[@]} 
  do
    # for M in ${Ms[@]} 
    # do
      for hybrid in ${hybr[@]}
      do
        mpirun -np ${hybrid} --hostfile hostfile ./bin/mc_asia_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${hybrid} >> ${results}
        echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
      done
    # done
  done
  echo "MC asia ${payoff_fun} w/ E=${E}: DONE -- N=$N M=$M"
}

# https://www.coggit.com/freetools aritmhetic asian option price

M=(200) # 1000)

C110=(13.71) # 13.74)
# for E in 110
for i in 0 #1
do
  common_cycle "call" 110 ${M[${m}]} ${C110[${i}]} #$(eval "echo \"\$C${E}\"")
done

# Since put and call give analogous results, only call optsion will be calculated from now on (2021-05-16)
# P90=(10.94 10.98)
# # for E in 110
# for i in 0 1
# do
#   common_cycle "put" 90 ${M[${m}]} ${P90[${i}]} #$(eval "echo \"\$P${E}\"")
# done

