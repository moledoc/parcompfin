# !/bin/sh

results=results/results_mc_amer.csv

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
  Ns=(10000 25000 32000 50000 64000 100000 250000 320000 500000 640000 1000000 2500000 3200000 5000000 6400000 10000000) # paths
  Ms=(200) # 1000) # steps in path
  thr=(1 5 10 25 32 50 64 100)
  proc=(1 5 10 25 32 50 64 100)
  hybr=(1 5 16 25 32 50)
  echo "#pragma once
double comparison = $(./bin/binom_vanilla_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 10000 | tr ',' '\t' | awk '{print $14}');" > include/comparison.h
  make mc_amer_bin
  for N in ${Ns[@]} 
  do
    for M in ${Ms[@]} 
    do
      ./bin/mc_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results} 
      echo "Serial N=${N} M=${M} -- DONE"
    done
  done

  for N in  ${Ns[@]}
  do
    for M in ${Ms[@]} 
    do
      for thread in ${thr[@]} 
      do
        ./bin/mc_amer_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}
        echo "OMP N=${N}, M=${M}, thread=${thread} -- DONE"
      done
    done
  done

  for N in ${Ns[@]}
  do
    for M in ${Ms[@]}
    do
      for p in ${proc[@]}
      do
        mpirun -np ${p} --hostfile hostfile ./bin/mc_amer_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}
        echo "MPI N=${N}, M=${M}, processes=${p} -- DONE"
      done
    done
  done

  for N in ${Ns[@]}
  do
    for M in ${Ms[@]}
    do
      for hybrid in ${hybr[@]}
      do
        mpirun -np ${hybrid} --hostfile hostfile ./bin/mc_amer_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${hybrid} >> ${results}
        echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
      done
    done
  done

  echo "MC amer ${payoff_fun} w/ E=${E}: DONE -- N=$N"
}

for E in 110 #90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

# Since put and call give analogous results, only call optsion will be calculated from now on (2021-05-16)
# for E in 90 #110 105 100 95 90 80 70 60
# do
#   common_cycle "put" ${E}
# done

