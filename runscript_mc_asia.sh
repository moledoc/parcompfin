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
  compare=$3
  Ns=(1000 2500 5000 7500 10000 25000 50000 75000 100000 250000 500000 750000 1000000)  # paths
  Ms=(200 1000) # steps in paths
  thr=(1 5 10 25 32 50 64 100 125)
  proc=(1 5 10 25 32 50 64 100 125)
  hybr=(1 5 16 25 32 50 60)
  echo "#pragma once
double comparison = ${compare};" > include/comparison.h
  make mc_asia_bin
  for N in ${Ns[@]}
  do
    for M in ${Ms[@]}
    do
      ./bin/mc_asia ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results} 
      echo "Serial N=${N}, M=${M} -- DONE"
    done
  done

  for N in ${Ns[@]} 
  do
    for M in ${Ms[@]} 
    do
      for thread in ${thr[@]}
      do
        ./bin/mc_asia_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}
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
        mpirun -np ${p} --hostfile hostfile ./bin/mc_asia_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}
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
        mpirun -np ${hybrid} --hostfile hostfile ./bin/mc_asia_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${hybrid} >> ${results}
        echo "Hybrid N=${N}, processes=${hybrid}, thread=${hybrid} -- DONE"
      done
    done
  done
  echo "MC asia ${payoff_fun} w/ E=${E}: DONE -- N=$N"
}

# https://www.coggit.com/freetools aritmhetic asian option price

C90=21.88
C95=19.53
C100=17.43
C105=15.54
C110=13.86
C120=11.02
C130=8.77


for E in 110 #90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E} $(eval "echo \"\$C${E}\"")
done

P60=1.49
P70=3.55
P80=6.75
P90=11.08
P95=13.64
P100=16.43
P105=19.45
P110=22.66

for E in 90 #110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E} $(eval "echo \"\$P${E}\"")
done

