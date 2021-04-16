# !/bin/sh

results=results/results_mc_amer.csv

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
double comparison = $(./bin/binom_vanilla_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} 1000 | tr ',' '\t' | awk '{print $13}');" > include/comparison.h
  make mc_amer_bin
  for N in  1024 8192 16384 65536 131072 524288 1048276 #1000 5000 10000 50000 #100000 #1000000 # paths
  do
    for M in 50 100 250 500 750 1000 10000 # steps in paths
    do
      ./bin/mc_amer ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results} 
      echo "Serial N=${N} M=${M} -- DONE"
      for thread in 1 2 4 #8 #16 32 64 #128
      do
        ./bin/mc_amer_omp ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} ${thread} >> ${results}
        echo "OMP N=${N}, M=${M}, thread=${thread} -- DONE"
      done
      for p in 1 2 4 8 16 32 64 #128
      do
        mpirun -n $p --hostfile hostfile ./bin/mc_amer_mpi ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} >> ${results}
        echo "MPI N=${N}, M=${M}, processes=${p} -- DONE"
        if [ ${p} == 1 ];
        then
          mpirun -n $p --hostfile hostfile ./bin/mc_amer_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} 1 >> ${results}
        elif [ ${p} -le 8 ];
        then
          for thread in 2 4 8 #16 32 64 #128
          do
            mpirun -n $p --hostfile hostfile ./bin/mc_amer_hybrid ${payoff_fun} ${S0} ${E} ${r} ${sigma} ${T} ${N} ${M} $thread >> ${results}
            ecHo "hybrid N=${N}, M=${M}, processes=${p}, thread=${thread} -- DONE"
          done
        fi
      done
      echo "MC amer ${payoff_fun} w/ E=${E}: DONE -- M=$M"
    done
    echo "MC amer ${payoff_fun} w/ E=${E}: DONE -- N=$N"
  done
}

for E in 90 95 100 105 110 120 130 140
do
  common_cycle "call" ${E}
done

for E in 110 105 100 95 90 80 70 60
do
  common_cycle "put" ${E}
done

