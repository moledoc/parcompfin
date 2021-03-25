#!/bin/sh

results=results/results_mc_amer.csv

echo "Method,T_overall,T_calculation,Result,Error,N,M,Parallel" > ${results}
for N in 1 10 100 #1000 # in thousands # paths
do
  for M in 1 10 #100 #1000 # in hundreds # steps in paths
  do
    ./bin/mc_amer ${N}000 ${M}00 >> ${results} 
    for thread in 1 2 4 8 16 32 64 #128
    do
      ./bin/mc_amer_omp ${N}000 ${M}00 ${thread} >> ${results}
    done
    for p in 1 2 4 8 16 32 64 #128
    do
      mpirun -n $p --hostfile hostfile ./bin/mc_amer_mpi ${N}000 ${M}00 >> ${results}
    done
    echo "DONE -- M=$M"
  done
  echo "DONE -- N=$N"
done

