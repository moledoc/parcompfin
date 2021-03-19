#!/bin/sh

results=results/results_mc_eur.csv

echo "Method,T_overall,T_calculation,Result,Error,N,Parallel" > ${results}
for N in 1 2 5 10 25 50 100 #1 10 100 #
do
 	./bin/mc_eur ${N}000000 >> ${results} 
  for thread in 1 2 4 8 16 32 64 #128
  do
    ./bin/mc_eur_omp ${N}000000 $thread >> ${results}
  done
  for p in 1 2 4 8 16 32 64 #128
  do
    mpirun -n $p --hostfile hostfile ./bin/mc_eur_mpi ${N}000000 >> ${results}
    # for thread in 2 4 8 16
    # do
    #   mpirun -n $p --hostfile hostfile ./bin/mc_eur_hybrid ${N}000000 $thread >> ${results}
    # done
  done
  echo "DONE -- $N"
done

