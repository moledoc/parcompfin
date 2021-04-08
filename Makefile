# CXX=clang++
# CXX=g++
CXX=mpic++
CXX_MPI=mpic++
# Set default compiler parameters
# -Wall 	shows all warnings when compiling, always use this!
# -std=c++11 	enables the C++17 standard mode
# CXXFLAGS = -Wall -std=c++17 -Iinclude
CXXFLAGS = -Wall -std=c++17 -Iinclude -Ofast

########################################################################################################################

all: init binom mc_eur mc_amer
all_bin: init binom_bin mc_eur_bin  mc_amer_bin 

init:
		if [ ! -d bin ]; then mkdir bin; fi
		if [ ! -d obj ]; then mkdir obj; fi
		if [ ! -d results ]; then mkdir results; fi
		echo "#pragma once" > include/comparison.h; 
		echo "double comparison=0;" >> include/comparison.h

########################################################################################################################

binom_prep: include/example_eur.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/binom_embar.cpp -o obj/binom_embar.o
	$(CXX) $(CXXFLAGS) -c src/binom_embar_omp.cpp -o obj/binom_embar_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/binom_embar_mpi.cpp -o obj/binom_embar_mpi.o
	$(CXX_MPI) $(CXXFLAGS) -c src/binom_embar_hybrid.cpp -o obj/binom_embar_hybrid.o -fopenmp

binom_bin: binom_prep
	$(CXX) $(CXXFLAGS) obj/binom_embar.o -o bin/binom_embar 
	$(CXX) $(CXXFLAGS) obj/binom_embar_omp.o -o bin/binom_embar_omp -fopenmp 
	$(CXX_MPI) $(CXXFLAGS) obj/binom_embar_mpi.o -o bin/binom_embar_mpi
	$(CXX_MPI) $(CXXFLAGS) obj/binom_embar_hybrid.o -o bin/binom_embar_hybrid -fopenmp 

binom_tst: binom_bin
	./bin/binom_embar call 100 95 0.02 0.75 1 1000
	./bin/binom_embar_omp call 100 95 0.02 0.75 1 1000 2
	mpirun -n 4 --hostfile hostfile ./bin/binom_embar_mpi call 100 95 0.02 0.75 1 1000
	mpirun -n 4 --hostfile hostfile ./bin/binom_embar_hybrid call 100 95 0.02 0.75 1 1000 2

binom: init #binom_bin
	#./runscript_binom_embar.sh
	dash runscript_binom_embar.sh


########################################################################################################################

mc_eur_prep: include/example_eur.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/mc_eur.cpp -o obj/mc_eur.o
	$(CXX) $(CXXFLAGS) -c src/mc_eur_omp.cpp -o obj/mc_eur_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_eur_mpi.cpp -o obj/mc_eur_mpi.o
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_eur_hybrid.cpp -o obj/mc_eur_hybrid.o -fopenmp

mc_eur_bin: mc_eur_prep
	$(CXX) $(CXXFLAGS) obj/mc_eur.o -o bin/mc_eur 
	$(CXX) $(CXXFLAGS) obj/mc_eur_omp.o -o bin/mc_eur_omp -fopenmp
	$(CXX_MPI) $(CXXFLAGS) obj/mc_eur_mpi.o -o bin/mc_eur_mpi
	$(CXX_MPI) $(CXXFLAGS) obj/mc_eur_hybrid.o -o bin/mc_eur_hybrid -fopenmp

mc_eur_tst: mc_eur_bin
	./bin/mc_eur 10000000
	./bin/mc_eur_omp 10000000 4
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_mpi 10000000
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_hybrid 10000000 4

mc_eur: init #mc_eur_bin
	./runscript_mc_eur.sh

########################################################################################################################

mc_amer_prep: include/example_amer.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/binom_vanilla_amer.cpp -o obj/binom_vanilla_amer.o # test american option value
	$(CXX) $(CXXFLAGS) -c src/mc_amer.cpp -o obj/mc_amer.o
	$(CXX) $(CXXFLAGS) -c src/mc_amer_omp.cpp -o obj/mc_amer_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_amer_mpi.cpp -o obj/mc_amer_mpi.o
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_amer_hybrid.cpp -o obj/mc_amer_hybrid.o -fopenmp

mc_amer_bin: mc_amer_prep
	$(CXX) $(CXXFLAGS) obj/binom_vanilla_amer.o -o bin/binom_vanilla_amer # test american option value
	$(CXX) $(CXXFLAGS) obj/mc_amer.o -o bin/mc_amer 
	$(CXX) $(CXXFLAGS) obj/mc_amer_omp.o -o bin/mc_amer_omp -fopenmp
	$(CXX_MPI) $(CXXFLAGS) obj/mc_amer_mpi.o -o bin/mc_amer_mpi
	$(CXX_MPI) $(CXXFLAGS) obj/mc_amer_hybrid.o -o bin/mc_amer_hybrid -fopenmp

mc_amer_tst: mc_amer_bin
	./bin/binom_vanilla_amer 1000
	./bin/mc_amer 10000 10000
	./bin/mc_amer_omp 10000 10000 4
	mpirun -n 4 --hostfile hostfile ./bin/mc_amer_mpi 10000 10000
	mpirun -n 4 --hostfile hostfile ./bin/mc_amer_hybrid 10000 10000 4

mc_amer: init #mc_amer_bin
	./runscript_mc_amer.sh

########################################################################################################################

clean:
	rm -rf bin obj #results
