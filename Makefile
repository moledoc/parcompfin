# CXX=clang++
CXX=g++
CXX_MPI=mpic++
# Set default compiler parameters
# -Wall 	shows all warnings when compiling, always use this!
# -std=c++11 	enables the C++11 standard mode
CXXFLAGS = -Wall -std=c++11 -Iinclude

all: init mc

init:
		if [ ! -d bin ]; then mkdir bin; fi
		if [ ! -d obj ]; then mkdir obj; fi
		if [ ! -d results ]; then mkdir results; fi

mc_prep: include/example_eur.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/mc_eur.cpp -o obj/mc_eur.o
	$(CXX) $(CXXFLAGS) -c src/mc_eur_omp.cpp -o obj/mc_eur_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_eur_mpi.cpp -o obj/mc_eur_mpi.o
	$(CXX_MPI) $(CXXFLAGS) -c src/mc_eur_hybrid.cpp -o obj/mc_eur_hybrid.o -fopenmp

mc_bin: mc_prep
	$(CXX) $(CXXFLAGS) obj/mc_eur.o -o bin/mc_eur 
	$(CXX) $(CXXFLAGS) obj/mc_eur_omp.o -o bin/mc_eur_omp -fopenmp
	$(CXX_MPI) $(CXXFLAGS) obj/mc_eur_mpi.o -o bin/mc_eur_mpi
	$(CXX_MPI) $(CXXFLAGS) obj/mc_eur_hybrid.o -o bin/mc_eur_hybrid -fopenmp

mc_tst: mc_bin
	./bin/mc_eur 10000000
	./bin/mc_eur_omp 10000000 4
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_mpi 10000000
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_hybrid 10000000 4

mc: init mc_bin
	./runscript_mc_eur.sh

mc_amer_prep: include/example_amer.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/mc_amer.cpp -o obj/mc_amer.o

mc_amer: mc_amer_prep
	$(CXX) $(CXXFLAGS) obj/mc_amer.o -o bin/mc_amer 


binom_prep: include/example_eur.h include/reporting.h
	$(CXX) $(CXXFLAGS) -c src/binom_vanilla.cpp -o obj/binom_vanilla.o
	$(CXX) $(CXXFLAGS) -c src/binom_vanilla_omp.cpp -o obj/binom_vanilla_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/binom_vanilla_mpi.cpp -o obj/binom_vanilla_mpi.o
	# $(CXX_MPI) $(CXXFLAGS) -c src/binom_vanilla_hybrid.cpp -o obj/binom_vanilla_hybrid.o -fopenmp
	$(CXX) $(CXXFLAGS) -c src/binom_embar.cpp -o obj/binom_embar.o
	$(CXX) $(CXXFLAGS) -c src/binom_embar_omp.cpp -o obj/binom_embar_omp.o -fopenmp
	$(CXX_MPI) $(CXXFLAGS) -c src/binom_embar_mpi.cpp -o obj/binom_embar_mpi.o

binom_bin: binom_prep
	$(CXX) $(CXXFLAGS) obj/binom_vanilla.o -o bin/binom_vanilla 
	$(CXX) $(CXXFLAGS) obj/binom_vanilla_omp.o -o bin/binom_vanilla_omp -fopenmp
	$(CXX_MPI) $(CXXFLAGS) obj/binom_vanilla_mpi.o -o bin/binom_vanilla_mpi
	# $(CXX_MPI) $(CXXFLAGS) obj/binom_vanilla_hybrid.o -o bin/binom_vanilla_hybrid -fopenmp
	$(CXX) $(CXXFLAGS) obj/binom_embar.o -o bin/binom_embar 
	$(CXX) $(CXXFLAGS) obj/binom_embar_omp.o -o bin/binom_embar_omp -fopenmp 
	$(CXX_MPI) $(CXXFLAGS) obj/binom_embar_mpi.o -o bin/binom_embar_mpi

binom_tst: binom_bin
	./bin/binom_vanilla 1000
	./bin/binom_vanilla_omp 1000 2
	# mpirun -n 4 --hostfile hostfile ./bin/binom_vanilla_mpi 10
	# mpirun -n 4 --hostfile hostfile ./bin/binom_vanilla_hybrid 50000 8
	./bin/binom_embar 1000
	./bin/binom_embar_omp 1000 2
	mpirun -n 4 --hostfile hostfile ./bin/binom_embar_mpi 1000

binom: init binom_bin
	# ./runscript_binom_vanilla.sh

clean:
	rm -rf bin obj #results
