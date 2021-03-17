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
	./bin/mc_eur_omp 10000000 8
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_mpi 10000000
	mpirun -n 4 --hostfile hostfile ./bin/mc_eur_hybrid 10000000 8

mc: init mc_bin
	./runscript_mc_eur

clean:
	rm -rf bin obj results
