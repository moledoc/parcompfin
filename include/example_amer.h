#pragma once
#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include <math.h>
#include <random>
#include <omp.h>
#include <mpi.h>
#include <iomanip>
#include <stdexcept> // exceptions
#include <algorithm>
#include <iterator>
#include <unordered_map>

#include <reporting.h>

double S0 = 100;
double E = 40;
double r = 0.05;
double sigma = 0.75;
double T = 0.5;
double analytical = 0;// TODO: //61.47789; 

double payoff(double St,double E){
  /* american call */
  if(St-E < 0) return 0;
  else return St-E;
}
