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
double analytical = 61.5; // not proven number

double payoff(double ST,double E){
  /* european call */
  if(ST-E < 0) return 0;
  else return ST-E;
}

double comb(int N,int i){
  if (i==0 || i==N) return 1;
  if (i==1 || i==(N-1)) return N;
  double result=1;
  int until;
  if (i >= (N-i)) until = N-i;
  else until = i;
  for(int j=1;j<=until;++j){
    result *= (double)(N-j+1)/(double)(until-j+1);
    /* result *= (double)(N-j+1)/(double)j; */
  };
  return result;
}
