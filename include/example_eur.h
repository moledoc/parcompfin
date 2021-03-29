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
double analytical = 61.47789; 

double payoff(double ST,double E){
  /* european call */
  if(ST-E < 0) return 0;
  else return ST-E;
}

double comb(int N,int i){
  if (i==0 || i==N) return log(1);
  if (i==1 || i==(N-1)) return log(N);
  double result=0;
  /* int until; */
  /* if (i > (N-i)) until = N-i; */
  /* else until = i; */
  /* for(int j=1;j<=until;++j){ */
  /*   result += log((double)(N-j+1))-log((double)(until-j+1)); */
  /*   /1* result *= (double)(N-j+1)/(double)j; *1/ */
  /* }; */

  for(int n=N;n>i;--n)
    result += log((double)n);

  /* if (i==(N-i)) until = i+1; */
  for(int j=2;j<=(N-i);++j)
    result -= log((double)j);

    /* result *= (double)(N-j+1)/(double)j; */
  return result;
}
