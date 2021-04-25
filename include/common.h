#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>
#include <ctime>
#include <vector>
#include <math.h>
#include <random>
#include <omp.h>
#include <mpi.h>
#include <stdexcept> // exceptions
#include <algorithm>
#include <iterator>
#include <unordered_map>

#define _USE_MATH_DEFINES
#include<cmath>

#include <comparison.h>


int getArg(char *argv[],int idx){
  std::size_t pos;
  std::string arg = argv[idx];
  int argi = std::stoi(arg,&pos);
  return argi;
}

double getArgD(char *argv[],int idx){
  /* std::size_t pos; */
  std::string arg = argv[idx];
  double argd = std::stod(arg);
  return argd;
}

double payoff_call(double St,double E){
  if(St-E < 0) return 0;
  else return St-E;
}

double payoff_put(double St,double E){
  if(E-St < 0) return 0;
  else return E-St;
}

double payoff(double St,double E,std::string payoff_fun){
  if(payoff_fun == "call") return payoff_call(St,E);
  if(payoff_fun == "put") return payoff_put(St,E);
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");
}

/* double payoff(double St,double E,int payoff_fun){ */
/*   return std::max(payoff_fun*(St-E),0); */
/* } */

double comb(int N,int i){
  if (i==0 || i==N) return 0;
  if (i==1 || i==(N-1)) return log(N);
  double result=0;
  for(int n=N;n>i;--n)
    result += log((double)n);
  for(int j=2;j<=(N-i);++j)
    result -= log((double)j);
  return result;
}

void reporting
(
  std::string method
  ,std::string payoff_fun
  ,double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,double time_overall
  ,double time
  ,double result
  ,double comparison
  ,int N
  ,int parallel=0
  ,int M=0
  ,int assets=1
)
{
  std::cout << std::setprecision(10) \
    << method << "," \
    << payoff_fun << "," \
    << S0 << "," \
    << E << "," \
    << r << "," \
    << sigma << "," \
    << T << "," \
    << N << "," \ 
    << M << "," \ 
    << parallel << ","\
    << assets << ","\
    << time_overall << "," \
    << time << "," \
    << result << "," \
    << abs(result-comparison) << "," \ 
    << result-comparison/* << ","*/ \ 
    << std::endl;
}
