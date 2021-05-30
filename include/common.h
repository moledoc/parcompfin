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


// make argument into integer
int getArg(char *argv[],int idx){
  std::size_t pos;
  std::string arg = argv[idx];
  int argi = std::stoi(arg,&pos);
  return argi;
}

// make argument into double
double getArgD(char *argv[],int idx){
  /* std::size_t pos; */
  std::string arg = argv[idx];
  double argd = std::stod(arg);
  return argd;
}

/* double payoff_call(double St,double E){ */
/*   if(St-E < 0) return 0; */
/*   else return St-E; */
/* } */

/* double payoff_put(double St,double E){ */
/*   if(E-St < 0) return 0; */
/*   else return E-St; */
/* } */

/* double payoff(double St,double E,std::string payoff_fun){ */
/*   if(payoff_fun == "call") return payoff_call(St,E); */
/*   if(payoff_fun == "put") return payoff_put(St,E); */
/*   if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function"); */
/* } */

double payoff(double St,double E,int payoff_fun){
  //  1 = call option payoff
  // -1 = put option payoff
  return std::max(payoff_fun*(St-E),0.0);
}

// for binom embar
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


// print vector
void vecprinter(std::vector<double> vec)
{
  for(int i = 0;i<vec.size();++i){
    std::cout<<vec[i]<< " " ;
  };
  std::cout<<std::endl;
}


// print matrix
void matprinter(std::vector<std::vector<double>> mat)
{
  for(int i = 0;i<mat.size();++i){
    for(int j = 0;j<mat[i].size();++j){
      std::cout<<mat[i][j] << " " ;
    };
    std::cout<<std::endl;
  };
}

// for mc amer
// in my case it is always 3x3 or 2x2 matrix
std::vector<std::vector<double>> inverse
(
 std::vector<std::vector<double>> x
 ,int size=3
)
{
  std::vector<std::vector<double>> inversed(size);
  for(int i=0;i<size;++i){
    inversed[i].resize(size);
  };
  double determinant=0;
  
  //finding determinant of the matrix
  for(int i=0; i<size;++i){
    determinant += (x[0][i] * (x[1][(i+1)%3] * x[2][(i+2)%3] - x[1][(i+2)%3] * x[2][(i+1)%3]));
  };
  //Condition to check if the derterminat is zero or not if zero than inverse dont exists
  if(determinant<=0){
    throw std::invalid_argument("Detereminant is not > 0");
  };
  for(int i=0;i<size;++i){
    for(int j=0;j<size;++j){
      inversed[j][i] = ((x[(j+1)%3][(i+1)%3] * x[(j+2)%3][(i+2)%3]) - (x[(j+1)%3][(i+2)%3] * x[(j+2)%3][(i+1)%3]))/determinant;
    };
   };
  return inversed;
}

// for mc amer
// matrix/vector multiplicationi function for current solution.
std::vector<double> mat_vec_mul
(
  std::vector<std::vector<double>> x
  ,std::vector<double> y
)
{
  std::vector<double> mat(x.size());
  for(int i=0;i<x.size();++i){
    for(int j=0;j<y.size();++j){
        mat[i]+=x[i][j]*y[j]; 
      };
  };
  return mat;
}

// for mc_amer
// user defined matrix transpose function
std::vector<std::vector<double>> transpose
(
 std::vector<std::vector<double>> y
)
{
  std::vector<std::vector<double>> transposed(y[0].size());
  for(int i=0;i<y[0].size();++i){
    transposed[i].resize(y.size());
  };
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) nowait
  for(int j=0;j<y[0].size();++j){
    for(int i=0;i<y.size();++i){
      transposed[j][i] = y[i][j];
    };
  };
  }
  return transposed;
}

// for mc amer
// user defined matrix transpose function
std::vector<std::vector<double>> pathsfinder
(
 double S0
 ,double E
 ,double r
 ,double sigma
 ,double T
 ,int N
 ,int M
 ,int parallel=0
)
{
  if (N%2!=0) throw std::invalid_argument("N needs to be divisible by 2 for finding paths");
  double dt = T/M;
  // matrix to store paths
  std::vector<std::vector<double>> paths(M+1);
  for(int i=0;i<M+1;++i){
    paths[i].resize(N);
  };
  // prepare generator.
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(dt)};
  gen.seed(time(&cur_time)+100*parallel);

  // generate paths
  for(int n=0;n<N/2;++n){
    // init new path
    paths[0][n] = S0;
    paths[0][n+N/2] = S0;
    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths[m][n] = paths[m-1][n]*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths[m][n+N/2] = paths[m-1][n+N/2]*exp((r-0.5*sigma*sigma)*dt-sigma*w);

    };
  };
  return paths;
}


// output calculation results
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
