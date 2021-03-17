#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>

int getN(char *argv[]){
  std::size_t pos;
  std::string arg = argv[1];
  int N = std::stoi(arg,&pos);
  return N;
}

int getThreads(char *argv[]){
  std::size_t pos;
  std::string arg = argv[2];
  int threads = std::stoi(arg,&pos);
  return threads;
}

void reporting
(
    double time_overall,
    double time,
    double result,
    double analytical,
    int N,
    int processes=1,
    int threads=1
    )
{
  /* std::cout << "time_overall,time_calculation,result,N,error,P,script" << std::endl; */
  std::cout << std::setprecision(10) \
    << time_overall << "," \
    << time << "," \
    << result << "," \
    /* << abs(result-analytical) << "," \ */ 
    << N << "," \ 
    << processes << "," \ 
    << threads \
    << std::endl;
}
