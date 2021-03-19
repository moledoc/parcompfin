#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

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
    std::string method,
    double time_overall,
    double time,
    double result,
    double analytical,
    int N,
    int parallel=0
    )
{
  /* std::cout << "time_overall,time_calculation,result,N,error,P,script" << std::endl; */
  std::cout << std::setprecision(10) \
    << method << "," \
    << time_overall << "," \
    << time << "," \
    << result << "," \
    << abs(result-analytical) << "," \ 
    << N << "," \ 
    << parallel \
    << std::endl;
}
