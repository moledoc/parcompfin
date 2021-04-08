#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

int getArg(char *argv[],int idx){
  std::size_t pos;
  std::string arg = argv[idx];
  int argi = std::stoi(arg,&pos);
  return argi;
}

double getArgD(char *argv[],int idx){
  std::size_t pos;
  std::string arg = argv[idx];
  int argd = std::stod(arg,&pos);
  return argd;
}

void reporting
(
    std::string method,
    double time_overall,
    double time,
    double result,
    double analytical,
    int N,
    int M=0,
    int parallel=0
    )
{
  /* std::cout << "time_overall,time_calculation,result,N,error,P,script" << std::endl; */
  std::cout << std::setprecision(10) \
    << method << "," \
    << time_overall << "," \
    << time << "," \
    << result << "," \
    /* << abs(result-analytical) << "," \ */ 
    << result-analytical << "," \ 
    << N << "," \ 
    << M << "," \ 
    << parallel \
    << std::endl;
}
