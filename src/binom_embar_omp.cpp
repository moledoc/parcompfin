
#include <example_eur.h>

double binom 
(
  double r, double sigma, double S0,
  double T, int N, double E
)
{
  double result;
#pragma omp parallel
  {
  double dt = (double)T/(double)N;
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;
  double tmp;
  int until;
  if (N%2!=0) until = (N+1)/2;
  else until = N/2;
#pragma omp for schedule(dynamic,1000) nowait reduction(+:result)
  for(int i=0;i<until;++i){
    tmp = comb(N,i);
    result += tmp * pow(p,i)*pow(q,N-i)*payoff(S0*pow(u,i)*pow(d,N-i),E);
    result += tmp * pow(p,N-i)*pow(q,i)*payoff(S0*pow(u,N-i)*pow(d,i),E);
  if(i==0&&N%2==0) result+=comb(N,N/2)*pow(p,N/2)*pow(q,N/2)*payoff(S0*pow(u,N/2)*pow(d,N/2),E);
  };
  }
  return exp(-r*T)*result;
}


int main (int argc, char *argv[]){
  int N = getN(argv);
  auto start_overall = std::chrono::system_clock::now();
  int threads = getThreads(argv);
  omp_set_num_threads(threads);
  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "OMP",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N,
      threads
      );
  return EXIT_SUCCESS;
}
