
#include <example_eur.h>

double binom (
    double r, double sigma, double S0,
    double T, int N, double E
    )
{
  double result;
#pragma omp parallel 
    {
  std::vector<double> v_ij(N);
  double dt = T/(double)N;
  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta-1); */
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;//beta - sqrt(beta*beta-1);
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;
#pragma omp for schedule(dynamic,1000) nowait
  for(int i=0;i<N;++i){
    v_ij[i] = payoff(S0*pow(u,i)*pow(d,N-i),E);
  };

  for (int n=N-1;n>=0;--n){
#pragma omp for schedule(dynamic,1000) nowait
    for(int i=0;i<n+1;++i){
      v_ij[i] = (p*v_ij[i+1] + q*v_ij[i])/R;
    };
    };
  result = v_ij[0];
  }
  return result;
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
