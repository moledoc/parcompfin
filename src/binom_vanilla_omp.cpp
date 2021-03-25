
#include <example_eur.h>

double binom (
    double r, double sigma, double S0,
    double T, int N, double E
    )
{
  double result;
  std::vector<double> v_ij(N);
  double dt = T/(double)N;
  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta-1); */
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;//beta - sqrt(beta*beta-1);
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;
#pragma omp parallel //private(dt,u,d,R,p,q,E,S0)
    {
#pragma omp for schedule(dynamic,1000) nowait
  for(int i=0;i<N;++i){
    v_ij[i] = payoff(S0*pow(u,i)*pow(d,N-i),E);
  };
    }

  std::vector<double> tmp(N);
  for (int n=N-1;n>=0;--n){
#pragma omp parallel //private(dt,u,d,R,p,q,E,S0)
    {
#pragma omp for schedule(dynamic,10) nowait
    for(int i=0;i<n+1;++i){
      /* v_ij[i] = (p*v_ij[i+1] + q*v_ij[i])/R; */
      tmp[i] = (p*v_ij[i+1] + q*v_ij[i])/R;
      };
/* #pragma omp for schedule(dynamic,1000) nowait */
/*     for(int i=0;i<n+1;++i){ */
/*       v_ij[i] = tmp[i]; */
/*     }; */
    v_ij.clear();
    copy(tmp.begin(),tmp.end(),back_inserter(v_ij));
    tmp.clear();
    }

    };
  /* result = v_ij[0]; */
  /* return result; */
    return v_ij[0];
}



int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int threads = getArg(argv,2);
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
