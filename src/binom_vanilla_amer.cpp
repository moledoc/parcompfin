
#include <example_amer.h>

double binom (
    double r, double sigma, double S0,
    double T, int N, double E
    )
{
  double dt = (double)T/(double)N;
  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta+1); */
  /* double d = beta - sqrt(beta*beta+1); */
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;

  std::vector<double> v_ij;//(N+1,0);

  for(int i=0;i<N;++i)
    v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E));
    /* v_ij[i] = payoff(S0*pow(u,i)*pow(d,N-i),E); */

  for (int n=N-1;n>=0;--n){
    for(int i=0;i<n+1;++i){
      double sij = payoff(S0*pow(u,i)*pow(d,n-i),E);
      double jatk = (p*v_ij[i+1] + q*v_ij[i])/R;
      v_ij[i] = std::max(jatk,sij);
    };
  };

  return v_ij[0];
}


int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N
      );
  return EXIT_SUCCESS;
}
