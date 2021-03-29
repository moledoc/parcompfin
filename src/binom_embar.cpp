
#include <example_eur.h>

double binom (
    double r, double sigma, double S0,
    double T, int N, double E
    )
{
  double dt = (double)T/(double)N;
  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta-1); */
  /* double d = beta - sqrt(beta*beta-1); */
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;

  double result;

  int until;
  if (N%2!=0) until = (N+1)/2;
  else until = N/2;
  if(N%2==0) result+=exp(comb(N,N/2) + log(pow(p,N/2)*pow(q,N/2)*payoff(S0*pow(u,N/2)*pow(d,N/2),E)));
  for(int i=0;i<until;++i){
    double tmp = comb(N,i);
    /* std::cout << tmp << std::endl; */
    result += exp(tmp + log(pow(p,i)*pow(q,N-i)*payoff(S0*pow(u,i)*pow(d,N-i),E)));
    result += exp(tmp + log(pow(p,N-i)*pow(q,i)*payoff(S0*pow(u,N-i)*pow(d,i),E)));
  };
  /* for(int i=0;i<N;++i) */
  /*   result += comb(N,i)* pow(p,i)*pow(q,N-i)*payoff(S0*pow(u,i)*pow(d,N-i),E); */

  return exp(-r*T)*result;
}


int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start;
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
