
#include <common.h>
#include <comparison.h>

double binom 
(
  double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,std::string payoff_fun
)
{
  double dt = (double)T/(double)N;
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;

  std::vector<double> v_ij;//(N+1,0);

  for(int i=0;i<N;++i)
    v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun));

  for (int n=N-1;n>=0;--n){
    for(int i=0;i<n+1;++i){
      v_ij[i] = (p*v_ij[i+1] + q*v_ij[i])/R;
    };
  };

  return v_ij[0];
}


int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  std::string payoff_fun =  argv[1];
  double S0 =               getArgD(argv,2);
  double E =                getArgD(argv,3);
  double r =                getArgD(argv,4);
  double sigma =            getArgD(argv,5);
  double T =                getArgD(argv,6);
  int N =                   getArg(argv,7);
  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(S0,E,r,sigma,T,N,payoff_fun);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial"
      ,payoff_fun
      ,S0
      ,E
      ,r
      ,sigma
      ,T
      ,elapsed_seconds_overall.count()
      ,elapsed_seconds.count()
      ,result
      ,comparison
      ,N
      );
  return EXIT_SUCCESS;
}
