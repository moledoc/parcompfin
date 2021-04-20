
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
  double result=0;
  double dt = (double)T/(double)N;
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;
  int until;

  /* for(int i=0;i<=N;++i){ */
  /*   comb_val = comb(N,i); */
  /*   double binom = comb_val + i*log(p) + (N-i)*log(q); */
  /*   result +=  exp(binom) * payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun); */
  /* }; */

  if (N%2!=0) until = (N+1)/2;
  else until = N/2;
  for(int i=0;i<until;++i){
    double comb_val = comb(N,i);
    double binom1 = comb_val + i*log(p) + (N-i)*log(q);
    double binom2 = comb_val + (N-i)*log(p) + i*log(q);
    result +=  exp(binom1) * payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun);
    result +=  exp(binom2) * payoff(S0*pow(u,N-i)*pow(d,i),E,payoff_fun);

    if(i==0 && N%2==0) {
      double binom_mid = comb(N,N/2) + N/2*log(p) + N/2*log(q);
      result+= exp(binom_mid) * payoff(S0*pow(u,N/2)*pow(d,N/2),E,payoff_fun);
    };
  };

  return exp(-r*T)*result;
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
