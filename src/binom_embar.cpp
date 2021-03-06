
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
  ,double payoff_fun
)
{
  // initialize result variable
  double result=0;
  // calculate parameters
  double dt = (double)T/(double)N;
  double beta = 0.5*(exp(-r*dt)+exp((r+pow(sigma,2))*dt));
  double u = beta + sqrt(pow(beta,2)-1);
  double d = beta - sqrt(pow(beta,2)-1);
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  /* double u = exp(sigma*sqrt(dt)); */
  /* double d = 1/u; */
  double q = 1-p;

  // We will reuse calculated combination value.
  // Handle the iterationi limit for combination value reuse.
  int until;
  if (N%2!=0) until = (N+1)/2;
  else until = N/2;
  for(int i=0;i<until;++i){
    double comb_val = comb(N,i);
    double binom1 = comb_val + i*log(p) + (N-i)*log(q);
    double binom2 = comb_val + (N-i)*log(p) + i*log(q);
    result +=  exp(binom1) * payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun);
    result +=  exp(binom2) * payoff(S0*pow(u,N-i)*pow(d,i),E,payoff_fun);

    // Handle the middle combination value, when mod(N,2)==0.
    if(i==0 && N%2==0) {
      double binom_mid = comb(N,N/2) + N/2*log(p) + N/2*log(q);
      result+= exp(binom_mid) * payoff(S0*pow(u,N/2)*pow(d,N/2),E,payoff_fun);
    };
  };

  // discount the result.
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

  double payoff_fun_d;
  if(payoff_fun=="put") payoff_fun_d=-1;
  if(payoff_fun=="call") payoff_fun_d=1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(S0,E,r,sigma,T,N,payoff_fun_d);
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
