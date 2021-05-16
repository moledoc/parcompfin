
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
  double dt = (double)T/(double)N;
  double beta = 0.5*(exp(-r*dt)+exp((r+pow(sigma,2))*dt));
  double u = beta + sqrt(pow(beta,2)-1);
  double d = beta - sqrt(pow(beta,2)-1);
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  /* double u = exp(sigma*sqrt(dt)); */
  /* double d = 1/u; */
  double q = 1-p;

  std::vector<double> v_ij;

  for(int i=0;i<N+1;++i)
    v_ij.push_back(std::max(((S0*pow(u,i)*pow(d,N-i))-E)*payoff_fun,(double)0.0));
    /* v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun)); */

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

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(S0,E,r,sigma,T,N,payoff_fun_d);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial_vanilla"
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
