
#include <common.h>
#include <comparison.h>
#include <mvn.h>

double mc_eur
(
  double S0 ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,double payoff_fun
  ,int assets
  ,double rho
)
{

  // pre-calculate random variables from N(Mu,Sigma);
  Eigen::MatrixXd Bt = mvnorm(N,assets,rho);

  // payoff of the basket option will depend on the  arithmetic average of prices at maturity T.
  double w_i = 1.0/(double)assets;

  double result=0;
  for(int n=0;n<N;++n){
    double result_n = 0;
    // calculate basket price.
    for(int asset=0;asset<assets;++asset){
      result_n += w_i*S0*exp((r-pow(sigma,2)/2)*T+sigma*Bt(asset,n));
    };
    result += payoff(result_n,E,payoff_fun);
  };
  return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  std::string payoff_fun =  argv[1];
  double S0 =               getArgD(argv,2);
  double E =                getArgD(argv,3);
  double r =                getArgD(argv,4);
  /* double sigma =       sqrt(getArgD(argv,5)); */
  double sigma =            getArgD(argv,5);
  double T =                getArgD(argv,6);
  int N =                   getArg(argv,7);
  int assets =              getArg(argv,8);
  double rho =              getArgD(argv,9);

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(S0,E,r,sigma,T,N,payoff_fun_d,assets,rho);
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
      ,0
      ,0
      ,assets
      );
  return EXIT_SUCCESS;
}
