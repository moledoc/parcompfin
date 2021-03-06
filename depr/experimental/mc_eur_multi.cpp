
#include <common.h>
#include <comparison.h>
#include <mvn.h>

double mc_eur
(
  double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,std::string payoff_fun
  ,int assets
  ,double rho
)
{

  Eigen::MatrixXd samples = sample(N,assets,sigma,rho);

  // payoff of the basket option will depend on the  arithmetic average of prices at maturity T.
  double w_i = 1.0/(double)assets;

  double result=0;
  for(int n=0;n<N;++n){
    double result_n = 0;
    for(int asset=0;asset<assets;++asset){
      // assuming same constant volatility for each underlying asset.
      result_n += S0*exp((r-pow(sigma,2)/2)*T+sigma*samples(asset,n));
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

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(S0,E,r,sigma,T,N,payoff_fun,assets,rho);
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
