
#include <example_eur.h>

double mc_eur
 (
 int N,double S0,double E,double r, double T,double sigma
 )
{
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time));
  std::normal_distribution<> norm{0,sqrt(T)};

  double result=0;
  for(int n=0;n<N;++n){
    result += payoff(S0*exp((r-sigma*sigma/2)*T+sigma*norm(gen)),E);
  };
  return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  int N = getN(argv);

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(N,S0,E,r,T,sigma);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start;
  reporting(
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N
      );
  return EXIT_SUCCESS;
}
