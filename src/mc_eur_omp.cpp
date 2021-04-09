
#include <common.h>
#include <comparison.h>

double mc_eur
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
#pragma omp parallel 
{
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time));
  std::normal_distribution<> norm{0,sqrt(T)};
  /* double local_result=0; */
#pragma omp for nowait schedule(dynamic,1000) reduction(+:result) 
  for(int n=0;n<N;++n){
    result += payoff(S0*exp((r-pow(sigma,2)/2)*T+sigma*norm(gen)),E,payoff_fun);
  };
}
  return (exp(-r*T)*result)/(double)N;
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
  int threads =             getArg(argv,8); 
  omp_set_num_threads(threads);

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(S0,E,r,sigma,T,N,payoff_fun);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "OMP"
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
      ,threads
      );
  return EXIT_SUCCESS;
}
