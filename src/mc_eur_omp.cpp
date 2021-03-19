
/* #include <mc_example_eur.h> */
#include <example_eur.h>

double mc_eur
 (
 int N,double S0,double E,double r, double T,double sigma
 )
{
  double result=0;
  /* int chunk = 100; */
#pragma omp parallel //private(T,r,sigma,E,S0,norm,gen)
{
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time));
  std::normal_distribution<> norm{0,sqrt(T)};
  /* double local_result=0; */
#pragma omp for nowait schedule(dynamic,1000) reduction(+:result) 
  for(int n=0;n<N;++n){
    result += payoff(S0*exp((r-pow(sigma,2)/2)*T+sigma*norm(gen)),E);
    /* local_result += payoff(S0*exp((r-sigma*sigma/2)*T+sigma*norm(gen)),E)/(double)N; */
  };
/* #pragma omp atomic */
/*     result += local_result; */
}
  return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  int N = getN(argv);
  auto start_overall = std::chrono::system_clock::now();
  int threads = getThreads(argv);
  omp_set_num_threads(threads);

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(N,S0,E,r,T,sigma);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "OMP",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N,
      threads
      );
  return EXIT_SUCCESS;
}
