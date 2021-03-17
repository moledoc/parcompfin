
/* #include <mc_example_eur.h> */
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
  int n;// int chunk = 100;
#pragma omp parallel private(n) shared(T,N,r,sigma,E,S0)
{
  /* double local_result=0; */
#pragma omp for nowait schedule(dynamic,10000) reduction(+:result) 
  for(n=0;n<N;++n){
    result += payoff(S0*exp((r-sigma*sigma/2)*T+sigma*norm(gen)),E);
    /* local_result += payoff(S0*exp((r-sigma*sigma/2)*T+sigma*norm(gen)),E)/(double)N; */
  };
/* #pragma omp atomic */
/*     result += local_result; */
}
  return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  int N = getN(argv);
  int threads = getThreads(argv);
  omp_set_num_threads(threads);

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
      N,
      1,
      threads
      );
  return EXIT_SUCCESS;
}
