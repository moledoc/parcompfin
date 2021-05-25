
#include <common.h>
#include <comparison.h>

double mc_asia
(
  double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,int M
  ,double payoff_fun
)
{
  double result=0;  

  // start parallel region 
#pragma omp parallel
  {
  double dt = (double)T/(double)M;

  int thr_num = omp_get_thread_num();

  // prepare generator
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time)*(1+thr_num));
  std::normal_distribution<> norm{0,sqrt(dt)};

#pragma omp for nowait schedule(dynamic,1000) reduction(+:result)
  for (int n=0;n<N;++n){
    double St = S0;
    double I = 0;
    // reuse variables I and St, so we could avoid making a vector.
    for (int m=0;m<M;++m){
      double dBi = norm(gen);
      I += St*(1+r*dt/2+sigma*dBi/2);
      St *= exp((r-pow(sigma,2)/2)*dt+sigma*dBi);
    };
    result += payoff(I/(double)M,E,payoff_fun);
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
  int M =                   getArg(argv,8);
  int threads =             getArg(argv,9); 
  omp_set_num_threads(threads);

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  auto start = std::chrono::system_clock::now();
  double result = mc_asia(S0,E,r,sigma,T,N,M,payoff_fun_d);
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
      ,M
      );
  return EXIT_SUCCESS;
}
