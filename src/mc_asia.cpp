
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
  ,std::string payoff_fun
)
{
  double dt = (double)T/(double)M;

  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time));
  int unif_prec = 100000000;
  std::uniform_int_distribution<> unif{1,unif_prec};


  time_t cur_time2;
  std::random_device rd2{};
  std::mt19937 gen2{rd2()};
  gen.seed(time(&cur_time2));
  std::normal_distribution<> norm{0,sqrt(dt)};

  double result=0;  

  for (int n=0;n<N;++n){
    double St = S0;
    double I = 0;

    for (int m=0;m<M;++m){
      double x = (double)unif(gen)/(double)unif_prec;
      double y = (double)unif(gen)/(double)unif_prec;
      double v = sqrt(-2*log(x))*sin(2*M_PI*y);

      I += St*(1+r*dt/2+sigma*norm(gen2)/2);
      St *= exp((r-pow(sigma,2)/2)*dt+sigma*v);

    };
    result += payoff(I/(double)M,E,payoff_fun);
  };

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

  auto start = std::chrono::system_clock::now();
  double result = mc_asia(S0,E,r,sigma,T,N,M,payoff_fun);
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
      ,M
      );
  return EXIT_SUCCESS;
}
