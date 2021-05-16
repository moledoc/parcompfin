
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
  double dt = (double)T/(double)M;

  time_t cur_time1;
  std::random_device rd1{};
  std::mt19937 gen1{rd1()};
  gen1.seed(time(&cur_time1));
  std::normal_distribution<> norm1{0,sqrt(dt)};

  time_t cur_time2;
  std::random_device rd2{};
  std::mt19937 gen2{rd2()};
  gen2.seed(time(&cur_time2));
  std::normal_distribution<> norm2{0,sqrt(dt)};

  double result=0;  

  for (int n=0;n<N;++n){
    double St = S0;
    double I = 0;

    for (int m=0;m<M;++m){
      I += St*(1+r*dt/2+sigma*norm2(gen2)/2);
      St *= exp((r-pow(sigma,2)/2)*dt+sigma*norm1(gen1));
    };
    result += payoff(I/(double)M,E,payoff_fun);
  };

  return exp(-r*T)*result/(double)N;
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
