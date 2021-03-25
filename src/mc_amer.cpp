
#include <example_amer.h>

double mc_amer
 (
 int N,int M,double S0,double E,double r, double T,double sigma
 )
{
  double result=0;
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time));
  std::normal_distribution<> norm{0,sqrt(T)};

  /* std::vector<double> walk; */
  /* std::vector<std::vector<double>> random_walks; */

  for(int n=0;n<N;++n){
    double maximum=S0;
    double prev=S0;
    /* walk.push_back(S0); */
    for (int m=1;m<M;++m){
      prev+=norm(gen);
      /* double st = walk[m-1]+norm(gen); */
      /* walk.push_back(st); */
      if(prev>maximum) maximum = prev;
    };
    result+=maximum;
    /* walk.clear(); */
  };
  return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int M = getArg(argv,2);

  auto start = std::chrono::system_clock::now();
  double result = mc_amer(N,M,S0,E,r,T,sigma);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N
      );
  return EXIT_SUCCESS;
}
