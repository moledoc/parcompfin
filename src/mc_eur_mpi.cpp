
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
  ,double payoff_fun
  ,time_t cur_time
  ,int size
  ,int rank
)
{

  /* Fix N if necessary */
  double N_p;
  if (N%size!=0) N_p = (N+size-N%size)/size;
  else N_p=N/size;

  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time)*(rank+1));
  std::normal_distribution<> norm{0,sqrt(T)};

  double result;
  double result_inter=0;

  for(int n=0;n<N_p;++n){
    result_inter += payoff(S0*exp((r-pow(sigma,2)/2)*T+sigma*norm(gen)),E,payoff_fun);
  };

  MPI_Reduce(&result_inter,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0) return (exp(-r*T)*result)/(double)N_p;
  else return 0;
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

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  time_t cur_time;
  /* Init MPI */
  int ierr = MPI_Init(&argc,&argv);
  if (ierr !=0){
    std::cout << "ERROR" << std::endl;
    exit(1);
  }

  /* Init size and rank */
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  auto start = std::chrono::system_clock::now();
  double result = mc_eur(S0,E,r,sigma,T,N,payoff_fun_d,cur_time,size,rank);
  auto end = std::chrono::system_clock::now();

  // close processes
  MPI_Finalize();
  auto end_overall = std::chrono::system_clock::now();

  if(rank==0){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end_overall-start_overall;
    reporting(
        "MPI"
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
        ,size
        );
  };
  return EXIT_SUCCESS;
}
