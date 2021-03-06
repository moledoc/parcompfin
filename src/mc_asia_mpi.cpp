
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
  ,int size
  ,int rank
)
{
  double result;  
  double result_inter=0;

  // divide nr of iterations between processes.
  // if N is not divisible by size, then fix N, so that it would - each process does one more iteration and when N is big, it does not matter that much.
  int N_p;
  if(N%size!=0) N_p = (N+size-N%size)/size;
  else N_p = N/size;
  double dt = (double)T/(double)M;

  // prepare generator
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  gen.seed(time(&cur_time)*(1+rank));
  std::normal_distribution<> norm{0,sqrt(dt)};

  for (int n=0;n<N_p;++n){
    double St = S0;
    double I = 0;
    // reuse variables I and St, so we could avoid making a vector.
    for (int m=0;m<M;++m){
      double dBi = norm(gen);
      I += St*(1+r*dt/2+sigma*dBi/2);
      St *= exp((r-pow(sigma,2)/2)*dt+sigma*dBi);
    };
    result_inter += payoff(I/(double)M,E,payoff_fun);
  };

  MPI_Reduce(&result_inter,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0) return (exp(-r*T)*result)/((double)N_p*size);
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
  int M =                   getArg(argv,8);

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  int ierr = MPI_Init(&argc,&argv);
  if (ierr !=0){
    std::cout << "ERROR" << std::endl;
    exit(1);
  }

  /* Init size and rank */
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /* if (N%size!=0) throw std::invalid_argument("N needs to be divisible by P"); */

  auto start = std::chrono::system_clock::now();
  double result = mc_asia(S0,E,r,sigma,T,N,M,payoff_fun_d,size,rank);
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
        ,M
        );
  };
  return EXIT_SUCCESS;
}
