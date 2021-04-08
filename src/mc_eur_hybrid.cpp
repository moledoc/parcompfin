
#include <example_eur.h>

double mc_eur(
    int N,double S0,double E,double r, double T,double sigma,
    time_t cur_time, int rank=0
    ){

    double result=0;
#pragma omp parallel //private(n) shared(T,N,r,sigma,E,S0)
{
    std::random_device rd{};
    std::mt19937 gen{rd()};
    gen.seed(time(&cur_time)+rank);
    std::normal_distribution<> norm{0,sqrt(T)};
#pragma omp for nowait schedule(dynamic,1000) reduction(+:result) 
    for(int n=0;n<N;++n){
      result += payoff(S0*exp((r-sigma*sigma/2)*T+sigma*norm(gen)),E);
    };
}
    return (exp(-r*T)*result)/(double)N;
}

int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int threads =getArg(argv,2); 
  omp_set_num_threads(threads);
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

  /* /1* Normalize start time *1/ */
  /* if (rank == 0){ */
  /*   MPI_Bcast(&cur_time,1,MPI_INT,0,MPI_COMM_WORLD); */
  /* }; */

  /* Fix N if necessary */
  if (N % size > 0){
    N = N + (size - (N % size));
  };

  auto start = std::chrono::system_clock::now();
  double inter_result = mc_eur(N/size,S0,E,r,T,sigma,cur_time,rank);
  double result;
  MPI_Reduce(&inter_result,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  result = result/size;
  auto end = std::chrono::system_clock::now();
  if (rank==0){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
    reporting(
        "Hybrid",
        elapsed_seconds_overall.count(),
        elapsed_seconds.count(),
        result,
        analytical,
        N,
        size*1000+threads
        );
  };
  MPI_Finalize();
  return EXIT_SUCCESS;
}