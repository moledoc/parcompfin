
#include <example_amer.h>

double mc_eur
(
 int N,int M,double S0,double E,double r, double T,double sigma, 
 /* time_t cur_time, int rank=0 */
 int rank=0
)
{
  double result=0;

  for(int n=0;n<N;++n){

    time_t cur_time;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    gen.seed(time(&cur_time)+rank);
    std::normal_distribution<> norm{0,sqrt(T)};

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

  /* time_t cur_time; */
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
  double N_fixed = N;
  if (N % size > 0){
    N_fixed = N + (size - (N % size));
  };

  auto start = std::chrono::system_clock::now();
  /* double inter_result = mc_eur(N/size,M,S0,E,r,T,sigma,cur_time,rank); */
  double inter_result = mc_eur(N_fixed/size,M,S0,E,r,T,sigma,rank);
  double result;
  MPI_Reduce(&inter_result,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  result = result/size;
  auto end = std::chrono::system_clock::now();
  if (rank==0){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
    reporting(
        "MPI",
        elapsed_seconds_overall.count(),
        elapsed_seconds.count(),
        result,
        analytical,
        N,
        M,
        size
        );
  };
  MPI_Finalize();
  return EXIT_SUCCESS;
}
