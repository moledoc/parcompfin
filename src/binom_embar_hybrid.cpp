
#include <example_eur.h>

double binom 
(
 double r, double sigma, double S0,
 double T, int N, double E,
 int size, int rank
 )
{
  double V0=0;

  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta-1); */
  /* double d = beta - sqrt(beta*beta-1); */
  /* double p = (exp(r*dt)-d)/(u-d); */
  /* double q = 1-p; */

#pragma omp parallel
  {
  double dt = T/(double)N;
  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;


  int until;
  if (N%2!=0) until = (N+1)/2;
  else until = N/2;

  int ni;
  int ni1;
  if (N%2==0){
    ni = rank * (N/2)/size;
    ni1 = (rank+1) * (N/2)/size;
    if (rank == 0) V0=comb(N,N/2)*pow(p,N/2)*pow(q,N/2)*payoff(S0*pow(u,N/2)*pow(d,N/2),E);
  } else {
    ni = rank * (N/2)/size;
    ni1 = (rank+1) * (N/2)/size;
  };

  ni1 = std::min(ni1,until);
#pragma omp for schedule(dynamic,100) nowait //reduction(+:V0)
  for(int i=ni;i<ni1;++i){
    double tmp = comb(N,i);
    V0+=tmp*pow(p,i)*pow(q,N-i)*payoff(S0*pow(u,i)*pow(d,N-i),E);
    V0+=tmp*pow(p,N-i)*pow(q,i)*payoff(S0*pow(u,N-i)*pow(d,i),E);
  };
  }
    
  double result;
  MPI_Reduce(&V0,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if(rank==0) return exp(-r*T)*result;
  else return 0;
}


int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int threads = getArg(argv,2);
  omp_set_num_threads(threads);

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

  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */
  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E,size,rank);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
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
