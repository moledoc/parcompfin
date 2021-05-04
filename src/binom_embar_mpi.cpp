
#include <common.h>
#include <comparison.h>

double binom 
(
  double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,int size
  ,int rank
  ,double payoff_fun
 )
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

  double V0=0;
  int ni;
  int ni1;
  if (N%2==0){
    ni = rank * (N/2)/size;
    ni1 = (rank+1) * (N/2)/size;
    if (rank == 0) {
      double binom_mid = comb(N,N/2) + N/2*log(p) + N/2*log(q);
      V0 = exp(binom_mid) * payoff(S0*pow(u,N/2)*pow(d,N/2),E,payoff_fun);
    };
  } else {
    ni = rank * (N/2)/size;
    ni1 = (rank+1) * (N/2)/size;
  };

  ni1 = std::min(ni1,until);
  for(int i=ni;i<ni1;++i){
    double comb_val = comb(N,i);
    double binom1 = comb_val + i*log(p) + (N-i)*log(q);
    double binom2 = comb_val + (N-i)*log(p) + i*log(q);
    V0 +=  exp(binom1) * payoff(S0*pow(u,i)*pow(d,N-i),E,payoff_fun);
    V0 +=  exp(binom2) * payoff(S0*pow(u,N-i)*pow(d,i),E,payoff_fun);
  };

  double result;
  MPI_Reduce(&V0,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if(rank==0) return exp(-r*T)*result;
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
  if(payoff_fun=="put") payoff_fun_d=-1;
  if(payoff_fun=="call") payoff_fun_d=1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

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
  double result = binom(S0,E,r,sigma,T,N,size,rank,payoff_fun_d);
  auto end = std::chrono::system_clock::now();
  
  // close processes
  MPI_Finalize();
  auto end_overall = std::chrono::system_clock::now();

  if (rank==0){
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
