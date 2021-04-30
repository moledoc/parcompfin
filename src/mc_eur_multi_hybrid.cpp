
#include <common.h>
#include <comparison.h>
#include <mvn.h>

double mc_eur
(
  double S0
  ,double E
  ,double r
  ,double sigma
  ,double T
  ,int N
  ,double payoff_fun
  ,int assets
  ,double rho
  ,int size
  ,int rank
)
{
  double result;
  double result_inter=0;

  int N_p;
  if(N%size!=0) N_p=(N+size-N%size)/size;
  else N_p=N/size;

  // pre-calculate random variables from N(Mu,Sigma);
  /* Eigen::MatrixXd samples = sample(N_p,assets,sigma,rho); */
  Eigen::MatrixXd samples;
#pragma omp sections
  {
#pragma omp section
  samples = sample(N_p,assets,sigma,rho);
  }

#pragma omp parallel
  {
  // payoff of the basket option will depend on the  arithmetic average of prices at maturity T.
  double w_i = 1.0/(double)assets;
#pragma omp for schedule(dynamic,1000) reduction(+:result_inter) nowait private(rho)
  for(int n=0;n<N_p;++n){
    double result_n = 0;
    for(int asset=0;asset<assets;++asset){
      // assuming same constant volatility for each underlying asset.
      result_n += w_i*S0*exp((r-pow(sigma,2)/2)*T+sigma*samples(asset,n));
    };
    result_inter += payoff(result_n,E,payoff_fun);
  };
  }

  MPI_Reduce(&result_inter,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if (rank==0) return (exp(-r*T)*result)/(double)N;
  else return 0;
}

int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  std::string payoff_fun =  argv[1];
  double S0 =               getArgD(argv,2);
  double E =                getArgD(argv,3);
  double r =                getArgD(argv,4);
  /* double sigma =       sqrt(getArgD(argv,5)); */
  double sigma =            getArgD(argv,5);
  double T =                getArgD(argv,6);
  int N =                   getArg(argv,7);
  int assets =              getArg(argv,8);
  double rho =              getArgD(argv,9);
  int threads =             getArg(argv,10);
  omp_set_num_threads(threads);

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
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

  auto start = std::chrono::system_clock::now();
  double result = mc_eur(S0,E,r,sigma,T,N,payoff_fun_d,assets,rho,size,rank);
  auto end = std::chrono::system_clock::now();

  if(rank ==0 ){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
    reporting(
        "Hybrid"
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
        ,size*1000+threads
        ,0
        ,assets
        );
  };
  MPI_Finalize();
  return EXIT_SUCCESS;
}
