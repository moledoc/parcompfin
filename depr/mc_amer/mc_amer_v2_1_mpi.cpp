
#include <common.h>
#include <comparison.h>
#include <mvn.h>

Eigen::MatrixXd pathsfinder
(
 double S0
 ,double E
 ,double r
 ,double sigma
 ,double T
 ,int N
 ,int M
 ,int rank
)
{
  double dt = T/M;
  // matrix to store paths
  Eigen::MatrixXd paths(M+1,N);
  //
  // make a generator from  N(0,sqrt(T))
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(dt)};
  // generate paths
/* #pragma omp parallel for */
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)+(n+1));
    // init new path
    paths(0,n) = S0;
    paths(0,n+N/2) = S0;
    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths(m,n) = paths(m-1,n)*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths(m,n+N/2) = paths(m-1,n+N/2)*exp((r-0.5*sigma*sigma)*dt-sigma*w);
    };
  };
  /* return paths.transpose(); */
  return paths;
}

double mc_amer
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
  double dt = T/M;
  double result;
  double result_p = 0;
  // calculate the number of paths process is responsible for
  int N_p;
  if(N%(2*size)!=0) N_p = N/size+1;
  else N_p = N/size;
  // calculate paths
  Eigen::MatrixXd paths = pathsfinder(S0,E,r,sigma,T,N_p,M,rank);
  
  // store each paths timestep value when option is exercised
  Eigen::VectorXd exercise_when(N_p);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  Eigen::VectorXd exercise_st(N_p);

  for(int n=0;n<N_p;++n){ 
    exercise_when(n) = M;
    exercise_st(n) = payoff(paths(M,n),E,payoff_fun);
  };

  for(int m=M-1;m>0;--m){
    Eigen::ArrayXXf info(3,N);
    
    for(int n=0;n<N_p;++n){
      double tmp = paths(m,n);
      info(0,n) = payoff(tmp,E,payoff_fun);
      info(1,n) = tmp;
      info(2,n) = exercise_when(n);
    };

    int in_money=(info.row(0)>0).count();
    if(in_money==0) continue;

    if(in_money==1){
      for(int n=0;n<N_p;++n){
        if(info(0,n)>0){
          double payoff_val = info(0,n);
          double discounted =  exp(-r*dt*(info(2,n)-m))*payoff(paths(info(2,n),n),E,payoff_fun);
          if (payoff_val > discounted) {
            exercise_when(n) = m;
            exercise_st(n) = payoff_val;
          };
          break;
        };
      };
      continue;
    };

    Eigen::MatrixXd x(in_money,std::min(in_money,3));
    Eigen::VectorXd y(in_money);

    int counter=0;
    for(int n=0;n<N_p;++n){
      if(info(0,n)>0){
        x(counter,0) = 1;
        x(counter,1) = info(1,n);
        if (in_money > 2) x(counter,2) = pow(info(1,n),2);
        y(counter) = exp(-r*dt*(info(2,n)-m))*payoff(paths(info(2,n),n),E,payoff_fun);
        ++counter;
      };
    };

    Eigen::MatrixXd xT = x.transpose();
    Eigen::VectorXd coef = (xT*x).inverse()*xT*y;

    counter=0;
    for(int n=0;n<N_p;++n){
      if(info(0,n)>0){
        double EYIX;
        double payoff_val;
        double poly=0;
        if(in_money>2){
          poly = coef(2)*pow(x(counter,1),2);
        };
        EYIX = coef(0) + coef(1)*x(counter,1) + poly;
        // exercise value at t_m
        payoff_val = payoff(x(counter,1),E,payoff_fun);

        if (payoff_val > EYIX) {
          exercise_when(n) = m;
          exercise_st(n) = payoff_val;
        };
        ++counter;
      };
    };
  };
  

  for(int n=0;n<N_p;++n){
    if(exercise_st(n)!=0) result_p+=exp(-r*exercise_when(n)*dt)*exercise_st(n);
  };

  MPI_Reduce(&result_p,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0) return std::max(payoff(S0,E,payoff_fun),result/(double)N);
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
  double result = mc_amer(S0,E,r,sigma,T,N,M,payoff_fun_d,size,rank);
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
        ,M
        );
  };
  return EXIT_SUCCESS;
}
