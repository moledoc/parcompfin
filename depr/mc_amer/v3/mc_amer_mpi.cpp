
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
    gen.seed(time(&cur_time)+(n+1)*(rank+1));
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
  double dt = T/(double)M;
  double result;
  double result_p = 0;
  int N_p;
  if(N%size!=0) N_p = (N+size-N%size)/size;
  else N_p = N/size;
  if(N_p%2!=0) ++N_p;
  double N_fixed=N_p*size;
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
    Eigen::VectorXd x(N_p);
    Eigen::VectorXd y(N_p);
    double sum_x = 0; double sum_x2 = 0; double sum_x3 = 0; double sum_x4 = 0; double sum_y = 0; double sum_yx = 0; double sum_yx2 = 0;

    double sum_x_p = 0; double sum_x2_p = 0; double sum_x3_p = 0; double sum_x4_p = 0; double sum_y_p = 0; double sum_yx_p = 0; double sum_yx2_p = 0;

    double x_length; double x_length_p=0;

    double fst_po;
    double fst_y;
    int fst_n;
    for(int n=0;n<N_p;++n){
      double payoff_val = payoff(paths(m,n),E,payoff_fun);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length_p;
        // stock price at time t_m
        double exer = paths(m,n);
        x(n) = exer;
        // discounted cashflow at time t_{m+1}
        double cont = exp(-r*dt*(exercise_when(n)-m))*payoff(paths(exercise_when(n),n),E,payoff_fun);
        y(n) = cont;
        // calc values for xTx and xTy.
        sum_x_p   += exer;
        sum_x2_p  += exer*exer;
        sum_x3_p  += exer*exer*exer;
        sum_x4_p  += exer*exer*exer*exer;
        sum_y_p   += cont;
        sum_yx_p  += cont*exer;
        sum_yx2_p += cont*exer*exer;
        if(x_length_p==1){
          fst_po=payoff_val;
          fst_y=cont;
          fst_n=n;
        };
      };
    };


    MPI_Reduce(&sum_x_p  ,&sum_x  ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x2_p ,&sum_x2 ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x3_p ,&sum_x3 ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x4_p ,&sum_x4 ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_y_p  ,&sum_y  ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx_p ,&sum_yx ,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx2_p,&sum_yx2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Allreduce(&x_length_p,&x_length,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    Eigen::MatrixXd xTx;
    Eigen::VectorXd xTy;
    // if no path was in the money, skip it, because we are not interested in it.
    // when M is big and dt is small, the step m=1 might not be in money.
    if(x_length==0){
      continue;
    } else if(x_length==1){
    // if only 1 paths in the money, then compare current and discounted price.
      if(fst_po>fst_y){
        exercise_when(fst_n) = m;
        exercise_st(fst_n) = fst_po;
      };
      continue;
    } else if(x_length==2){
    // if only 2 paths in the money, then do linear regression
      xTx.resize(2,2); xTy.resize(2);
      xTx << x_length,sum_x,sum_x2,
             sum_x,sum_x2,sum_x3;
      xTy << sum_y,sum_yx;
   // more than 2 paths in the money
    } else if(x_length>2){
      xTx.resize(3,3); xTy.resize(3);
      xTx << x_length,sum_x,sum_x2,
             sum_x,sum_x2,sum_x3,
             sum_x2,sum_x3,sum_x4;
      xTy << sum_y,sum_yx,sum_yx2;
    }; 
     
    Eigen::VectorXd coef = xTx.inverse()*xTy;
    for(int n=0;n<N_p;++n){
      if(x(n)>0){
        double poly=0;
        if (x_length>2) poly=coef[2]*pow(x(n),2);
        double EYIX = coef(0) + coef(1)*x(n) + poly;
        // exercise value at t_m
        double payoff_val = payoff(x(n),E,payoff_fun);
        if (payoff_val > EYIX) {
          exercise_when(n) = m;
          exercise_st(n) = payoff_val;
        };
      };
    };
  };

  for(int n=0;n<N_p;++n){
    if(exercise_st(n)!=0) result_p+=exp(-r*exercise_when(n)*dt)*exercise_st(n);
  };

  MPI_Reduce(&result_p,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0) return std::max(payoff(S0,E,payoff_fun),result/(double)N_fixed);
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
