
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
  double dt = T/(double)M;
  // matrix to store paths
  Eigen::MatrixXd paths(N,M+1);
  //
  // make a generator from  N(0,sqrt(T))
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(dt)};
  // generate paths
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)+(n+1)*(rank+1));
    // init new path
    paths(n,0) = S0;
    paths(n+N/2,0) = S0;
    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths(n,m) = paths(n,m-1)*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths(n+N/2,m) = paths(n+N/2,m-1)*exp((r-0.5*sigma*sigma)*dt-sigma*w);
    };
  };
  return paths.transpose();
  /* return transpose(paths); */
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
  double result_p = 0;
  double result;
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
    Eigen::VectorXd x_p=Eigen::VectorXd::Zero(N_p);
    Eigen::VectorXd y_p=Eigen::VectorXd::Zero(N_p);
    Eigen::VectorXd is_in_money_p=Eigen::VectorXd::Zero(N_p);
    int in_money_p=0;
    
    for(int n=0;n<N_p;++n){
      double tmp = paths(m,n);
      is_in_money_p(n) = payoff(tmp,E,payoff_fun);
      if(is_in_money_p(n)>0) {
        ++in_money_p;
        x_p(n) = tmp;
        y_p(n) = exp(-r*dt*(exercise_when(n)-m))*payoff(paths(exercise_when(n),n),E,payoff_fun);
      };
    };

    // find how many paths where in money in each process
    int in_money;
    MPI_Allreduce(&in_money_p,&in_money,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(in_money==0) continue;

    Eigen::VectorXd x(N_p*size);
    Eigen::VectorXd y(N_p*size);
    Eigen::VectorXd is_in_money(N_p*size);
    
    /* MPI_Allgather(x_p.data(),N_p,MPI_DOUBLE,x.data(),N_p,MPI_DOUBLE,MPI_COMM_WORLD); */
    /* MPI_Allgather(y_p.data(),N_p,MPI_DOUBLE,y.data(),N_p,MPI_DOUBLE,MPI_COMM_WORLD); */
    /* MPI_Allgather(is_in_money_p.data(),N_p,MPI_DOUBLE,is_in_money.data(),N_p,MPI_DOUBLE,MPI_COMM_WORLD); */

    MPI_Gather(x_p.data(),N_p,MPI_DOUBLE,x.data(),N_p,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(y_p.data(),N_p,MPI_DOUBLE,y.data(),N_p,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(is_in_money_p.data(),N_p,MPI_DOUBLE,is_in_money.data(),N_p,MPI_DOUBLE,0,MPI_COMM_WORLD);

    Eigen::VectorXd coef(std::min(in_money,3));

    if(rank==0){
    Eigen::MatrixXd xmat=Eigen::MatrixXd::Zero(in_money,std::min(in_money,3));
    Eigen::MatrixXd yvec=Eigen::VectorXd::Zero(in_money);
    int counter = 0;
    for(int i=0;i<N_p*size;++i){
      if(is_in_money(i)>0){
        xmat(counter,0)=1;
        if (in_money > 1) xmat(counter,1)=x(i);
        if (in_money > 2) xmat(counter,2)=pow(x(i),2);
        yvec(counter)=y(i);
        ++counter;
      };
    };

    Eigen::MatrixXd xmatT = xmat.transpose();
    coef = (xmatT*xmat).inverse()*xmatT*yvec;
    };
    MPI_Bcast(coef.data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD);


    for(int i=0;i<N_p;++i){
      if(is_in_money_p(i)>0){
        double EYIX;
        double payoff_val;
        double poly=0;
        if(in_money!=1){
          if(in_money>2){
            poly = coef(2)*pow(x_p(i),2);
          };
          EYIX = coef(0) + coef(1)*x_p(i) + poly;
          // exercise value at t_m
          payoff_val = payoff(x_p(i),E,payoff_fun);
        }else {
          EYIX = coef(0);
          // exercise value at t_m
          payoff_val = payoff(x_p(i),E,payoff_fun);
        };
        // exercise value at t_m
        if (payoff_val > EYIX) {
          exercise_when(i) = m;
          exercise_st(i) = payoff_val;
        };
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
