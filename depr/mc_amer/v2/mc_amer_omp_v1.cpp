
#include <common.h>
#include <comparison.h>
#include <mvn.h>

Eigen::MatrixXd merge(Eigen::MatrixXd A,Eigen::MatrixXd B){
  if (A.isZero(0)){
    return B;
  }
  Eigen::MatrixXd C(A.rows(),A.cols()+B.cols());
  C << A,B;
  return C;
};

Eigen::MatrixXd pathsfinder
(
 double S0
 ,double E
 ,double r
 ,double sigma
 ,double T
 ,int N
 ,int M
 ,int threads
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
    gen.seed(time(&cur_time)+(n+1)*(threads+1));
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
  ,int threads
 )
{
  double dt = T/(double)M;
  double result = 0;

  int N_p;
  if(N%threads!=0) N_p=(N+threads-N%threads)/threads;
  else N_p=N/threads;
  if(N_p%2!=0) ++N_p;
  int N_fixed = N_p*threads;
  // calculate paths
  Eigen::MatrixXd paths(M+1,N_p);

#pragma omp declare reduction (merge: Eigen::MatrixXd: omp_out=merge(omp_out,omp_in))//\
     initializer(omp_priv=Eigen::MatrixXd::Zero(omp_orig.rows(),omp_orig.cols()))
     /* initializer(omp_priv(omp_orig)) */

#pragma omp parallel
  {
#pragma omp for reduction(merge:paths) schedule(dynamic,1)  //nowait
  for(int i=0;i<threads;++i){
    paths = pathsfinder(S0,E,r,sigma,T,N_p,M,omp_get_thread_num());
  };
  }

  // store each paths timestep value when option is exercised
  Eigen::VectorXd exercise_when(N_fixed);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  Eigen::VectorXd exercise_st(N_fixed);

#pragma omp parallel for 
  for(int n=0;n<N_fixed;++n){ 
    exercise_when(n) = M;
    exercise_st(n) = payoff(paths(M,n),E,payoff_fun);
  };

  for(int m=M-1;m>0;--m){
    Eigen::ArrayXXf info(3,N_fixed);
    
#pragma omp parallel
    {
#pragma omp for nowait schedule(dynamic,100) //private(N)
    for(int n=0;n<N_fixed;++n){
      double tmp = paths(m,n);
      info(0,n) = payoff(tmp,E,payoff_fun);
      info(1,n) = tmp;
      info(2,n) = exercise_when(n);
    };
    }

    int in_money=(info.row(0)>0).count();
    if(in_money==0) continue;

    if(in_money==1){
      for(int n=0;n<N_fixed;++n){
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

    // handle cases where in money < 3 with:
    // * 1 - constant
    // * 2 - linear regression
    Eigen::MatrixXd x(in_money,std::min(in_money,3));
    Eigen::VectorXd y(in_money);

    int counter=0;
#pragma omp parallel 
    {
#pragma omp for schedule(dynamic,1000) nowait //private(N) //E,payoff_fun,m,r,dt,info)
    for(int n=0;n<N_fixed;++n){
      if(info(0,n)>0){
        x(counter,0) = 1;
        x(counter,1) = info(1,n);
        if (in_money > 2) x(counter,2) = pow(info(1,n),2);
        y(counter) = exp(-r*dt*(info(2,n)-m))*payoff(paths(info(2,n),n),E,payoff_fun);
        ++counter;
      };
    };
    }

    Eigen::MatrixXd xTx;
    Eigen::MatrixXd xTy;
    Eigen::MatrixXd xTx_inv;

    Eigen::MatrixXd xT = x.transpose();
/* #pragma omp sections */
/*     { */
/* #pragma omp section */
    xTx = xT*x;
/* #pragma omp section */
    xTy = xT*y;
    /* } */
    xTx_inv = xTx.inverse();

    Eigen::VectorXd coef = xTx_inv*xTy;

    counter=0;
/* #pragma omp parallel */
/*     { */
/* #pragma omp for schedule(dynamic,1000) nowait //private(E,payoff_fun,m,coef,x,exercise_when,exercise_st) */
    for(int n=0;n<N_fixed;++n){
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
    /* } */

  };

  int counter=0;
#pragma omp parallel
  {
#pragma omp for nowait reduction(+:result) schedule(dynamic,1000) //private(r,dt)
  for(int n=0;n<N_fixed;++n){
    if(exercise_st(n)!=0) {++counter;result+=exp(-r*exercise_when(n)*dt)*exercise_st(n);};
  };
  }
  std::cout << counter << std::endl;

  return std::max(payoff(S0,E,payoff_fun),result/(double)N_fixed);
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
  int threads =             getArg(argv,9);
  omp_set_num_threads(threads);

  double payoff_fun_d;
  if (payoff_fun=="call") payoff_fun_d = 1;
  if (payoff_fun=="put") payoff_fun_d = -1;
  if(payoff_fun != "call" && payoff_fun != "put") throw std::invalid_argument("Unknown payoff function");

  auto start = std::chrono::system_clock::now();
  double result = mc_amer(S0,E,r,sigma,T,N,M,payoff_fun_d,threads);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "OMP"
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
      ,threads
      ,M
      );
  return EXIT_SUCCESS;
}
