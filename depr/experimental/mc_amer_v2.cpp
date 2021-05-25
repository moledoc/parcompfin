
#include <common.h>
#include <comparison.h>
#include <mvn.h>

/* void vecprinter(Eigen::VectorXd vec) */
/* { */
/*   for(int i = 0;i<vec.size();++i){ */
/*     std::cout<<vec(i)<< " " ; */
/*   }; */
/*   std::cout<<std::endl; */
/* } */


/* void matprinter(Eigen::MatrixXd mat) */
/* { */
/*   for(int i = 0;i<mat.size();++i){ */
/*     for(int j = 0;j<mat(i).size();++j){ */
/*       std::cout<<mat(i,j) << " " ; */
/*     }; */
/*     std::cout<<std::endl; */
/*   }; */
/* } */


/* Eigen::MatrixXd transpose */
/* ( */
/*  Eigen::MatrixXd y */
/* ) */
/* { */
/*   Eigen::MatrixXd transposed(y.cols(),y.rows()); */
/*   for(int j=0;j<y.cols();++j){ */
/*     for(int i=0;i<y.rows();++i){ */
/*       transposed(j,i) = y(i,j); */
/*     }; */
/*   }; */
/*   return transposed; */
/* } */

Eigen::MatrixXd pathsfinder
(
 double S0
 ,double E
 ,double r
 ,double sigma
 ,double T
 ,int N
 ,int M
)
{
  double dt = T/M;
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
    gen.seed(time(&cur_time)+n+1);
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

Eigen::VectorXd mat_vec_mul
(
  Eigen::MatrixXd x
  ,Eigen::VectorXd y
)
{
  Eigen::VectorXd mat(x.size());
  for(int i=0;i<x.size();++i){
    for(int j=0;j<y.size();++j){
        mat(i)+=x(i,j)*y(j);
      };
  };
  return mat;
}

/* // in my case it is always 3x3 matrix */
/* Eigen::MatrixXd inverse */
/* ( */
/*  Eigen::MatrixXd x */
/* ) */
/* { */
/*   Eigen::MatrixXd inversed(3,3); */

/*   double determinant=0; */
/*   //finding determinant of the matrix */
/*   for(int i=0; i<3;++i) */
/*     determinant += (x(0,i) * (x(1,(i+1)%3) * x(2,(i+2)%3) - x(1,(i+2)%3) * x(2,(i+1)%3))); */
/*   //Condition to check if the derterminat is zero or not if zero than inverse dont exists */
/*   if(determinant<=0){ */
/*     throw std::invalid_argument("Detereminant is not > 0"); */
/*   }; */
/*   for(int i=0;i<3;++i){ */
/*     for(int j=0;j<3;++j){ */
/*       inversed(j,i) = ((x((j+1)%3,(i+1)%3) * x((j+2)%3,(i+2)%3)) - (x((j+1)%3,(i+2)%3) * x((j+2)%3,(i+1)%3)))/determinant; */
/*     }; */
/*    }; */
/*   return inversed; */
/* } */


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
 )
{
  double dt = T/M;
  double result = 0;
  // calculate paths
  Eigen::MatrixXd paths = pathsfinder(S0,E,r,sigma,T,N,M);
  
  // store each paths timestep value when option is exercised
  Eigen::VectorXd exercise_when(N);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  Eigen::VectorXd exercise_st(N);

  for(int n=0;n<N;++n){ 
    exercise_when(n) = M;
    exercise_st(n) = payoff(paths(M,n),E,payoff_fun);
  };

  for(int m=M-1;m>M-2;--m){
    Eigen::ArrayXXf info(3,N);
    
    for(int n=0;n<N;++n){
      double tmp = paths(m,n);
      info(0,n) = payoff(tmp,E,payoff_fun);
      info(1,n) = tmp;
      info(2,n) = exercise_when(n);
    }

    int in_money=(info.row(0)>0).count();
    if(in_money==0) continue;

    Eigen::MatrixXd x(in_money,3);
    Eigen::VectorXd y(in_money);

    int counter=0;
    for(int n=0;n<N;++n){
      if(info(0,n)>0){
        x(counter,0) = 1;
        x(counter,1) = info(1,n);
        x(counter,2) = pow(info(1,n),2);
        y(counter) = exp(-r*dt*(info(2,n)-m))*payoff(paths(info(2,n),n),E,payoff_fun);
        ++counter;
      };
    };

    Eigen::MatrixXd xT = x.transpose();
    /* Eigen::MatrixXd xTx = xT*x; */
    /* Eigen::MatrixXd xTy = xT*y; */
    /* Eigen::MatrixXd xTx_inv = xTx.inverse(); */
    /* Eigen::VectorXd coef = xTx_inv*xTy; */

    Eigen::VectorXd coef = (xT*x).inverse()*xT*y;

    counter=0;
    for(int n=0;n<N;++n){
      if(info(0,n)>0){
        double tmp_x = x(counter,1);
        double EYIX = coef(0) + coef(1)*x(counter,1) + coef(2)*pow(x(counter,1),2);
        // exercise value at t_m
        double payoff_val = payoff(x(counter,1),E,payoff_fun);
        if (payoff_val > EYIX) {
          exercise_when(n) = m;
          exercise_st(n) = payoff_val;
        };
        ++counter;
      };
    };
  };
  

  for(int n=0;n<N;++n){
    if(exercise_st(n)!=0) result+=exp(-r*exercise_when(n)*dt)*exercise_st(n);
  };

  return std::max(payoff(S0,E,payoff_fun),result/(double)N);
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

  /* std::cout << pathsfinder(S0,E,r,sigma,T,N,M) << std::endl; */ 
  
  auto start = std::chrono::system_clock::now();
  double result = mc_amer(S0,E,r,sigma,T,N,M,payoff_fun_d);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial"
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
      ,0
      ,M
      );
  return EXIT_SUCCESS;
}
