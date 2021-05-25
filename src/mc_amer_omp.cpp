
#include <common.h>
#include <comparison.h>

/* void vecprinter(std::vector<double> vec) */
/* { */
/*   for(int i = 0;i<vec.size();++i){ */
/*     std::cout<<vec[i]<< " " ; */
/*   }; */
/*   std::cout<<std::endl; */
/* } */

/* void matprinter(std::vector<std::vector<double>> mat) */
/* { */
/*   for(int i = 0;i<mat.size();++i){ */
/*     for(int j = 0;j<mat[i].size();++j){ */
/*       std::cout<<mat[i][j] << " " ; */
/*     }; */
/*     std::cout<<std::endl; */
/*   }; */
/* } */


std::vector<std::vector<double>> transpose
(
 std::vector<std::vector<double>> y
)
{
  std::vector<std::vector<double>> transposed(y[0].size());
  for(int i=0;i<y[0].size();++i){
    transposed[i].resize(y.size());
  };
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) nowait
  for(int j=0;j<y[0].size();++j){
    for(int i=0;i<y.size();++i){
      transposed[j][i] = y[i][j];
    };
  };
  }
  return transposed;
}

std::vector<std::vector<double>> pathsfinder
(
 double S0
 ,double E
 ,double r
 ,double sigma
 ,double T
 ,int N
 ,int M
 ,int thread
)
{
  double dt = T/M;
  // matrix to store paths
  std::vector<std::vector<double>> paths(M+1);
  for(int i=0;i<M+1;++i){
    paths[i].resize(N);
  };
  // make a generator from  N(0,sqrt(T))
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(dt)};
  // generate paths
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)*(n+1)*(thread+1));
    // init new path
    paths[0][n] = S0;
    paths[0][n+N/2] = S0;
    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths[m][n] = paths[m-1][n]*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths[m][n+N/2] = paths[m-1][n+N/2]*exp((r-0.5*sigma*sigma)*dt-sigma*w);
    };
  };
  /* return transpose(paths); */
  return paths;
}

std::vector<double> mat_vec_mul
(
  std::vector<std::vector<double>> x
  ,std::vector<double> y
)
{
  std::vector<double> mat(x.size());
  for(int i=0;i<x.size();++i){
    for(int j=0;j<y.size();++j){
        mat[i]+=x[i][j]*y[j]; 
      };
  };
  return mat;
}

/* std::vector<std::vector<double>> merge(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B){ */
/*   bool is_zero=1; */
/*   for(int i=0;i<A[0].size();++i){ */
/*     if (A[0][i] !=0) {is_zero=0;break;}; */
/*   }; */
/*   std::vector<std::vector<double>> C; */
/*   std::vector<double> tmp; */
/*   if(is_zero!=1){ */
/*     for(int i=0;i<A.size();++i){ */
/*       copy(A[i].begin(),A[i].end(),back_inserter(tmp)); */
/*       copy(B[i].begin(),B[i].end(),back_inserter(tmp)); */
/*       C.push_back(tmp); */
/*       tmp.clear(); */
/*      }; */
/*   }else{ */
/*     for(int i=0;i<B.size();++i){ */
/*       copy(B[i].begin(),B[i].end(),back_inserter(tmp)); */
/*       C.push_back(tmp); */
/*       tmp.clear(); */
/*      }; */
/*   }; */
/*   return C; */
/* }; */

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
  double dt = T/M;
  double result = 0;
  // calculate paths
  std::vector<std::vector<double>> paths = pathsfinder(S0,E,r,sigma,T,N,M,threads);

  /* int N_p; */
  /* if(N%threads!=0) N_p=(N+threads-N%threads)/threads; */
  /* else N_p=N/threads; */
  /* if(N_p%2!=0) ++N_p; */

  /* std::vector<std::vector<double>> paths(M+1); */
  /* for(int i=0;i<M+1;++i) paths[i].resize(N_p); */
/* #pragma omp declare reduction (merge: std::vector<std::vector<double>>: omp_out=merge(omp_out,omp_in)) */
/* #pragma omp parallel */
  /* { */
/* #pragma omp for reduction(merge:paths) schedule(dynamic,1)  //nowait */
  /* for(int i=0;i<threads;++i){ */
  /*   paths = pathsfinder(S0,E,r,sigma,T,N_p,M,omp_get_thread_num()); */
  /* }; */
  /* } */

  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when(N,M);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(N);

  for(int n=0;n<N;++n) exercise_st[n] = payoff(paths[M][n],E,payoff_fun);

  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.

  for(int m=M-1;m>0;--m){
    std::vector<double> x(N,-1);
    std::vector<double> y(N,-1);
    double sum_x = 0; double sum_x2 = 0; double sum_x3 = 0; double sum_x4 = 0; double sum_y = 0; double sum_yx = 0; double sum_yx2 = 0;
    double x_length=0;
    // save first non-zero payoff incase there is only one path in the money.
    double fst_po;
    double fst_y;
    int fst_n;
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) private(E,r,dt) nowait reduction(+:sum_x,sum_x2,sum_x3,sum_x4,sum_y,sum_yx,sum_yx2,x_length)
    for(int n=0;n<N;++n){
      double payoff_val = payoff(paths[m][n],E,payoff_fun);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length;
        // stock price at time t_m
        double exer = paths[m][n];
        x[n] = exer;
        // discounted cashflow at time t_{m+1}
        double cont = exp(-r*dt*(exercise_when[n]-m))*payoff(paths[exercise_when[n]][n],E,payoff_fun);
        y[n] = cont;

        // calc values for xTx and xTy.
        sum_x   += exer;
        sum_x2  += exer*exer;
        sum_x3  += exer*exer*exer;
        sum_x4  += exer*exer*exer*exer;
        sum_y   += cont;
        sum_yx  += cont*exer;
        sum_yx2 += cont*exer*exer;
        if(x_length==1){
          fst_po=payoff_val;
          fst_y=cont;
          fst_n=n;
        };
      };
    };
  }
    std::vector<std::vector<double>> xTx;
    std::vector<double> xTy;
    
    // if no path was in the money, skip it, because we are not interested in it.
    // when M is big and dt is small, the step m=1 might not be in money.
    if (x_length==0){ 
      continue;
    } else if(x_length==1){
    // if only 1 paths in the money, then compare current and discounted price.
      if(fst_po>fst_y){
        exercise_when[fst_n] = m;
        exercise_st[fst_n] = fst_po;
      };
      continue;
    }else if(x_length==2){
      // if only 2 paths in the money, then do linear regression
      // compose xTx and xTy
      xTx.resize(2); xTy.resize(2);
      for(int i=0;i<2;++i) xTx[i].resize(2);
      xTx[0][0] = x_length; xTx[0][1] = sum_x ;
      xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2;
      xTy[0]    = sum_y   ; xTy[1]    = sum_yx;
    }else if(x_length>2){
     // more than 2 paths in the money
      xTx.resize(3); xTy.resize(3);
      for(int i=0;i<3;++i) xTx[i].resize(3);
      xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
      xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
      xTx[2][0] = sum_x2  ; xTx[2][1] = sum_x3; xTx[2][2] = sum_x4 ;
      xTy[0]    = sum_y   ; xTy[1]    = sum_yx; xTy[2]    = sum_yx2;
    };
    
    std::vector<double> coef = mat_vec_mul(inverse(xTx),xTy);
    for(int i=0;i<N;++i){
      if(x[i]!=-1){
        double poly=0;
        if(coef.size()>2) poly=coef[2]*pow(x[i],2);
        double EYIX = coef[0] + coef[1]*x[i] + poly;
        // exercise value at t_m
        double payoff_val = payoff(x[i],E,payoff_fun);
        if (payoff_val > EYIX) {
          exercise_when[i] = m;
          exercise_st[i] = payoff_val;
        };
      };
    };
  };

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) nowait reduction(+:result)
  for(int n=0;n<N;++n){
    if(exercise_when[n]!=0) result+=exp(-r*exercise_when[n]*dt)*exercise_st[n];
  };
  }  

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
