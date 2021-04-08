
#include <example_amer.h>

void vecprinter(std::vector<double> vec)
{
  for(int i = 0;i<vec.size();++i){
    std::cout<<vec[i]<< " " ;
  };
  std::cout<<std::endl;
}

void matprinter(std::vector<std::vector<double>> mat)
{
  for(int i = 0;i<mat.size();++i){
    for(int j = 0;j<mat[i].size();++j){
      std::cout<<mat[i][j] << " " ;
    };
    std::cout<<std::endl;
  };
}


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
 int N,int M,double S0
)
{
  double dt = T/M;
  // matrix to store paths
  std::vector<std::vector<double>> paths(N);
  for(int i=0;i<N;++i){
    paths[i].resize(M+1);
  };
  // make a generator from  N(0,sqrt(T))
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(dt)};
    
/* #pragma omp parallel for */
  // generate paths
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)*(n+1));

    // init new path
    paths[n][0] = S0;
    paths[n+N/2][0] = S0;

    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths[n][m] = paths[n][m-1]*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths[n+N/2][m] = paths[n+N/2][m-1]*exp((r-0.5*sigma*sigma)*dt-sigma*w);
    };
  };

  return transpose(paths);

}

std::vector<double> mat_vec_mul
(
  std::vector<std::vector<double>> x,
  std::vector<double> y
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

// in my case it is always 3x3 matrix
std::vector<std::vector<double>> inverse
(
 std::vector<std::vector<double>> x
)
{
  std::vector<std::vector<double>> inversed(3);
  for(int i=0;i<3;++i){
    inversed[i].resize(3);
  };
  double determinant=0;
  //finding determinant of the matrix
  for(int i=0; i<3;++i)
    determinant += (x[0][i] * (x[1][(i+1)%3] * x[2][(i+2)%3] - x[1][(i+2)%3] * x[2][(i+1)%3]));
  //Condition to check if the derterminat is zero or not if zero than inverse dont exists
  if(determinant<=0){
    throw std::invalid_argument("Detereminant is not > 0");
  };
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      inversed[j][i] = ((x[(j+1)%3][(i+1)%3] * x[(j+2)%3][(i+2)%3]) - (x[(j+1)%3][(i+2)%3] * x[(j+2)%3][(i+1)%3]))/determinant;
    };
   };
  return inversed;
}


double mc_amer
 (
 int N,int M,double S0,double E,double r, double T,double sigma
 )
{
  double dt = T/M;
  double result = 0;
std::vector<std::vector<double>> paths;
#pragma omp parallel
{
#pragma omp sections nowait
  {
  // calculate paths
#pragma omp section
  paths = pathsfinder(N,M,S0);
  }
}
  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when(N,M);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(N);

#pragma omp parallel for
  for(int n=0;n<N;++n) exercise_st[n] = payoff(paths[M][n],E);

  std::vector<std::vector<double>> xTx(3);
  for(int i=0;i<3;++i) xTx[i].resize(3);
  /* std::vector<std::vector<double>> xTy(1); */
  /* xTy[0].resize(3); */
  std::vector<double> xTy(3);

  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.

  for(int m=M-1;m>0;--m){
    std::vector<double> x(N,-1);
    std::vector<double> y(N,-1);
    double sum_x = 0; double sum_x2 = 0; double sum_x3 = 0; double sum_x4 = 0; double sum_y = 0; double sum_yx = 0; double sum_yx2 = 0;
    double x_length=0;
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) private(E,r,dt) nowait reduction(+:sum_x,sum_x2,sum_x3,sum_x4,sum_y,sum_yx,sum_yx2,x_length)
    for(int n=0;n<N;++n){
      double payoff_val = payoff(paths[m][n],E);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length;
        // stock price at time t_m
        double exer = paths[m][n];
        x[n] = exer;
        // discounted cashflow at time t_{m+1}
        double cont = exp(-r*dt*(exercise_when[n]-m))*payoff(paths[exercise_when[n]][n],E);
        y[n] = cont;

        // calc values for xTx and xTy.
        sum_x   += exer;
        sum_x2  += exer*exer;
        sum_x3  += exer*exer*exer;
        sum_x4  += exer*exer*exer*exer;
        sum_y   += cont;
        sum_yx  += cont*exer;
        sum_yx2 += cont*exer*exer;
      };
    };
  }
    
    // compose xTx and xTy
    xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
    xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
    xTx[2][0] = sum_x2  ; xTx[2][1] = sum_x3; xTx[2][2] = sum_x4 ;
    xTy[0]    = sum_y   ; xTy[1]    = sum_yx; xTy[2]    = sum_yx2;

/* #pragma omp parallel */
/*     { */
    std::vector<double> coef = mat_vec_mul(inverse(xTx),xTy);
/* #pragma omp for schedule(dynamic,1000) nowait */ 
    for(int i=0;i<N;++i){
      if(x[i]!=-1){
        double EYIX = coef[0] + coef[1]*x[i] + coef[2]*pow(x[i],2);
        // exercise value at t_m
        double payoff_val = payoff(x[i],E);
        if (payoff_val > EYIX) {
          exercise_when[i] = m;
          exercise_st[i] = payoff_val;
        };
      };
    };
  };
  /* } */

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) nowait reduction(+:result)
  for(int n=0;n<N;++n){
    if(exercise_when[n]!=0) result+=exp(-r*exercise_when[n]*dt)*exercise_st[n];
  };
  }  

  return std::max(payoff(S0,E),result/(double)N);
}

int main (int argc, char *argv[]){

  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int M = getArg(argv,2);
  int threads = getArg(argv,3);
  omp_set_num_threads(threads);

  

  /* matprinter(pathsfinder(10,10,100)); */
  /* std::vector<std::vector<double>> mat;//= pathsfinder(10,10,100); */
  /* int k=1; */
  /* for(int i=0;i<3;++i){ */
  /*   std::vector<double> tmp; */
  /*   for(int j=0;j<3;++j){ */
  /*     if(k==5) tmp.push_back(0); */
  /*     /1* else if(k==8) tmp.push_back(0); *1/ */
  /*     else tmp.push_back(k); */
  /*     ++k; */
  /*   }; */
  /*   mat.push_back(tmp); */
  /* }; */

  /* matprinter(mat); */
  /* std::cout << "==============================" << std::endl; */
  /* /1* matprinter(transpose(mat)); *1/ */
  /* vecprinter(transpose(mat)[0]); */
  /* std::cout << "==============================" << std::endl; */

  /* /1* std::vector<std::vector<double>> tst; *1/ */
  /* /1* std::vector<double> tst2 = mat_vec_mul(mat,transpose(mat)[0]); *1/ */
  /* /1* tst.push_back(tst2); *1/ */
  /* /1* matprinter(tst); *1/ */
  /* vecprinter(mat_vec_mul(mat,transpose(mat)[0])); */
  /* std::cout << "==============================" << std::endl; */


  auto start = std::chrono::system_clock::now();
  double result;
#pragma omp parallel
  {
#pragma omp sections
  {
#pragma omp section
  result = mc_amer(N,M,S0,E,r,T,sigma);
  }
  }
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "OMP",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N,
      M,
      threads
      );
  return EXIT_SUCCESS;
}
