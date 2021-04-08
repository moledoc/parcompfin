
#include <example_amer.h>

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
 ,int rank
)
{
  if (N%2!=0) throw std::invalid_argument("N needs to be divisible by 2 for finding paths");
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

  /* gen.seed(time(&cur_time)); */
    
  // generate paths
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)+(rank+1)*(n+1));

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
        mat[i]+=x[i][j]*y[j]; // yT*xT
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
 ,int rank,int size
 )
{
  double dt = T/M;
  double result_p = 0;
  double result;
  // calculate paths
  std::vector<std::vector<double>> paths;
#pragma omp parallel
  {
#pragma omp sections
    {
#pragma omp section
  paths = pathsfinder(N/size,M,S0,rank);
    }
  }

  // calculate length of process responsibility
  /* int n_start = rank*N/size; */
  /* int n_p = ((rank+1)*N/size)-(rank*N/size); */
  int n_p = N/size;

  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when(n_p,M);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(n_p);
  for(int n=0;n<n_p;++n) exercise_st[n] = payoff(paths[M][n],E);
  
  std::vector<std::vector<double>> xTx(3);
  for(int i=0;i<3;++i) xTx[i].resize(3);
  std::vector<double> xTy(3);


  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.
  for(int m=M-1;m>0;--m){
    std::vector<double> x(n_p,-1);
    std::vector<double> y(n_p,-1);
    double sum_x; double sum_x2; double sum_x3; double sum_x4; double sum_y; double sum_yx; double sum_yx2;
    double sum_x_p = 0; double sum_x2_p = 0; double sum_x3_p = 0; double sum_x4_p = 0; double sum_y_p = 0; double sum_yx_p = 0; double sum_yx2_p = 0;

    double x_length; double x_length_p=0;

/* #pragma omp parallel */
/*   { */
/* #pragma omp for schedule(dynamic,1000) private(E,r,dt) nowait reduction(+:sum_x_p,sum_x2_p,sum_x3_p,sum_x4_p,sum_y_p,sum_yx_p,sum_yx2_p,x_length_p) */
    for(int n=0;n<n_p;++n){
      double payoff_val = payoff(paths[m][n],E);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length_p;
        // stock price at time t_m
        double exer = paths[m][n];
        x[n] = exer;
        double cont = exp(-r*dt*(exercise_when[n]-m))*payoff(paths[exercise_when[n]][n],E);
        y[n] = cont;

        // calc values for xTx and xTy.
        sum_x_p   += exer;
        sum_x2_p  += exer*exer;
        sum_x3_p  += exer*exer*exer;
        sum_x4_p  += exer*exer*exer*exer;
        sum_y_p   += cont;
        sum_yx_p  += cont*exer;
        sum_yx2_p += cont*exer*exer;
      };
    };
  /* } */

    std::vector<double> coef(3);

    MPI_Reduce(&sum_x_p,  &sum_x,  1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x2_p, &sum_x2, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x3_p, &sum_x3, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x4_p, &sum_x4, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_y_p,  &sum_y,  1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx_p, &sum_yx, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx2_p,&sum_yx2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    MPI_Reduce(&x_length_p,&x_length,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    if(rank==0){
      // compose xTx and xTy
      xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
      xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
      xTx[2][0] = sum_x2  ; xTx[2][1] = sum_x3; xTx[2][2] = sum_x4 ;
      xTy[0]    = sum_y   ; xTy[1]    = sum_yx; xTy[2]    = sum_yx2;

      coef = mat_vec_mul(inverse(xTx),xTy);
    };


    MPI_Bcast(coef.data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
    for(int i=0;i<n_p;++i){
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

#pragma omp parallel
  {
#pragma omp for schedule(dynamic,100) nowait reduction(+:result_p)
  for(int n=0;n<n_p;++n){
    if(exercise_st[n]!=0) result_p+=exp(-r*exercise_when[n]*dt)*exercise_st[n];
    /* if(exercise_when[n]!=0) result+=exp(-r*dt)*exercise_st[n]; */
  };
  }

  MPI_Reduce(&result_p,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0) return std::max(payoff(S0,E),result/(double)N);
  else return 0;
}

int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int M = getArg(argv,2);
  int threads = getArg(argv,3);
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

  auto start = std::chrono::system_clock::now();
  double result = mc_amer(N,M,S0,E,r,T,sigma,rank,size);
  auto end = std::chrono::system_clock::now();

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
        M,
        size*1000+threads
        );
  };
  MPI_Finalize();
  return EXIT_SUCCESS;
}
