
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
  for(int j=0;j<y[0].size();++j){
    for(int i=0;i<y.size();++i){
      transposed[j][i] = y[i][j];
    };
  };
  return transposed;
}

std::vector<std::vector<double>> pathsfinder
(
 int N,int M,double S0,
 int rank
)
{
  if (N%2!=0) throw std::invalid_argument("N needs to be divisible by 2");
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
    
  // generate paths
  for(int n=0;n<N/2;++n){
    /* // for each path use different seed */
    gen.seed(time(&cur_time)+(rank+1)*(n+1));

    // init new path
    /* std::vector<double> path(M+1); */
    paths[n][0] = S0;
    paths[n+N/2][0] = S0;

    // fill path
    for(int m=1;m<M+1;++m){
      double w = norm(gen);
      paths[n][m] = paths[n][m-1]*exp((r-0.5*sigma*sigma)*dt+sigma*w);
      paths[n+N/2][m] = paths[n+N/2][m-1]*exp((r-0.5*sigma*sigma)*dt-sigma*w);;
    };
    // add path to paths vector.
    /* paths[n] = path; */
  };

  return transpose(paths);
  /* return paths; */

}

/* std::vector<std::vector<double>> mat_vec_mul */
std::vector<double> mat_vec_mul
(
  std::vector<std::vector<double>> x,
  /* std::vector<std::vector<double>> y */
  std::vector<double> y
)
{
  /* std::vector<std::vector<double>> mat(x.size()); */
  std::vector<double> mat(x.size());
  /* for(int i=0;i<x.size();++i){ */
  /*   mat[i].resize(y.size()); */
  /* }; */
  /* std::vector<std::vector<double>> yT = transpose(y); */
  for(int i=0;i<x.size();++i){
    /* for(int k=0;k<y.size();++k){ */
      double tmp=0;
      for(int j=0;j<x[i].size();++j){
        /* tmp+=x[i][j]*y[k][j]; */
        tmp+=x[i][j]*y[j];
      };
      mat[i]=tmp;
    /* }; */
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
    std::vector<double> tmp;
    for(int j=0;j<3;++j){
      inversed[i][j] = ((x[(j+1)%3][(i+1)%3] * x[(j+2)%3][(i+2)%3]) - (x[(j+1)%3][(i+2)%3] * x[(j+2)%3][(i+1)%3]))/determinant;
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
  
  // calculate N/size paths in each process
  std::vector<std::vector<double>> paths = pathsfinder(N/size,M,100,rank);

  /* // allgather the paths, so that each process has all paths. */
  /* for(int i=0;i<M+1;++i){ */
  /*   std::vector<double> step;step.resize(N); */
  /*   MPI_Allgather((paths[i]).data(),N/size,MPI_DOUBLE,step.data(),N/size,MPI_DOUBLE,MPI_COMM_WORLD); */
  /*   paths[i] = step; */
  /* }; */

  // calculate length of process responsibility
  int n_start = rank*N/size;
  int n_p = ((rank+1)*N/size)-(rank*N/size);

  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when(n_p,M);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(n_p);
  for(int n=0;n<n_p;++n) exercise_st.push_back(payoff(paths[M][n_start+n],E));
  
  std::vector<std::vector<double>> xTx(3);
  for(int i=0;i<3;++i) xTx[i].resize(3);
  /* std::vector<std::vector<double>> xTy(1); */
  /* xTy[0].resize(3); */
  std::vector<double> xTy(3);


  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.
  for(int m=M-1;m>0;--m){
    std::vector<double> x_p(N/size,-1);
    std::vector<double> y_p(N/size,-1);
    double sum_x; double sum_x2; double sum_x3; double sum_x4; double sum_y; double sum_yx; double sum_yx2;
    double sum_x_p = 0; double sum_x2_p = 0; double sum_x3_p = 0; double sum_x4_p = 0; double sum_y_p = 0; double sum_yx_p = 0; double sum_yx2_p = 0;

    double x_length; double x_length_p=0;
    for(int n=0;n<n_p;++n){
      /* int n_global = n_start+n; */
      /* double payoff_val = payoff(paths[m][n_global],E); */
      double payoff_val = payoff(paths[m][n],E);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length_p;
        // stock price at time t_m
        double tmp_x = paths[m][n_start+n];
        /* x[n] = paths[m][n]; */
        x_p[n] = tmp_x;
        // discounted cashflow at time t_{m+1}
        /* double tmp_y = exp(-r*dt)*payoff(paths[exercise_when[n]][n_global],E); */
        double tmp_y = exp(-r*dt)*payoff(paths[exercise_when[n]][n],E);
        /* y[n] = exp(-r*dt)*payoff(paths[m+1][n],E); */
        y_p[n] = tmp_y;

        // calc values for xTx and xTy.
        sum_x_p   += tmp_x;
        sum_x2_p  += tmp_x*tmp_x;
        sum_x3_p  += tmp_x*tmp_x*tmp_x;
        sum_x4_p  += tmp_x*tmp_x*tmp_x*tmp_x;
        sum_y_p   += tmp_y;
        sum_yx_p  += tmp_y*tmp_x;
        sum_yx2_p += tmp_y*tmp_x*tmp_x;
      };
    };

    std::vector<double> coef(3);
    /* std::vector<double> x;x.resize(N); */
    /* std::vector<double> y;y.resize(N); */
    /* MPI_Allgather(x.data(),N/size,MPI_DOUBLE,x_p.data(),N/size,MPI_DOUBLE,MPI_COMM_WORLD); */
    /* MPI_Allgather(y.data(),N/size,MPI_DOUBLE,y_p.data(),N/size,MPI_DOUBLE,MPI_COMM_WORLD); */

    // TODO: sync x, y, matrix sums,  x_length

    /* MPI_Allreduce(&sum_x_p,  &sum_x,  1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_x2_p, &sum_x2, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_x3_p, &sum_x3, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_x4_p, &sum_x4, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_y_p,  &sum_y,  1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_yx_p, &sum_yx, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
    /* MPI_Allreduce(&sum_yx2_p,&sum_yx2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */

    /* MPI_Allreduce(&x_length_p,&x_length,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */

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

      // calc coefficients
      /* std::vector<double> coef = mat_vec_mul(inverse(xTx),xTy); */
      coef = mat_vec_mul(inverse(xTx),xTy);
    };
    MPI_Bcast(coef.data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for(int n=0;n<n_p;++n){
      if(x_p[n]!=-1){
        /* int n_global = n_start + n; */
        double EYIX = coef[0] + coef[1]*x_p[n] + coef[2]*pow(x_p[n],2);
        // exercise value at t_m
        /* double payoff_val = payoff(paths[m][n_global],E); */
        double payoff_val = payoff(paths[m][n],E);
        if (payoff_val > EYIX) {
          exercise_when[n] = m;
          exercise_st[n] = payoff_val;
        };
      };
    };
  };
  
/* TODO: parallelize and then reduce */
  for(int n=0;n<n_p;++n){
    if(exercise_st[n]!=0) result_p+=exp(-r*exercise_when[n])*exercise_st[n];
  };
    /* std::cout <<result_p<< std::endl; */

/* TODO: vb teha result ainult p=0. */
  MPI_Reduce(&result_p,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(rank==0){
    return result/(double)N;
  } else return 0;
  /* MPI_Allreduce(&result_p,&result,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); */
  /* return result/(double)N; */
}

int main (int argc, char *argv[]){
  /* std::vector<std::vector<double>> mat= pathsfinder(1000,100,100); */
  /* std::vector<std::vector<double>> mat;//= pathsfinder(10,10,100); */
  /* std::vector<double> tmp1(3); */
  /* std::vector<double> tmp2(3); */
  /* std::vector<double> tmp3(3); */
  /* mat.push_back(tmp1); */
  /* mat.push_back(tmp2); */
  /* mat.push_back(tmp3); */
  /* mat[0][0] = 5; */
  /* mat[0][1] = 7; */
  /* mat[0][2] = 9; */
  /* mat[1][0] = 4; */
  /* mat[1][1] = 3; */
  /* mat[1][2] = 8; */
  /* mat[2][0] = 7; */
  /* mat[2][1] = 5; */
  /* mat[2][2] = 6; */
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
  /* std::vector<double> tmp2; */
  /* tmp.push_back(0); */
  /* tmp.push_back(2); */
  /* tmp2.push_back(3); */
  /* tmp2.push_back(4); */
  /* mat.push_back(tmp);mat.push_back(tmp2); */

  /* matprinter(mat); */
  /* std::cout << "==============================" << std::endl; */

  /* std::vector<std::vector<double>> mat_t = transpose(mat); */
  /* matprinter(mat_t); */
  /* std::cout << "==============================" << std::endl; */

  /* std::vector<std::vector<double>> mat_x = matmul(mat,mat_t); */
  /* matprinter(mat_x); */
  /* std::cout << "==============================" << std::endl; */

  /* std::vector<std::vector<double>> mat_inv = inverse(mat); */
  /* matprinter(mat_inv); */
  /* std::cout << "==============================" << std::endl; */


  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
  int M = getArg(argv,2);
  /* if((M+1)%2!=0) throw std::invalid_argument("M+1 needs to be divisible by 2"); */

  time_t cur_time;
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
  MPI_Status status;


  /* auto start = std::chrono::system_clock::now(); */
  /* // calculate N/size paths in each process */
  /* std::vector<std::vector<double>> mat= pathsfinder(N/size,M,100,rank); */

  /* // allgather the paths, so that each process has all paths. */
  /* for(int i=0;i<M+1;++i){ */
  /*   std::vector<double> step;step.resize(N); */
  /*   MPI_Allgather((mat[i]).data(),N/size,MPI_DOUBLE,step.data(),N/size,MPI_DOUBLE,MPI_COMM_WORLD); */
  /*   mat[i] = step; */
  /* }; */


  /* auto end = std::chrono::system_clock::now(); */
  /* std::chrono::duration<double> elapsed_seconds = end-start; */

  /* if(rank==0)std::cout << elapsed_seconds.count() << std::endl; */

  /* if(rank==0){ */
  /*   auto start = std::chrono::system_clock::now(); */
  /*   std::vector<std::vector<double>> mat= pathsfinder(N,M,100,rank); */
  /*   auto end = std::chrono::system_clock::now(); */
  /*   std::chrono::duration<double> elapsed_seconds = end-start; */
  /*   std::cout << elapsed_seconds.count() << std::endl; */
  /* }; */



  auto start = std::chrono::system_clock::now();
  double result = mc_amer(N,M,S0,E,r,T,sigma,rank,size);
  auto end = std::chrono::system_clock::now();

  if(rank==0){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
    reporting(
        "MPI",
        elapsed_seconds_overall.count(),
        elapsed_seconds.count(),
        result,
        analytical,
        N,
        M,
        size
        );
  };
  MPI_Finalize();
  return EXIT_SUCCESS;
}
