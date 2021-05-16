
#include <common.h>
#include <comparison.h>

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
 ,int rank
 ,int thread
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

  // generate paths
  for(int n=0;n<N/2;++n){
    // for each path use different seed
    gen.seed(time(&cur_time)+(rank+1)*(n+1)*(thread+1));

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

std::vector<std::vector<double>> merge(std::vector<std::vector<double>> A,std::vector<std::vector<double>> B){
  bool is_zero=1;
  for(int i=0;i<A[0].size();++i){
    if (A[0][i] !=0) {is_zero=0;break;};
  };
  std::vector<std::vector<double>> C;
  std::vector<double> tmp;
  if(is_zero!=1){
    for(int i=0;i<A.size();++i){
      copy(A[i].begin(),A[i].end(),back_inserter(tmp));
      copy(B[i].begin(),B[i].end(),back_inserter(tmp));
      C.push_back(tmp);
      tmp.clear();
     };
  }else{
    for(int i=0;i<B.size();++i){
      copy(B[i].begin(),B[i].end(),back_inserter(tmp));
      C.push_back(tmp);
      tmp.clear();
     };
  };
  return C;
};

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
  ,int threads
 )
{
  double dt = T/M;
  double result_p = 0;
  double result;
  // calculate the number of paths process is responsible for
  int N_p;
  
  // calculate paths
  // make N_p for paths: N is divided by processe with are then divided by threads
  
  int parallel=(threads*size);
  if(N%parallel!=0) N_p=(N+parallel-N%parallel)/parallel;
  else N_p=N/parallel;
  if(N_p%2!=0) ++N_p ;

  std::vector<std::vector<double>> paths(M+1);
  for(int i=0;i<M+1;++i) paths[i].resize(N_p);
#pragma omp declare reduction (merge: std::vector<std::vector<double>>: omp_out=merge(omp_out,omp_in))
#pragma omp parallel
  {
#pragma omp for reduction(merge:paths) schedule(dynamic,1)  //nowait
  for(int i=0;i<threads;++i){
    paths = pathsfinder(S0,E,r,sigma,T,N_p,M,rank,omp_get_thread_num());
  };
  }

  // make N_p for rest of program, where N is divided by processes.
  if(N%size!=0) N_p=(N+size-N%size)/size;
  else N_p=N/size;
  if(N_p%2!=0) ++N_p ;
  /* std::vector<std::vector<double>> paths = pathsfinder(S0,E,r,sigma,T,N_p,M,rank,threads); */

  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when(N_p,M);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(N_p);
  for(int n=0;n<N_p;++n) exercise_st[n] = payoff(paths[M][n],E,payoff_fun);

  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.
  for(int m=M-1;m>0;--m){
    std::vector<double> x(N_p,-1);
    std::vector<double> y(N_p,-1);
    double sum_x; double sum_x2; double sum_x3; double sum_x4; double sum_y; double sum_yx; double sum_yx2;
    double sum_x_p = 0; double sum_x2_p = 0; double sum_x3_p = 0; double sum_x4_p = 0; double sum_y_p = 0; double sum_yx_p = 0; double sum_yx2_p = 0;

    double x_length; double x_length_p=0;

    // save first non-zero payoff incase there is only one path in the money.
    double fst_po;
    double fst_y;
    int fst_n;
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) private(E,r,dt) nowait reduction(+:sum_x_p,sum_x2_p,sum_x3_p,sum_x4_p,sum_y_p,sum_yx_p,sum_yx2_p,x_length_p)
    for(int n=0;n<N_p;++n){
      double payoff_val = payoff(paths[m][n],E,payoff_fun);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length_p;
        // stock price at time t_m
        double exer = paths[m][n];
        x[n] = exer;
        double cont = exp(-r*dt*(exercise_when[n]-m))*payoff(paths[exercise_when[n]][n],E,payoff_fun);
        y[n] = cont;

        // calc values for xTx and xTy.
        sum_x_p   += exer;
        sum_x2_p  += exer*exer;
        sum_x3_p  += exer*exer*exer;
        sum_x4_p  += exer*exer*exer*exer;
        sum_y_p   += cont;
        sum_yx_p  += cont*exer;
        sum_yx2_p += cont*exer*exer;
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

    MPI_Reduce(&sum_x_p,  &sum_x,  1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x2_p, &sum_x2, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x3_p, &sum_x3, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_x4_p, &sum_x4, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_y_p,  &sum_y,  1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx_p, &sum_yx, 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sum_yx2_p,&sum_yx2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    /* MPI_Reduce(&x_length_p,&x_length,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); */
    //
    // Use allreduce, becase if some process continues, then following bcast will deadlock.
    // if every x_length_p==0 then collectively go to next iteration.
    MPI_Allreduce(&x_length_p,&x_length,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    std::vector<double> coef;
    
    // if no path was in the money, skip it, because we are not interested in it.
    // when M is big and dt is small, the step m=1 might not be in money.
    if (x_length==0){ 
      continue;
    } else if(x_length==1){
      // if only 1 paths in the money, then compare current and discounted price.
      if(x_length_p) continue;
      else if(fst_po>fst_y){
        exercise_when[fst_n] = m;
        exercise_st[fst_n] = fst_po;
      };
      continue;
    }else if(x_length==2){
      // if only 2 paths in the money, then do linear regression
      // compose xTx and xTy
      coef.resize(2);
      if(rank==0){
        xTx.resize(2); xTy.resize(2);
        for(int i=0;i<2;++i) xTx[i].resize(2);
        xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
        xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
        xTy[0]    = sum_y   ; xTy[1]    = sum_yx;
      };
    }else if(x_length>2){
      coef.resize(3);
     // more than 2 paths in the money
      if(rank==0){
        xTx.resize(3); xTy.resize(3);
        for(int i=0;i<3;++i) xTx[i].resize(3);
        xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
        xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
        xTx[2][0] = sum_x2  ; xTx[2][1] = sum_x3; xTx[2][2] = sum_x4 ;
        xTy[0]    = sum_y   ; xTy[1]    = sum_yx; xTy[2]    = sum_yx2;
      };
    };
    
    if(rank==0){
      coef = mat_vec_mul(inverse(xTx),xTy);
    };


    MPI_Bcast(coef.data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
    for(int i=0;i<N_p;++i){
      if(x[i]!=-1){
        double EYIX = coef[0] + coef[1]*x[i] + coef[2]*pow(x[i],2);
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
#pragma omp for schedule(dynamic,100) nowait reduction(+:result_p)
  for(int n=0;n<N_p;++n){
    if(exercise_st[n]!=0) result_p+=exp(-r*exercise_when[n]*dt)*exercise_st[n];
  };
  }

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
  int threads =             getArg(argv,9);
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
  double result = mc_amer(S0,E,r,sigma,T,N,M,payoff_fun_d,size,rank,threads);
  auto end = std::chrono::system_clock::now();
  
  // close processes
  MPI_Finalize();
  auto end_overall = std::chrono::system_clock::now();

  if(rank==0){
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::chrono::duration<double> elapsed_seconds_overall = end_overall-start_overall;
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
        ,M
        );
  };
  return EXIT_SUCCESS;
}
