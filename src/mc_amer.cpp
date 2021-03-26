
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


std::vector<std::vector<double>> pathsfinder
(
 int N,int M,double S0
)
{
  std::vector<std::vector<double>> paths;
  for(int n=0;n<N;++n){

    // for each path create a new generator
    time_t cur_time;
    std::random_device rd{};
    std::mt19937 gen{rd()};
    gen.seed(time(&cur_time)+n);
    std::normal_distribution<> norm{0,sqrt(T)};

    // init new path
    std::vector<double> path;
    path.push_back(S0);

    // fill path
    for(int m=1;m<M+1;++m){
      path.push_back(path[m-1]+norm(gen));

    };
    // add path to paths vector.
    paths.push_back(path);
  };

  return paths;
}

std::vector<std::vector<double>> transpose
(
 std::vector<std::vector<double>> y
)
{
  std::vector<std::vector<double>> transposed;
  for(int j=0;j<y[0].size();++j){
    std::vector<double> tmp;
    for(int i=0;i<y.size();++i){
      tmp.push_back(y[i][j]);
    };
    transposed.push_back(tmp);
  };
  return transposed;
}

std::vector<std::vector<double>> matmul
(
  std::vector<std::vector<double>> x,
  std::vector<std::vector<double>> y
)
{
  std::vector<std::vector<double>> matmul;
  std::vector<std::vector<double>> yT = transpose(y);
  for(int i=0;i<x.size();++i){
    std::vector<double> tmp_vec;
    for(int k=0;k<x.size();++k){
      double tmp=0;
      for(int j=0;j<x[i].size();++j){
        tmp+=x[i][j]*yT[k][j];
      };
      tmp_vec.push_back(tmp);
    };
    matmul.push_back(tmp_vec);
  };
  return matmul;
}

std::vector<std::vector<double>> inverse
(
 std::vector<std::vector<double>> x
)
{
  std::vector<std::vector<double>> inversed;
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
      tmp.push_back(((x[(j+1)%3][(i+1)%3] * x[(j+2)%3][(i+2)%3]) - (x[(j+1)%3][(i+2)%3] * x[(j+2)%3][(i+1)%3]))/determinant);
    };
    inversed.push_back(tmp);
   };
  return inversed;
}


std::vector<double> coefficients
(
  std::vector<double> x, std::vector<double> y, int N
)
{
  std::vector<double> coef;

  std::vector<std::vector<double>> x_mat_transposed;

  // fill x_mat^T
  std::vector<double> ones(N,1);
  std::vector<double> x_mat_trans_tmp;

  // fill with ones
  x_mat_transposed.push_back(ones);
  // fill with x
  x_mat_transposed.push_back(x);

  // fill with x^2
  for(int n=0;n<N;++n){
    /* std::vector<double> x_mat_tmp; x_mat_tmp.push_back(1); x_mat_tmp.push_back(x[n]); x_mat_tmp.push_back(pow(x[n],2)); x_mat.push_back(x_mat_tmp); */
    x_mat_trans_tmp.push_back(pow(x[n],2));
  };
  x_mat_transposed.push_back(x_mat_trans_tmp);


  std::vector<std::vector<double>> x_mat = transpose(x_mat_transposed); // make x matrix
  std::vector<std::vector<double>> xTx = matmul(x_mat_transposed,x_mat); // calc x^T * x
  std::vector<std::vector<double>> y_mat; y_mat.push_back(y); //make y vector into a matrix
  std::vector<std::vector<double>> xTy = matmul(x_mat_transposed,y_mat); // calc x^T*y
  std::vector<std::vector<double>> xTx_inv = inverse(xTx); // calc (x^T*x)^{-1}
  std::vector<std::vector<double>> coef_mat = matmul(xTx_inv,xTy); // calculate coefficient matrix

  for(int i=0;i<coef_mat.size();++i) coef.push_back(coef_mat[i][0]);//extract coefficients to a vector;
  
  return coef;

}

double mc_amer
 (
 int N,int M,double S0,double E,double r, double T,double sigma
 )
{
  double dt = T/M;
  double result = 0;
  // calculate paths
  std::vector<std::vector<double>> paths = pathsfinder(N,M,S0);

  // store each paths timestep value when option is exercised
  std::vector<double> exercise_when;
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st;
  
  // Calculate timestep M payoff values for each path.
  // Store corresponding when and st values.
  for(int n=0;n<N;++n){
    double exercise = payoff(paths[n].back(),E);
    if(exercise>0){
      exercise_when.push_back(M);
      exercise_st.push_back(exercise);
    } else{
      exercise_when.push_back(0);
      exercise_st.push_back(0);
    };
  };

  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.
  for(int m=M-1;m>0;--m){
    std::vector<double> x;
    std::vector<double> y;
    for(int n=0;n<N;++n){
      // stock price at time t_m
      x.push_back(paths[n][m]);
      // discounted cashflow at time t_{m+1}
      y.push_back(exp(-r*dt)*payoff(paths[n][m+1],E));

    };
    for(int n=0;n<N;++n){
      std::vector<double> coef = coefficients(x,y,N);
      double EYIX = coef[0] + coef[1] * x[n] + coef[2]*pow(x[n],2);
      // exercise value at t_m
      double exercise = payoff(paths[n][m],E);
      if (exercise > EYIX) {
        exercise_when[n] = m;
        exercise_st[n] = exercise;
      };
    };
  };
  
  for(int n=0;n<N;++n){
    result+=exp(-r*exercise_when[n])*payoff(exercise_st[n],E);
  };

  return result/(double)N;
}

int main (int argc, char *argv[]){
  /* std::vector<std::vector<double>> mat= pathsfinder(10,10,100); */
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

  auto start = std::chrono::system_clock::now();
  double result = mc_amer(N,M,S0,E,r,T,sigma);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  reporting(
      "Serial",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N,
      M
      );
  return EXIT_SUCCESS;
}
