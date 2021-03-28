
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
 int N,int M,double S0
)
{
  // matrix to store paths
  std::vector<std::vector<double>> paths(N);
  for(int i=0;i<N;++i){
    paths[i].resize(M+1);
  };

  // make a generator from  N(0,sqrt(T))
  time_t cur_time;
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> norm{0,sqrt(T)};
    
  // generate paths
  for(int n=0;n<N;++n){
    /* // for each path use different seed */
    gen.seed(time(&cur_time)+n);

    // init new path
    /* std::vector<double> path(M+1); */
    paths[n][0] = S0;

    // fill path
    for(int m=1;m<M+1;++m){
      paths[n][m] = paths[n][m-1]+norm(gen);

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


std::vector<double> coefficients
(
  std::vector<double> x, std::vector<double> y
)
{
  /* std::vector<double> coef; */
    std::vector<double> coef(3,1);

  // calc (x^T*x) - analytically calculated the 3x3 matrix form

  /* double sum_x = 0; */
  /* double sum_x2 = 0; */
  /* double sum_x3 = 0; */
  /* double sum_x4 = 0; */
  /* double sum_y = 0; */
  /* double sum_yx = 0; */
  /* double sum_yx2 = 0; */

  /* std::cout << x.size() << std::endl; */
  /* std::for_each(x.begin(), x.end(), [&] (double xi) { */
  /*     sum_x += xi; */
  /*     sum_x2 += xi*xi; */
  /*     sum_x3 += xi*xi*xi; */
  /*     sum_x4 += xi*xi*xi*xi; */
  /* }); */

  /* for(int i=0;i<x.size();++i){ */
  /*   sum_x   += x[i]; */
  /*   sum_x2  += x[i]*x[i]; */
  /*   sum_x3  += x[i]*x[i]*x[i]; */
  /*   sum_x4  += x[i]*x[i]*x[i]*x[i]; */
  /*   sum_y   += y[i]; */
  /*   sum_yx  += y[i]*x[i]; */
  /*   sum_yx2 += y[i]*x[i]*x[i]; */
  /* }; */

  /* for(int i=0;i<x.size();++i){ */
  /*   double tmp = pow(x[i],2); */
  /*   sum_x += x[i]; */
  /*   /1* sum_x2 += pow(x[i],2); *1/ */
  /*   sum_x2 += tmp; */
  /*   /1* sum_x3 += pow(x[i],3); *1/ */
  /*   sum_x3 += tmp*x[i]; */
  /*   sum_x4 += pow(tmp,2); */
  /*   sum_y += y[i]; */
  /*   sum_yx += y[i]*x[i]; */
  /*   sum_yx2 += y[i]*tmp; */
  /* }; */

  /* std::vector<double> row(3); */
  /* std::vector<std::vector<double>> xTx; */
  /* std::vector<std::vector<double>> xTy; */

  /* row[0] = x.size(); row[1] = sum_x ; row[2] =  sum_x2; xTx.push_back(row); */
  /* row[0] =    sum_x; row[1] = sum_x2; row[2] =  sum_x3; xTx.push_back(row); */
  /* row[0] =   sum_x2; row[1] = sum_x3; row[2] =  sum_x4; xTx.push_back(row); */

  /* row[0] =    sum_y; row[1] = sum_yx; row[2] = sum_yx2; xTy.push_back(row); */
  /* std::vector<std::vector<double>> coef_mat = matmul(inverse(xTx),xTy); */

/* ----------------------------------------------- */
  /* std::vector<std::vector<double>> x_mat_transposed; */

  /* // fill x_mat^T */
  /* std::vector<double> x_mat_trans_tmp(x.size(),1); */
  /* // fill with ones */
  /* x_mat_transposed.push_back(x_mat_trans_tmp); */
  /* // fill with x */
  /* x_mat_transposed.push_back(x); */
  /* // fill with x^2 */
  /* for(int i=0;i<x.size();++i){ */
  /*   /1* std::vector<double> x_mat_tmp; x_mat_tmp.push_back(1); x_mat_tmp.push_back(x[n]); x_mat_tmp.push_back(pow(x[n],2)); x_mat.push_back(x_mat_tmp); *1/ */
  /*   x_mat_trans_tmp[i] = pow(x[i],2); */
  /* }; */
  /* x_mat_transposed.push_back(x_mat_trans_tmp); */


  /* std::vector<std::vector<double>> x_mat = transpose(x_mat_transposed); // make x matrix */
  /* std::vector<std::vector<double>> xTx = matmul(x_mat_transposed,x_mat); // calc x^T * x */
  /* std::vector<std::vector<double>> y_mat(1,y);// y_mat.push_back(y); //make y vector into a matrix */
  /* std::vector<std::vector<double>> xTy = matmul(x_mat_transposed,y_mat); // calc x^T*y */
  /* std::vector<std::vector<double>> xTx_inv = inverse(xTx); // calc (x^T*x)^{-1} */
  /* std::vector<std::vector<double>> coef_mat = matmul(xTx_inv,xTy); // calculate coefficient matrix */

  /* /1* matmul(x_mat_transposed,transpose(x_mat_transposed));// calc x^T* x *1/ */
  /* std::vector<std::vector<double>> y_mat(1,y); //make y vector into a matrix */
  /* /1* xTy = matmul(x_mat_transposed,y_mat); // calc x^T*y *1/ */
  /* /1* xTx_inv = inverse(matmul(x_mat_transposed,transpose(x_mat_transposed));); // calc (x^T*x)^{-1} *1/ */
  /* std::vector<std::vector<double>> coef_mat = matmul(inverse(matmul(x_mat_transposed,transpose(x_mat_transposed))),matmul(x_mat_transposed,y_mat)); // calculate coefficient matrix */
/* ----------------------------------------------- */


  /* for(int i=0;i<coef_mat.size();++i) coef.push_back(coef_mat[i][0]);//extract coefficients to a vector; */
  

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
  std::vector<double> exercise_when(N,0);
  // store each paths payoff value at timestep, when option is exercised. Value is 0 when it's not exercised
  std::vector<double> exercise_st(N,0);
  
  std::vector<std::vector<double>> xTx(3);
  for(int i=0;i<3;++i) xTx[i].resize(3);
  /* std::vector<std::vector<double>> xTy(1); */
  /* xTy[0].resize(3); */
  std::vector<double> xTy(3);

  // Calculate timestep M payoff values for each path.
  // Store corresponding when and st values.
  for(int n=0;n<N;++n){
    double payoff_val = payoff(paths[M][n],E);
    if(payoff_val>0){
      exercise_when[n] = M;
      exercise_st[n] = payoff_val;
      /* exercise_when.push_back(M); */
      /* exercise_st.push_back(payoff_val); */
    /* } else{ */
    /*   exercise_when.push_back(M); */
    /*   exercise_st.push_back(0); */
    };
  };

  // Find timesteps at each path when the option is exercised.
  // Store corresponding when and st values. Update them when earier exercise timestep is found.
  for(int m=M-1;m>0;--m){
    std::vector<double> x(N,0);
    std::vector<double> y(N,0);
    double sum_x = 0; double sum_x2 = 0; double sum_x3 = 0; double sum_x4 = 0; double sum_y = 0; double sum_yx = 0; double sum_yx2 = 0;

    for(int n=0;n<N;++n){
      double payoff_val = payoff(paths[m][n],E);

      // keep only paths that are in the money
      if(payoff_val>0){
        // stock price at time t_m
        double tmp_x = paths[m][n];
        /* x[n] = paths[m][n]; */
        x[n] = tmp_x;
        // discounted cashflow at time t_{m+1}
        double tmp_y = exp(-r*dt)*payoff(paths[m+1][n],E);
        /* y[n] = exp(-r*dt)*payoff(paths[m+1][n],E); */
        y[n] = tmp_y;

        sum_x   += tmp_x;
        sum_x2  += tmp_x*tmp_x;
        sum_x3  += tmp_x*tmp_x*tmp_x;
        sum_x4  += tmp_x*tmp_x*tmp_x*tmp_x;
        sum_y   += tmp_y;
        sum_yx  += tmp_y*tmp_x;
        sum_yx2 += tmp_y*tmp_x*tmp_x;
      };
    };

    // delete those elements, that are not in the money.
    remove(x.begin(),x.end(),0);
    remove(y.begin(),y.end(),0);
    // calculate coef
    // and calculate whether we continue or exercise

    xTx[0][0] = x.size(); xTx[0][1] = sum_x ; xTx[0][2] =  sum_x2;
    xTx[1][0] =    sum_x; xTx[1][1] = sum_x2; xTx[1][2] =  sum_x3;
    xTx[2][0] =   sum_x2; xTx[2][1] = sum_x3; xTx[2][2] =  sum_x4;
    xTy[0] =    sum_y; xTy[1] = sum_yx; xTy[2] = sum_yx2;

    std::vector<double> coef = mat_vec_mul(inverse(xTx),xTy);
    /* std::vector<double> coef(3,1); */

    for(int i=0;i<x.size();++i){


/* ---------------------------------------------------------------------------------------------------- */
  /* double sum_x  = 1; */
 /* double sum_x2  = 0; */
 /* double sum_x3  = 30; */
 /* double sum_x4  = 14; */
 /* double sum_y   = 0; */
 /* double sum_yx  = 60; */
 /* double sum_yx2 = 7; */

      /* for(int i=0;i<x.size();++i){ */
      /*   sum_x   += x[i]; */
      /*   sum_x2  += x[i]*x[i]; */
      /*   sum_x3  += x[i]*x[i]*x[i]; */
      /*   sum_x4  += x[i]*x[i]*x[i]*x[i]; */
      /*   sum_y   += y[i]; */
      /*   sum_yx  += y[i]*x[i]; */
      /*   sum_yx2 += y[i]*x[i]*x[i]; */
      /* }; */

      /* std::vector<std::vector<double>> xTx(3); */
      /* std::vector<double> row(3); */
      /* std::vector<std::vector<double>> xTx; */
      /* std::vector<std::vector<double>> xTy; */

      /* xTx[0][0] = x.size(); xTx[0][1] = sum_x ; xTx[0][2] =  sum_x2; */
      /* xTx[1][0] =    sum_x; xTx[1][1] = sum_x2; xTx[1][2] =  sum_x3; */
      /* xTx[2][0] =   sum_x2; xTx[2][1] = sum_x3; xTx[2][2] =  sum_x4; */

      /* /1* xTy[0][0] =    sum_y; xTy[0][1] = sum_yx; xTy[0][2] = sum_yx2; *1/ */
      /* xTy[0] =    sum_y; xTy[1] = sum_yx; xTy[2] = sum_yx2; */

      /* std::vector<std::vector<double>> coef_mat = mat_vec_mul(inverse(xTx),xTy); */
      /* for(int i=0;i<coef_mat.size();++i) coef[i] = coef_mat[i][0];//extract coefficients to a vector; */
/* ---------------------------------------------------------------------------------------------------- */


      /* std::vector<double> coef = coefficients(x,y); */
      double EYIX = coef[0] + coef[1] * x[i] + coef[2]*pow(x[i],2);
      // exercise value at t_m
      double exercise = payoff(paths[m][i],E);
      if (exercise > EYIX) {
        exercise_when[i] = m;
        exercise_st[i] = exercise;
      };
    };
  };
  
  for(int n=0;n<N;++n){
    if(exercise_st[n]!=0) result+=exp(-r*exercise_when[n])*exercise_st[n];
  };

  return result/(double)N;
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
