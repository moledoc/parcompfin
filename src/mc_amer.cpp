
#include <common.h>
#include <comparison.h>

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
  std::vector<std::vector<double>> paths = pathsfinder(S0,E,r,sigma,T,N,M);
  
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

    // save first/two non-zero payoff incase there is only one/two path in the money.
    std::vector<double> fst_po;
    std::vector<double> fst_y;
    std::vector<int> fst_n;
    for(int n=0;n<N;++n){
      double payoff_val = payoff(paths[m][n],E,payoff_fun);
      // keep only paths that are in the money
      if(payoff_val>0){
        ++x_length;
        // stock price at time t_m
        double exer = paths[m][n]-E;
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
        if(x_length<=2){
          fst_po.push_back(payoff_val);
          fst_y.push_back(cont);
          fst_n.push_back(n);
        };
      };
    };
    
    std::vector<std::vector<double>> xTx;
    std::vector<double> xTy;

    // if no path was in the money, skip it, because we are not interested in it.
    // when M is big and dt is small, the step m=1 might not be in money.
    if (x_length==0){ 
      continue;
    } else if(x_length<=2 && x_length>0){
      // if only 1-2 paths in the money, then compare current and discounted price.
      for(int i=0;i<fst_n.size();++i){
        if(fst_po[i]>fst_y[i]){
          exercise_when[fst_n[i]] = m;
          exercise_st[fst_n[i]] = fst_po[i];
        };
      };
      continue;
    }else if(x_length>2){
     // more than 2 paths in the money
      xTx.resize(3); xTy.resize(3);
      for(int i=0;i<3;++i) xTx[i].resize(3);
      xTx[0][0] = x_length; xTx[0][1] = sum_x ; xTx[0][2] = sum_x2 ;
      xTx[1][0] = sum_x   ; xTx[1][1] = sum_x2; xTx[1][2] = sum_x3 ;
      xTx[2][0] = sum_x2  ; xTx[2][1] = sum_x3; xTx[2][2] = sum_x4 ;
      xTy[0]    = sum_y   ; xTy[1]    = sum_yx; xTy[2]    = sum_yx2;
    };
    
    // calculate LSM coefficients
    std::vector<double> coef = mat_vec_mul(inverse(xTx),xTy);
    // save asset price and step index, when we use the options rights
    for(int i=0;i<N;++i){
      if(x[i]!=-1){
        double Y_hat = coef[0] + coef[1]*x[i] + coef[2]*pow(x[i],2);
        double payoff_val = payoff(x[i],E,payoff_fun);
        if (payoff_val > Y_hat) {
          exercise_when[i] = m;
          exercise_st[i] = payoff_val;
        };
      };
    };
  };

  for(int n=0;n<N;++n){
    if(exercise_st[n]!=0) result+=exp(-r*exercise_when[n]*dt)*exercise_st[n];
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
