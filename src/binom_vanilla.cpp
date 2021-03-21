
#include <example_eur.h>

double binom (
    double r, double sigma, double S0,
    double T, int N, double E
    )
{
  double dt = (double)T/(double)N;
  double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt));
  double u = beta + sqrt(beta*beta-1);
  /* double d = beta - sqrt(beta*beta-1); */
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;

  std::vector<double> v_ij;//(N+1,0);

  for(int i=0;i<N+1;++i)
    v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E));
    /* v_ij[i] = payoff(S0*pow(u,i)*pow(d,N-i),E); */

  for (int n=N;n>=0;--n){
    for(int i=0;i<n+1;++i){
      v_ij[i] = (p*v_ij[i+1] + q*v_ij[i])/R;
    };
  };

  /* std::vector<double> V_ij_1; */
  /* std::vector<double> V_ij_2; */

  /* for(int j=0;j<N+1;++j){ */
  /*   V_ij_2.push_back(payoff(S0*pow(u,j)*pow(d,N-j),E)); */
  /* }; */

  /* int until; */
  /* if (N%2==0) until = (N-1)/2+1; */
  /* else until = (N-1)/2; */

  /* for(int n=0;n<until;++n){ */
  /*   for(int i=0;i<V_ij_2.size()-1;++i){ */
  /*     V_ij_1.push_back(p*V_ij_2[i+1] + q*V_ij_2[i]); */
  /*   }; */
  /*   V_ij_2.clear(); */
  /*   if (V_ij_1.size() == 1){ */
  /*     V_ij_2 = V_ij_1; */
  /*     break; */
  /*   }; */
  /*   for(int i=0;i<V_ij_1.size()-1;++i){ */
  /*     V_ij_2.push_back(p*V_ij_1[i+1] + q*V_ij_1[i]); */
  /*   }; */
  /*   V_ij_1.clear(); */
  /* }; */

  /* double result = V_ij_2.front(); */
  /* return exp(-r*dt)*result; */
  /* return exp(-r*dt)*v_ij[0]; */
  return v_ij[0];
}


int main (int argc, char *argv[]){
  int N = getN(argv);
  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */

  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start;
  reporting(
      "Serial",
      elapsed_seconds_overall.count(),
      elapsed_seconds.count(),
      result,
      analytical,
      N
      );
  return EXIT_SUCCESS;
}
