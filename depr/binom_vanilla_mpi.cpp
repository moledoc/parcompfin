
#include <example_eur.h>


double binom 
(
 double r, double sigma, double S0,
 double T, int N, double E,
 int size, int rank
 )
{
  double dt = T/(double)N;

  /* double beta = 0.5*(exp(-r*dt)+exp((r+sigma*sigma)*dt)); */
  /* double u = beta + sqrt(beta*beta-1); */
  /* double p = (exp(r*dt)-d)/(u-d); */
  /* double q = 1-p; */

  double u = exp(sigma*sqrt(dt));
  double d = 1/u;
  double R = exp(r*dt);
  double p = (R-d)/(u-d);
  double q = 1-p;

  std::vector<double> v_ij;
  int fallout=0;

  std::vector<double> v_ij_1;
  std::vector<double> v_ij_2;
  double result;
  MPI_Status status;
  MPI_Request request;
  if(rank!=0) result = 0;
  
  for (int n=N-1;n>0;--n){
    /* int n = N; */
    int n_begin = rank * (n+1)/(size-fallout);
    int n_end = (rank+1) * (n+1)/(size-fallout);
    if(rank != 0) n_begin-=1;
    if(n_end == v_ij.size()) v_ij.pop_back();

    /* std::cout << "rank: " << rank << " " << n_begin << ", " << n_end << std::endl; */

    if(n/(size-fallout)<2) ++fallout;
    if(rank == 0 || (rank!=0 && n/std::max(rank,1) > 2)){
      if (v_ij.size()==0){
        for(int i=n_begin;i<n_end;++i)
          v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E));
      };

      for(int i=0;i<n_end-n_begin-1;++i){
        v_ij[i] = (p*v_ij[i+1] + q*v_ij[i])/R;
      };
      v_ij.pop_back();
    };


    std::vector<double> tmp;
    if (rank == 0){
      tmp.resize(size);
    };
    MPI_Gather(&(v_ij.back()),1,MPI_DOUBLE,tmp.data(),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    /* std::cout << "rank: " << rank << " " << tmp.size() << std::endl; */
    if(rank==0){
      v_ij.push_back(tmp[1]);
      for (int p = size-1-fallout;p>1;--p){
        MPI_Send(&(tmp[p]),1,MPI_DOUBLE,p-1,0,MPI_COMM_WORLD);
      };
    } else if(rank < (size-1-fallout)) {
      double extra;
      MPI_Recv(&extra,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
    };
  };


  // igal juhul anna eelmisele.
  // kui eelmise n_end == vec.size(); siis ta viskab selle minema

  /* /1* /2* deal with odd N/p *2/ *1/ */
  /* /1* if (N%size==0 && rank == (size-1)) n_end -= 1; *1/ */

  /* /1* /2* First iteration *2/ *1/ */
  /* /1* for(int i=n_begin;i<n_end;++i){ *1/ */
  /* /1*   v_ij_1.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E)); *1/ */
  /* /1* }; *1/ */

  /* /1* v_ij_2.resize(N); *1/ */
  /* /1* // allgather gives me errors, use gather + bcast instead *1/ */
  /* /1* /2* MPI_Allgather(v_ij_2.data(),N+1,MPI_DOUBLE,v_ij_1.data(),N+1,MPI_DOUBLE,MPI_COMM_WORLD); *2/ *1/ */
  /* /1* MPI_Gather(v_ij_1.data(),N/size,MPI_DOUBLE,v_ij_2.data(),N/size,MPI_DOUBLE,0,MPI_COMM_WORLD); *1/ */
  /* /1* MPI_Bcast(v_ij_2.data(),N/size,MPI_DOUBLE,0,MPI_COMM_WORLD); *1/ */

  /* int fallout=0; */
  /* bool cond = rank == 0; */
  /* bool cond_1 = true; */

  /* for (int n=N;n>0;--n){ */
  /*   if (rank!=0) cond_1 = double(n)/(double)(size-fallout) >= 4; */
  /*   if (!cond_1) ++fallout; */
  /*   /1* if (fallout == size-rank) break; *1/ */
    
  /*   /1* if (n==1 && rank==0) {v_ij_1.clear();v_ij_1[0] = p*v_ij_2[1] + q*v_ij_2[0];break;}; *1/ */

  /*   /1* if (n/(rank+1) > 2){ *1/ */
  /*   /1* if (cond || cond_1){ *1/ */
  /*   if(cond_1){ */
  /*     n_begin = rank * n/size; */
  /*     n_end = (rank+1) * n/size; */

  /*     /1* deal with last process endpoint *1/ */
  /*     if (n%size==0 && rank == (size-1)) n_end -= 1; */
  /*     if (n==N && rank==(size-1)) n_end +=1; */
  /*   /1* std::cout << "n," << n<< "rank: " << rank << " n_begin " <<  n_begin << " n_end " <<  n_end << std::endl; *1/ */
  /*     /1* if (n%size!=0 && rank == (size-1)) n_end += 1; *1/ */

  /*     v_ij_1.clear(); */
  /*     if (n==N) { */
  /*       for(int i=n_begin;i<=n_end;++i){ */
  /*         v_ij_2.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E)); */
  /*       }; */
  /*     }; */
  /*   /1* std::cout << "rank: " << rank << " end2 " <<  n_end << std::endl; *1/ */
  /* /1* std::cout << "rank: " << rank << " " <<  v_ij_1.size() << std::endl; *1/ */
  /*     for(int i=n_begin;i<n_end;++i)//{ */
  /*       v_ij_1.push_back(p*v_ij_2[i+1] + q*v_ij_2[i]); */
  /*       /1* }; *1/ */
  /*     /1* }; *1/ */
  /*     v_ij_2.clear(); */
  /*     /1* v_ij_2.resize(n-1); *1/ */

  /* /1* std::cout << "n," << n<< " cond " << cond << " rank: " << rank << " fallout " <<  fallout << std::endl; *1/ */

  /*   }; */
  /*     double missing; */
  /*     double tmp[size]; */
  /*     /1* /2* tmp.resize(size); *2/ *1/ */
  /*     MPI_Gather(&(v_ij_1[0]),1,MPI_DOUBLE,&tmp,1,MPI_DOUBLE,0,MPI_COMM_WORLD); */

  /*     if(rank==0){ */
  /*       v_ij_1.push_back(tmp[1]); */
  /*       for(int i=1;i<size-fallout;++i) */
  /*         MPI_Isend(&(tmp[i]),1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&request); */
  /*     } else { */
  /*       if(cond_1) { */
  /*         MPI_Recv(&missing,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status); */
  /*         v_ij_1.push_back(missing); */
  /*       }; */
  /*     }; */
  /*     copy(v_ij_1.begin(),v_ij_1.end(),back_inserter(v_ij_2)); */


  /* }; */
  /* std::cout << "rank: " << rank << " " <<  v_ij_1.size() << std::endl; */

  /* if (rank == 0) result = v_ij_1[0]; */
  /* /1* if (rank == 0) result = p*v_ij_1[1] + q*v_ij_1[0]; *1/ */
  /* if (rank == 0){ */
  /*   auto it = v_ij_1.begin(); */
  /*   while(it!=v_ij_1.end()){ */
  /*     std::cout << (*it) << std::endl; */
  /*     it++; */
  /*   }; */
  /* }; */

  /* int v_ij_1_size = v_ij_1.size(); */
  /* std::cout << "rank: " << rank << " " <<  v_ij_1_size << std::endl; */

    /* for (int i=0;i<N+1;++i) */ 
    /*   std::cout << "rank: " << rank << "i: " << i << " :: "  <<  v_ij_1[i] << std::endl; */
  
  /* /1* Rest of the iterations *1/ */
  /* /1* for(int i=N-1;i>2*size;++i){ *1/ */
  /*   int i = N-1; */
    /* calc n ranges*/
    /* n_begin = rank * i/size; */
    /* n_end = (rank+1) * i/size; */
    /* /1* deal with odd N/p *1/ */
    /* if (N%size!=0 && rank == (size-1)) n_end += 1; */

    /* for(int i=n_begin;i<=n_end;++i){ */
    /*   v_ij_2.push_back(p*v_ij_1[i+1] + q*v_ij_1[i]); */
    /* }; */

    /* if(rank==0) v_ij_1_size = i/size +1; */
    /* else v_ij_1_size = v_ij_2.size(); */
    /* v_ij_1.clear(); */
    /* v_ij_1.resize(v_ij_1_size); */
    /* if(i%size!=0) v_ij_1_size = i/size + 1; */
    /* else v_ij_1_size = i/size; */
    /* v_ij_1.resize(v_ij_1_size); */
    /* std::cout << "rank: " << rank << ", i:  " << i << " " <<  v_ij_1_size << std::endl; */
    /* /1* break; *1/ */ 
    /* MPI_Allgather(v_ij_2.data(),i,MPI_DOUBLE,v_ij_1.data(),i,MPI_DOUBLE,MPI_COMM_WORLD); */
    /* v_ij_2.clear(); */
    /* break; */
  /* }; */

  /* ijf (rank == 0){ */
  /*   auto it = v_ij_1.begin(); */
  /*   while (it != v_ij_1.end()){ */
  /*     std::cout << (*it) << std::endl; */
  /*     ++it; */
  /*   }; */
  /* }; */
  /* std::cout<< "rank: " << rank << " "  << "kana" << std::endl; */


  /* do initial calc on N'th timestep */
  /*if (N%2==0 && (N%size!=0 || (N/2)%2!=0)){*/ 
  /*  for(int i=n_begin;i<=n_end;++i){*/
  /*    v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E));*/
  /*  };*/
  /*} else {*/
  /*  for(int i=n_begin;i<n_end;++i){*/
  /*    v_ij.push_back(payoff(S0*pow(u,i)*pow(d,N-i),E));*/
  /*  };*/
  /*};*/

  
  /*std::vector<double> v_ij_1 = v_ij_iter(v_ij,v_ij.size(),p,q);*/
  /*std::cout << "rank: " << rank << " v_ij_1: " << v_ij_1[0]<< std::endl;*/

  /*std::vector<double> v_ij_gather;*/
  /*if(rank==0) v_ij_gather.resize(4);*/
  /*// https://stackoverflow.com/questions/52784890/mpi-scatter-mpi-gather-with-stdvector*/
  /*MPI_Gather(v_ij_1.data(),1,MPI_DOUBLE,v_ij_gather.data(),1,MPI_DOUBLE,0,MPI_COMM_WORLD);*/
  /*if (rank == 0){*/
  /*  result = v_ij_iter(v_ij_gather,size,p,q)[0];*/
  /*  std::vector<double> tst = v_ij_iter(v_ij_gather,size,p,q);*/
  /*  std::cout << "suurus" << tst.size() << std::endl;*/
  /*};*/

  if(rank==0) return v_ij[0];
  else return 0;
}


int main (int argc, char *argv[]){
  auto start_overall = std::chrono::system_clock::now();
  int N = getArg(argv,1);
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
  
  if (size < 2){
    MPI_Finalize();
    std::cerr << "Needs at least 2 processes!" << std::endl;
  };

  /* bencmarking code found at: https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c */
  auto start = std::chrono::system_clock::now();
  double result = binom(r,sigma,S0,T,N,E,size,rank);
  auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::chrono::duration<double> elapsed_seconds_overall = end-start_overall;
  /* std::time_t end_time = std::chrono::system_clock::to_time_t(end); */

  /* std::cout << "finished computation at " << std::ctime(&end_time) */
  /*           << "elapsed time: " << elapsed_seconds.count() << "s\n"; */

  if (rank==0) std::cout << elapsed_seconds_overall.count() << "," << elapsed_seconds.count() << ","<< result << std::endl;
  MPI_Finalize();
  return EXIT_SUCCESS;
}
