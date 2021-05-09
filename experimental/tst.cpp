#include <common.h>

double dumb(int i){
  double result=0;
  for(int j=0;j<100;++j)
    result+=exp(log(i));
  return result;
}


int main(int argc,char *argv[]){
  double result;
  double result2;
  int count=100000;
  omp_set_num_threads(4);

  auto start = std::chrono::system_clock::now();
  result=0;
  result2=0;
  for(int i=0;i<count;++i)
    result+=dumb(i);

  for(int i=0;i<count;++i)
    result+=dumb(2*i);
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "serial: " << elapsed_seconds.count() <<", result: " << result2-result << std::endl;

  start = std::chrono::system_clock::now();
  result=0;
  result2=0;
#pragma omp parallel
  {
#pragma omp for schedule(dynamic,1000) nowait reduction(+:result)
  for(int i=0;i<count;++i){
    result+=dumb(i);
  };

#pragma omp for nowait //reduction(+:result)
  for(int i=0;i<count;++i){
    result+=dumb(2*i);
  };
  }
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "parallel for: " << elapsed_seconds.count() << ", result: " << result2-result <<std::endl;

  start = std::chrono::system_clock::now();
  result=0;
  result2=0;
#pragma omp sections
  {
#pragma omp section
  for(int i=0;i<count;++i)
    result+=dumb(i);

#pragma omp section
  for(int i=0;i<count;++i)
    result+=dumb(2*i);
  }
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "parallel sections: " << elapsed_seconds.count() << ", result: " << result2-result << std::endl;

  //---------------------------------------------------------------------------------------------
  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
  start = std::chrono::system_clock::now();
  result=0;
  result2=0;
  omp_set_nested(1);
#pragma omp parallel for collapse(1)
    for(int i=0;i<count/1000;++i)
      for(int j=0;j<count/1000;++j){
        result2+=dumb(2*j);
        result2-=dumb(j);
        result+=dumb(i);
      };
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "parallel not double for: " << elapsed_seconds.count() << ", result: " << result2-result << std::endl;
  
  start = std::chrono::system_clock::now();
  result=0;
  result2=0;
  omp_set_nested(1);
#pragma omp parallel for 
    for(int i=0;i<count/1000;++i)
#pragma omp parallel for
      for(int j=0;j<count/1000;++j){
        result2+=dumb(2*j);
        result2-=dumb(j);
        result+=dumb(i);
      };
  end = std::chrono::system_clock::now();
  elapsed_seconds = end-start;
  std::cout << "parallel nested: " << elapsed_seconds.count() << ", result: " << result2-result << std::endl;

}
