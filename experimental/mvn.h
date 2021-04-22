#include <Eigen/Dense>
#include <Eigen/LU>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    


// taken from https://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
// at 2021-04-22


/*
  We need a functor that can pretend it's const,
  but to be a good random number generator 
  it needs mutable state.
*/
namespace Eigen {
namespace internal {
template<typename Scalar> 
struct scalar_normal_dist_op 
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

/*
  Draw N samples from a assets-dimensional normal distribution
  with a specified mean and covariance
*/
Eigen::MatrixXd sample(int N,int assets,double sigma,double rho) 
/* void sample(int N,int assets,double sigma,double rho) */ 
{
  time_t cur_time;
  /* int assets = 2; // Dimensionality (rows) */
  /* int N=5;     // How many samples (columns) to draw */
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(time(&cur_time)); // Seed the rng with cur time 

  // Define mean and covariance of the distribution
  Eigen::VectorXd mean(assets);       
  Eigen::MatrixXd covar(assets,assets);

  for(int asset=0;asset<assets;++asset)
    mean(asset) = 0;
  for(int i=0;i<assets;++i){
    for(int j=0;j<assets;++j){
      if(i!=j) covar(i,j)=rho;
      else covar(i,j)=sigma;
    };
  };
  // mean << 0,0;
  /* covar <<  1, .5, */
  /*          .5,  1; */

  Eigen::MatrixXd normTransform(assets,assets);

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if 
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors() 
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform 
                           * Eigen::MatrixXd::NullaryExpr(assets,N,randN)).colwise() 
                           + mean;
  /* std::cout << "Mean\n" << mean << std::endl; */
  /* std::cout << "Covar\n" << covar << std::endl; */
  /* std::cout << "Samples\n" << samples.transpose() << std::endl; */
  return samples;
}

/* int main(int argc, char **argv) */
/* { */
/*   std::cout << sample(10,4,0.75,0.5) << std::endl; */
/*   return 0; */
/* } */

