#include "EigenSolver.h"
#include "SparseMatrix.h"
#include "eigen/Eigen/Sparse"
#include "ErrorHandler.h"
#include "Exception.h"

/**
 * @brief implementation of csmp-eigen interface solver class
 * @author Luat Khoa Tran
 */
 
using namespace std;
 
namespace csmp {

void EigenSolver::SolveMatrixEquation( SparseMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
  {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    // 0. initial checks
    if ( b.size() >= 20000U )
      csmp_error.Note( WARNING, "EigenSolver::SolveMatrixEquation",
                      "default setting is a direct solver and system has more than 20k DOF; solution may take time, consider using iterative solver");
  
    // 1. convert sparse matrix to eigen sparse matrix
    // generic eigen sparse matrix is used
    Eigen::SparseMatrix<double> mat;
    const size_t n_dof{ A.Rows() };
    mat.resize(A.Rows(), A.Cols());
    mat.reserve(Eigen::VectorXi::Constant(A.Rows(), 24));
    for (size_t i{0U}; i < A.Rows(); ++i) {
        const auto rowEnd{ A.RowEnd(i) };
        for ( auto pair = A.RowBegin(i); pair!=rowEnd; ++pair ) {
#ifdef DEBUG
           if ( (*pair).first >= n_dof )
               csmp_error.Note( ERROR, "EigenSolver::SolveMatrixEquation", to_string( (*pair).first ),
                               "column index retrieved from the SparseMatrix is out of bound.");
#endif
               mat.insert( i, (*pair).first ) = (*pair).second;
            }
        }
    mat.makeCompressed();

    // 2. conver rhs vector to eigen vector
    Eigen::VectorXd rhs;
    rhs.resize(A.Rows());
    for (size_t i = 0; i < x.size(); ++i)
      rhs[i] = b[i];

    // 3. solve using direct LU
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);

    // 4. copy back to ordinary vector
    Eigen::VectorXd solution = solver.solve(rhs);
    for (size_t i(0); i < x.size(); ++i)
      x[i] = solution[i];
      
} // end SolveMatrixEquation
  
  
  
  
/**
Not implemented yet!
*/
  void EigenSolver::SolveMatrixEquation( CompressedRowMatrix& A,
                                         vector<double>& b,
                                         vector<double>& x,
                                         size_t no_unknowns)
  {
      throw csmp::Exception( ERROR, "EigenSolver::SolveMatrixEquation",
                             "compressed row-matrix based version not implemented yet");
  }


} // end csmp
