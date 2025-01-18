#include "EigenSolver.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"
#include "eigen/Eigen/Sparse"
#include "ErrorHandler.h"
#include "Exception.h"

/**
 * @brief implementation of csmp-eigen interface solver class
 * @author Luat Khoa Tran
 */
 
using namespace std;
 
namespace csmp {

// Lut's original version
void EigenSolver::SolveMatrixEquation( SparseMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
  {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    // 0. initial checks
    if ( b.size() >= 10000U )
      csmp_error.Note( WARNING, "EigenSolver::SolveMatrixEquation",
                      "default setting is a direct solver and system has more than 10k DOF; solution may take time, consider using an iterative solver");
  
    // 1. convert sparse matrix to eigen sparse matrix
    // generic eigen sparse matrix is used
    Eigen::SparseMatrix<double> mat;
    const auto n_dof{ A.Rows() };
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

    // 2. convert rhs vector to eigen vector
    Eigen::VectorXd rhs;
    rhs.resize(A.Rows());
    for (size_t i = 0; i < x.size(); ++i)
      rhs[i] = b[i];

    // 3. solve using direct LU
    auto start = chrono::high_resolution_clock::now();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);
    cout <<"\n"<<"EigenSolver::SolveMatrixEquation (Sparse LU decomposition): solving "<< A.Rows() <<" x "<< A.Cols() <<" system..."<< endl;
    Eigen::VectorXd solution = solver.solve(rhs);
    auto stop = chrono::high_resolution_clock::now();

    // 4. copy back to ordinary vector
    for (size_t i(0); i < x.size(); ++i)
      x[i] = solution[i];
      
    cout <<"\n\t"<<"completed direct solution in "<< chrono::duration_cast<chrono::seconds>(stop-start).count() <<" seconds."<< endl;

} // end SolveMatrixEquation





// SKM modified for iterative solver (but slower than direct solver!)
/*
void EigenSolver::SolveMatrixEquation( SparseMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns )
  {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );
    
    // 1. convert sparse matrix to eigen sparse matrix
    // generic eigen sparse matrix is used
    Eigen::SparseMatrix<double> mat;
    const auto n_dof{ A.Rows() };
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

    // 2. convert rhs vector to eigen vector
    Eigen::VectorXd rhs;
    rhs.resize(A.Rows());
    for (size_t i = 0; i < x.size(); ++i)
      rhs[i] = b[i];

    // 3. solve system using BICGSTAB with preconditioning
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> >  solver;
    solver.preconditioner().setDroptol(1.0e-9);
    solver.preconditioner().setFillfactor(100);
    solver.compute(mat);
    Eigen::VectorXd solution = solver.solve(rhs);
    cout <<"\n"<<"EigenSolver::SolveMatrixEquation (BICGSTAB + preconditioner): solved "<< A.Rows() <<" x "<< A.Cols() <<" system."<< endl;
    cout << "\t"<<"#iterations:     " << solver.iterations() << std::endl;
    cout << "\t"<<"estimated error: " << solver.error()      << std::endl;
 
    // 4. copy back to ordinary vector
    for (size_t i(0u); i < x.size(); ++i)
      x[i] = solution[i];
      
} // end SolveMatrixEquation
*/





  

void EigenSolver::SolveMatrixEquation( CompressedRowMatrix& A,
                                       vector<double>& b,
                                       vector<double>& x,
                                       size_t no_unknowns)
{
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    // 0. initial checks
    if ( b.size() >= 20000U )
      csmp_error.Note( WARNING, "EigenSolver::SolveMatrixEquation",
                       "default setting is a direct solver and system has more than 20k DOF; solution may take time, consider using iterative solver");

    // 1. convert sparse matrix to eigen sparse matrix
    // generic eigen sparse matrix is used
    Eigen::SparseMatrix<double> mat;
    mat.resize(A.Rows(), A.Cols());
    mat.reserve(Eigen::VectorXi::Constant(A.Rows(), 24));
    for (size_t row{0U}; row < A.Rows(); ++row) {
      for (uint32_t index = A.ia[row]; index < A.ia[row + 1]; index++) {
        auto col = A.ja[index];
        auto value = A.a[index];
        mat.insert( row, col ) = value;
      }
    }
    mat.makeCompressed();

    // 2. conver rhs vector to eigen vector
    Eigen::VectorXd rhs;
    rhs.resize(A.Rows());
    for (size_t i = 0U; i < x.size(); ++i)
      rhs[i] = b[i];

    // 3. solve using direct LU
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);

    // 4. copy back to ordinary vector
    Eigen::VectorXd solution = solver.solve(rhs);
    for (size_t i(0); i < x.size(); ++i)
      x[i] = solution[i];
}


} // end csmp
