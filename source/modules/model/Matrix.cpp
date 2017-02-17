#include "Matrix.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

using namespace std;

namespace csmp {

// operator=( DM )
// ---------------------------------------

Matrix& Matrix::operator=( const Matrix& mat )
 {
    if ( &mat != this ) {
         data.resize( mat.Rows(), vector<double64>(mat.Cols(),0.) );
         rows  = mat.rows;
         cols  = mat.cols;

         for ( size_t i=0U; i<rows; i++ )
           for ( size_t j=0U; j<cols; j++ ) 
             data[i][j] = mat.data[i][j]; 
      }
      
    return *this;
 }




Matrix& Matrix::operator=( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] = val; 
      
    return *this;
 }



// copy constructor
// ---------------------------------------

Matrix::Matrix( const Matrix& mat )
 : rows(mat.rows), cols(mat.cols), data(mat.data)
 {
 }



// constructor (i,j, value)
// ---------------------------------------

Matrix::Matrix( size_t m, size_t n, double64 val )
 : rows(m), cols(n), data(m,vector<double64>(n,val))
 {
 }


// destructor
// ---------------------------------------

Matrix::~Matrix()
 {
 }





Matrix&  Matrix::operator*=( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] *= val;
       
    return *this;
 }



// matrix - tensor multiplication
// ------------------------------ 

Matrix& 
  Matrix::operator*=( const TensorVariable<3U>& ts )
 {
    if ( cols != 3U ) {
         std::cout <<"\nMatrix<double>::operator*=(TensorVar.3D): matrix.cols != tensor.dim\n";
         throw std::range_error("Matrix<>::operator*=(TensorVar.3D): matrix.cols != tensor.dim");
      }
      
    // when loop unrolling is not possible
    if ( rows != cols ) {
         // C++ named return value optimization
         Matrix  temp( rows, 3U );
    
         for ( size_t i=0U; i<rows; i++ )
           for ( size_t j=0U; j<3U; j++ ) {
                 temp.data[i][j] = static_cast<double64>(0.0);
                 for ( size_t k=0U; k<3U; k++ ) 
                   temp.data[i][j] += data[i][k] * ts(k,j);
             }
         return *this = temp;
      }

    TensorVariable<3U> temp;
    // Unrolled loop(i...3,j...3)
    // row 0
    temp(0,0) = data[0][0] * ts(0,0) + data[0][1] * ts(1,0) + data[0][2] * ts(2,0);
    temp(0,1) = data[0][0] * ts(0,1) + data[0][1] * ts(1,1) + data[0][2] * ts(2,1);
    temp(0,2) = data[0][0] * ts(0,2) + data[0][1] * ts(1,2) + data[0][2] * ts(2,2);
     // row 1
    temp(1,0) = data[1][0] * ts(0,0) + data[1][1] * ts(1,0) + data[1][2] * ts(2,0);
    temp(1,1) = data[1][0] * ts(0,1) + data[1][1] * ts(1,1) + data[1][2] * ts(2,1);
    temp(1,2) = data[1][0] * ts(0,2) + data[1][1] * ts(1,2) + data[1][2] * ts(2,2);
    // row 2
    temp(2,0) = data[2][0] * ts(0,0) + data[2][1] * ts(1,0) + data[2][2] * ts(2,0);
    temp(2,1) = data[2][0] * ts(0,1) + data[2][1] * ts(1,1) + data[2][2] * ts(2,1);
    temp(2,2) = data[2][0] * ts(0,2) + data[2][1] * ts(1,2) + data[2][2] * ts(2,2);
  
    return *this = temp; 
 }




// operator*=  matrix multiplication
// ---------------------------------------
/// @warning creates a temporary

Matrix&  Matrix::operator*=( const Matrix& mat )
 {
    if ( cols != mat.rows ) {
         std::cerr <<"\nMatrix<double>::operator*=: Matrices cannot be multiplied ";
         std::cerr <<"because of incompatible sizes (A *= B, see matrices below): "<< endl;
         Out();
         mat.Out();
         throw length_error("Matrix<mn_max>::operator*=");
      }

    //cout << "cols vs. rows."<<cols<<" "<<rows<<endl;
    // C++ named return value optimization
    //Matrix  temp( rows, mat.cols );
    Matrix  temp( cols, mat.rows ,0.0);

    for ( size_t i=0U; i<cols; i++ )
    {
      for ( size_t j=0U; j<mat.rows; j++ ) {
        //temp.data[i][j] = static_cast<double64>(0.0);
        for ( size_t k=0U; k<mat.rows; k++ )
        {
          temp.data[i][j] += this->data[i][k] * mat.data[k][j];
        }
      }
    }

    *this = temp;  
    //(*this).Out(3);
    return *this;
 }






Matrix& 
Matrix::operator*=( const TensorVariable<2U>& ts )
 {
    if ( cols != 2U ) {
         std::cout <<"\nMatrix<double>::operator*=(TensorVar.2D): matrix.cols != tensor.dim\n";
         throw std::range_error("Matrix<>::operator*=(TensorVar.2D): matrix.cols != tensor.dim");
      }
      
    // when loop unrolling is not possible
    if ( rows != cols ) {
         // C++ named return value optimization
         Matrix  temp( rows, 2U );
    
         for ( size_t i=0U; i<rows; i++ )
           for ( size_t j=0U; j<2U; j++ ) {
                 temp.data[i][j] = static_cast<double64>(0.0);
                 for ( size_t k=0U; k<2U; k++ ) 
                   temp.data[i][j] += data[i][k] * ts(k,j);
             }
         return *this = temp;
      }

    TensorVariable<3U> temp;
    // temp(i,j) += data[i][k]*ts(k][j];
    temp(0,0) = data[0][0] * ts(0,0) + data[0][1] * ts(1,0);
    temp(0,1) = data[0][0] * ts(0,1) + data[0][1] * ts(1,1);
    temp(1,0) = data[1][0] * ts(0,0) + data[1][1] * ts(1,0);
    temp(1,1) = data[1][0] * ts(0,1) + data[1][1] * ts(1,1);
       
    return *this = temp;
 }
 


// RowCondenseTo
// ---------------------------------------
/// vec = Mat * unity vector

void Matrix::RowCondenseTo( std::vector<double64>& vec ) const
 {
    vec.resize(rows);
    vector<double64>( vec ).swap( vec );
    for ( size_t i=0; i<rows; i++ ) vec[i] = RowSum(i);
 } 



// operator+=
// ---------------------------------------

Matrix& Matrix::operator+=( const Matrix& mat )
 {
    CheckSizes( mat, "Matrix::operator+=" );
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] += mat.data[i][j]; 
      
    return *this;
 }
 



// operator-=
// ---------------------------------------

Matrix& Matrix::operator-=( const Matrix& mat )
 {
    CheckSizes( mat, "Matrix::operator-=" );
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] -= mat.data[i][j]; 
      
    return *this;
 }





// operator+  
// ---------------------------------------

Matrix  Matrix::operator+( const Matrix& mat ) const
 {
    CheckSizes( mat, "Matrix::operator+" );
    Matrix  temp(cols,rows);
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        temp.data[i][j] = data[i][j] + mat.data[i][j]; 
      
    return temp;
 }

 
// operator-  
// ---------------------------------------

Matrix Matrix::operator-( const Matrix& mat ) const
 {
    CheckSizes( mat, "Matrix::operator-" );
    Matrix  temp(cols,rows);
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        temp.data[i][j] = data[i][j] - mat.data[i][j]; 
      
    return temp;
 }




Matrix& 
  Matrix::operator*=( const VectorVariable<2U>& vc )
 {   
    assert ( cols == 2U );
    double64 sum;
    
    // when loop unrolling is not possible
    if ( rows != cols ) {
        for ( size_t i=0; i<rows; i++ ) {
               sum = static_cast<double64>(0.0);
               for ( size_t j=0; j<cols; j++ ) 
                 sum +=  data[i][j] * vc[j];
               data[i][0] = sum;
            }
         cols = 1U;
         return *this;
      }
    
    sum = data[0][0] * vc[0] + data[0][1] * vc[1];
    data[0][0] = sum;
    sum = data[1][0] * vc[0] + data[1][1] * vc[1];
    data[1][0] = sum;
    
    cols = 1U;

    return *this;
 }




Matrix& 
  Matrix::operator*=( const VectorVariable<3U>& vc )
 {
    assert ( cols == 3U );
    double64 sum;
    
    // when loop unrolling is not possible
    if ( rows != cols ) {
        for ( size_t i=0; i<rows; i++ ) {
               sum = static_cast<double64>(0.0);
               for ( size_t j=0; j<cols; j++ ) 
                 sum +=  data[i][j] * vc[j];
               data[i][0] = sum;
            }
         cols = 1U;
         return *this;
      }
    
    sum = data[0][0] * vc[0] + data[0][1] * vc[1] + data[0][2] * vc[2];
    data[0][0] = sum;
    sum = data[1][0] * vc[0] + data[1][1] * vc[1] + data[1][2] * vc[2];
    data[1][0] = sum;
    sum = data[2][0] * vc[0] + data[2][1] * vc[1] + data[2][2] * vc[2];
    data[2][0] = sum;

    cols = 1U;

    return *this;
 }




// operator*  matrix multiplication
// ---------------------------------------

Matrix  Matrix::operator*( const Matrix& mat ) const
 {
//    cout << "cols vs. rows."<<cols<<" "<<rows<<endl;
    if ( cols != mat.rows ) {
         std::cout <<"\nMatrix<double>::operator*: Matrices cannot be multiplied "; 
         std::cout <<"because of incompatible sizes (A *= B, see matrices below): "<< std::endl;
         Out(3);
         mat.Out(3);
         throw std::length_error("Matrix<mn_max>::operator*");
      }

    // use STANDARD return value optimization
    Matrix  temp( rows, mat.cols );
    //temp.Out(3);
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<mat.cols; j++ ) {
           temp.data[i][j] = static_cast<double64>(0.0);
           for ( size_t k=0U; k<mat.rows; k++ ) 
             temp.data[i][j] += data[i][k]*mat.data[k][j];
        }  
    return temp;
 }




// operator*  matrix with STL vector<double64> multiplication
// ------------------------------------------------

Matrix& Matrix::operator*=( const std::vector<double64>& vecT )
 {
    if ( cols != vecT.size() ) {
         std::cout <<"\nMatrix<double>::operator*=: ";
         std::cout <<"Matrix and vector<double64> cannot be multiplied "; 
         std::cout <<"because of incompatible sizes (A(cols != vecT(size)): "<< vecT.size() << std::endl;
         throw std::length_error("Matrix<mn_max>::operator*=");
      }
    
    size_t  i, j;
    double64            sum;
    
    for ( i=0; i<rows; i++ ) {
         for ( sum=static_cast<double64>(0.0), j=0; j<vecT.size(); j++ ) sum += data[i][j] * vecT[j];
         data[i][0U] = sum;
      }
      
    // resizing the matrix
    cols = 1U;  
      
    return *this;
 }


// operator*  matrix with C-array multiplication
// ---------------------------------------------
/// Clearly, all hell will break loose if the C-array has not been initialized properly

Matrix& Matrix::operator*=( const double64* vecT ) 
 {
    if ( vecT == NULL ) {
         std::cout <<"\nMatrix<double>::operator*=: Matrix and vector<double64> cannot be multiplied "; 
         std::cout <<"because vecT C-array is not initialized." << std::endl;
         throw std::length_error("Matrix<mn_max>::operator*=");
      }
    size_t  i, j;
    double64      sum;
    
    for ( i=0U; i<rows; i++ ) {
         for ( sum=static_cast<double64>(0.0), j=0; j<rows; j++ ) sum += data[i][j] * vecT[j];
         data[i][0] = sum;
      }
      
    // resizing the matrix
    cols = 1;  
      
    return *this;
 }
    



/// RES = A B^T
// ---------------------------------------

void 
Matrix::MultiplyWithTransposedOf( const Matrix& B, Matrix& RES ) const
 {
    if ( cols != B.cols ) {
         cout <<"\nMatrix<double>::MultiplyWithTransposedOf: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out();
         B.Out();
         throw length_error("Matrix<mn_max>::MultiplyWithTransposedOf");
      }
    
    RES.Resize( rows, B.rows );
    double64  sum;
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<B.rows; j++ ) {
           sum = static_cast<double64>(0.0);
           for ( size_t k=0U; k<cols; k++ ) 
             sum += data[i][k] * B.data[j][k];
           RES.data[i][j] = sum;  
        }
 }




/// RES = A^T B
// ---------------------------------------

void 
Matrix::MultiplyTransposedOfWith( const Matrix& B, Matrix& RES ) const
 {
    if ( rows != B.rows ) {
         cout <<"\nMatrix<double>::MultiplyTransposedOfWith: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out();
         B.Out();
         throw length_error("Matrix<mn_max>::MultiplyTransposedOfWith");
      }
    
    RES.Resize( cols, B.cols );
    double64  sum;
    
    for ( size_t i=0U; i<cols; i++ )
      for ( size_t j=0U; j<B.cols; j++ ) {
           sum = static_cast<double64>(0.0);
           for ( size_t k=0U; k<B.rows; k++ ) 
             sum += data[k][i] * B.data[k][j];
           RES.data[i][j] = sum;
        }
 }




/// RES = A^T A
// ---------------------------------------

void 
Matrix::TransposedProduct( Matrix& RES ) const
 {
    RES.Resize( cols, cols );
    double64  sum;
    
    for ( size_t i=0U; i<cols; i++ )
      for ( size_t j=0U; j<cols; j++ ) {
           sum = static_cast<double64>(0.0);
           for ( size_t k=0U; k<rows; k++ ) 
             sum += data[k][i] * data[k][j];
           RES.data[j][i] = sum;  
        }
 }




/// unscaled L_infinity matrix norm ||A||inf = max_i * sum_j |aij|

double64 Matrix::NormL_Infinity() const
 {
	if ( rows == 0 || cols == 0 ) {
         cerr <<"\nMatrix<double>::NormL_Infinity: empty matrix!\n";
         return std::numeric_limits<double64>::quiet_NaN();
      }

	double64  maxval(static_cast<double64>(0.0)), sum;

	for ( size_t i=0; i<rows; i++ ) {
		 sum = static_cast<double64>(0.0);
		 for ( size_t j=0; j<cols; j++ )
		   sum += fabs(data[i][j]);
		 maxval = max(maxval,sum);
	  }

	return maxval;
 }

/*

// operator +
// ---------------------------------------

Matrix  operator+( const Matrix& a, const Matrix& b )
 {
    Matrix temp(a); 
    temp += b;
    return temp; 
 }


// operator -
// ---------------------------------------

Matrix  operator-( const Matrix& a, const Matrix& b )
 {
    Matrix temp(a);  
    temp -= b;
    return temp; 
 }


// operator *
// ---------------------------------------

Matrix  operator*( const Matrix& a, const Matrix& b )
 {
#ifndef NDEBUG
    if ( a.Cols() != b.Rows() ) {
         cout <<"\nMatrix<double> operator*: Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A * B, see matrices below): "<< endl;
         a.Out(3L);
         b.Out(3L);
         throw length_error("Matrix<mn_max>::operator*");
      }
#endif
    Matrix  temp( a.Rows(), b.Cols() );
    
    for ( size_t i=0U; i<a.Rows(); i++ )
      for ( size_t j=0U; j<b.Cols(); j++ ) {
           temp(i,j) = static_cast<double64>(0.0);
           for ( size_t k=0U; k<b.Rows(); k++ ) 
             temp(i,j) += a(i,k) * b(k,j);
        }
      
    return temp;
 }

*/

/// unscaled L1 matrix norm ||A1|| = max_j * sum_i |aij|

double64 Matrix::NormL1() const
 {
	if ( rows == 0 || cols == 0 ) {
         cerr <<"\nMatrix<double>::NormL1: empty matrix!\n";
         return std::numeric_limits<double64>::quiet_NaN();
      }

	double64  maxval(static_cast<double64>(0.0)), sum;

	for ( size_t j=0; j<cols; j++ ) {
		 sum = static_cast<double64>(0.0);
		 for ( size_t i=0; i<rows; i++ )
		   sum += fabs(data[i][j]);
		 maxval = max(maxval,sum);
	  }

	return maxval;
 }


// -------------------------------------------------------------------------------------------
// OPERATOR FUNCTIONS
// -------------------------------------------------------------------------------------------

/// vector^T = Matrix * vector^T
// ---------------------------------------

vector<double64>  operator*( const Matrix& mat, const vector<double64>& vec )
 {
    assert ( mat.Cols() == vec.size() );
    
    vector<double64> temp(mat.Rows(),static_cast<double64>(0.0));
    
    for ( size_t i=0; i<mat.Rows(); i++ )
      for ( size_t j=0; j<mat.Cols(); j++ ) temp[i] += mat(i,j) * vec[j];
    
    return temp;
      
 } // end operator*



/// Matrix = vector^T * Matrix
// ---------------------------------------

Matrix  operator*( const vector<double64>& vec, const Matrix& mat )
 {
    if ( vec.size() != mat.Rows() ) {
         cout <<"\noperator*: vector<double64> cannot be multiplied with matrix"; 
         cout <<"because of incompatible sizes (v * M): "<< endl;
         throw length_error("Matrix<mn_max>::operator*");
      }
    Matrix  temp( vec.size(), mat.Cols() );
    
    for ( size_t i=0U; i<vec.size(); i++ )
      for ( size_t j=0U; j<mat.Cols(); j++ ) {
           temp(i,j) = static_cast<double64>(0.0);
           for ( size_t k=0U; k<mat.Rows(); k++ ) 
             temp(i,j) += vec[k] * mat(k,j);
        }
      
    return temp;
    
 } // end operator




// Resize()
// ---------------------------------------------

void Matrix::Resize( size_t m, size_t n )
 {
    // resize matrix but keep storage as is if the new matrix 
    // is smaller than the old matrix
    if ( m <= rows && n <= cols ) {
        rows = m;
        cols = n;        
        return;
      }
    
    // increase matrix size
    rows = m;
    cols = n;    
    data.resize(m);
    vector<vector<double64> >( data ).swap( data );
    for ( size_t i=0U; i<m; i++ ) {
         data[i].resize(n);
         vector<double64>( data[i] ).swap( data[i] );
      }
 }


// Identity()
// ---------------------------------------------

void Matrix::Identity()
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        if ( i == j ) data[i][j] = static_cast<double64>(1.0);
        else          data[i][j] = static_cast<double64>(0.0);         
 }
 

// Zero()
// ---------------------------------------------

void Matrix::Zero()
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        data[i][j] = static_cast<double64>(0.0);
 }
 

// ZeroRow()
// ---------------------------------------------

void Matrix::ZeroRow( size_t row )
 {
    if ( row < rows ) for ( size_t j=0; j<cols; j++ ) 
      data[row][j] = static_cast<double64>(0.0);
      
    else {
        std::cout <<"\nMatrix<double>::ZeroRow: ";
        std::cout <<"Target rows does not exist: "<< row << std::endl;
        throw std::range_error("Matrix<double64>::ZeroRow");
     }
 }


// ZeroCol()
// ---------------------------------------------

void Matrix::ZeroCol( size_t col )
 {
    if ( col < cols ) for ( size_t i=0; i<rows; i++ ) data[i][col] = static_cast<double64>(0.0);
    else
    std::cout <<"\nMatrix<double>::ZeroCol: Target column does not exist: "<< col << std::endl;
 }



// Fill()
// ---------------------------------------------

void Matrix::Fill( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] = val;
 }


// FillRow()
// ---------------------------------------------

void Matrix::FillRow( size_t row, double64 val )
 {
    if ( row < rows ) for ( size_t j=0U; j<cols; j++ ) data[row][j] = val;
    else {
      std::cout <<"\nMatrix<double>::FillRow: ";
      std::cout <<"Target rows does not exist: "<< row << std::endl;
     }
 }
 
 
// FillCol()
// ---------------------------------------------

void Matrix::FillCol( size_t col, double64 val )
 {
    if ( col < cols ) for ( size_t i=0U; i<rows; i++ ) data[i][col] = val;
    else {
      std::cout <<"\nMatrix<double>::FillCol: ";
      std::cout <<"Target column does not exist: "<< col << std::endl;
     }
 }


// Transposed()
// ---------------------------------------------

void Matrix::Transposed( Matrix& M ) const
 {
    M.Resize(cols,rows);
 
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) M.data[j][i] = data[i][j];
 }




// RowSum()
// ---------------------------------------------

double64   Matrix::RowSum( size_t row ) const
 {
    double64  sum(0.);
 
    if ( row < rows ) {
         for ( size_t j=0U; j<cols; j++ ) sum += data[row][j];
         return sum;
      }

    std::cout <<"\nMatrix<double>::RowSum: ";
    std::cout <<"Target rows does not exist: "<< row << std::endl;
    throw std::length_error("Matrix<mn_max>::RowSum");

    return static_cast<double64>(0.0);
 }
 

// ColSum()
// ---------------------------------------------

double64   Matrix::ColSum( size_t col ) const
 {
    double64  sum(0.0);
 
    if ( col < cols ) {
         for ( size_t i=0U; i<rows; i++ ) sum += data[i][col];
         return sum;
      }

    std::cout <<"\nMatrix<double>::ColSum: ";
    std::cout <<"Target column does not exist: "<< col << std::endl;
    throw std::length_error("Matrix<mn_max>::ColSum");

    return static_cast<double64>(0.0);
 }



// In()
// ---------------------------------------

void Matrix::In()
 {
    cout <<"\nMatrix<double>: Please enter matrix entries as prompted for: "<< endl;
    cout <<"Enter number of rows and columns: ";
    cin >> rows >> cols;
    
    for ( size_t i=0; i<rows; i++ )
      {
         cout <<"\nEnter entries of row "<< i+1 <<": ";
         for ( size_t j=0; j<cols; j++ ) cin >> data[i][j];
      }

    cout <<"\nThank you."<< endl;    
 }



// Out( digits )
// ---------------------------------------

void Matrix::Out( long digits ) const
 {
    long  prec(cout.precision(digits));
    cout <<"\nMatrix<double>::Out(): m="<< rows <<", n="<< cols << endl;
    if ( digits != 0U ) cout.setf(ios::scientific);
    size_t row_break, split_adouble64er(10U);
     
    for ( size_t i=0; i<rows; i++ )
      {
         row_break = 1;
         for ( size_t j=0; j<cols; j++, row_break++ )
           {
              if ( data[i][j] >= 0. ) cout <<" ";
              cout << data[i][j] <<" ";
              if ( row_break == split_adouble64er )
                {
                   cout << endl;
                   row_break = 0U;
                }
           }
         cout << endl;
      }
      
    if ( digits != 0U ) {
         cout.unsetf( ios::scientific );
         cout.precision(prec);
      }
      
    cout << endl;    
    
 } // end Out()



// operator= ( TensorVariable 2D)
// ---------------------------------------
Matrix& Matrix::operator=( const TensorVariable<2U>& ts )
 {
    rows = cols = 2U;
    data[0][0] = ts(0,0);
    data[0][1] = ts(0,1);
    data[1][0] = ts(1,0);
    data[1][1] = ts(1,1);
      
    return *this; 
 }



// operator= ( TensorVariable 3D)
// ---------------------------------------
Matrix&  Matrix::operator=( const TensorVariable<3U>& ts )
 {
    rows = cols = 3U;
    data[0][0] = ts(0,0);
    data[0][1] = ts(0,1);
    data[0][2] = ts(0,2);
    data[1][0] = ts(1,0);
    data[1][1] = ts(1,1);
    data[1][2] = ts(1,2);
    data[2][0] = ts(2,0);
    data[2][1] = ts(2,1);
    data[2][2] = ts(2,2); 
      
    return *this; 
 }



/// default constructor
Matrix::Matrix()
 : rows(0), cols(0)
 {
 }


/// constructor (i,j)
Matrix::Matrix( size_t m, size_t n )
 : rows(m), cols(n), data(m,std::vector<double64>(n))
 {
 }





/// operator (i,j)
double64& Matrix::operator()( size_t m, size_t n )
 {
#ifndef NDEBUG
    CheckRange( m, n, "Matrix::operator()");
#endif
    return data[m][n];
 }



/// operator (i,j) const
const double64& Matrix::operator()( size_t m, size_t n ) const
 {
#ifndef NDEBUG
    CheckRange( m, n, "Matrix::operator()");
#endif
    return data[m][n];
 }





/// returns number of rows
size_t Matrix::Rows() const { return rows; }


/// returns number of columns
size_t Matrix::Cols() const { return cols; }




#ifndef NDEBUG
/// checks whether row and column numbers are compatible with the size of the matrix
bool Matrix::CheckRange( size_t m, size_t n,
                                const char* originator ) const
 {
    if ( m >= rows ) {
         std::cerr <<"\n"<< originator <<" row index violation, index="<< m;
         std::cerr <<" versus, row-max=" << rows << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckRange");
         return false;
      }
    if ( n >= cols ) {
         std::cerr <<"\n"<< originator <<" column index violation, index="<< n;
         std::cerr <<" versus, column-max=" << cols << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckRange");
         return false;
      }
    return true;
 }
#endif



/// checks whether the righthand matrix has the same size as the lefthand one in the calculation
bool Matrix::CheckSizes( const Matrix& mat, 
                                const char* originator ) const
 {
    if ( rows != mat.rows ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; rows1="<< rows;
         std::cerr <<" versus, rows2=" << mat.rows << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckSizes");
         return false;
      }
    if ( cols != mat.cols ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; columns1="<< cols;
         std::cerr <<" versus, columns2=" << mat.cols << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckSizes");
         return false;
      }
    return true;
 }


// AssignToDiagonal ( VectorVariable 2D)
// ---------------------------------------
void Matrix::AssignToDiagonal( const VectorVariable<2U>& vc )
 {
    assert ( rows == cols && cols == 2U );
    data[0][0] = vc(0);
    data[0][1] = static_cast<double64>(0.);
    data[1][0] = static_cast<double64>(0.);
    data[1][1] = vc(1);
 }



// AssignToDiagonal ( TensorVariable 3D)
// ---------------------------------------
void Matrix::AssignToDiagonal( const VectorVariable<3U>& vc )
 {
    assert ( rows == cols && cols == 3U );
    data[0][0] = vc(0);
    data[0][1] = static_cast<double64>(0.);
    data[0][2] = static_cast<double64>(0.);
    data[1][0] = static_cast<double64>(0.);
    data[1][1] = vc(1);
    data[1][2] = static_cast<double64>(0.);
    data[2][0] = static_cast<double64>(0.);
    data[2][1] = static_cast<double64>(0.);
    data[2][2] = vc(2); 
 }



// assign matrix to tensor variable
// --------------------------------
void Matrix::ExportTo( TensorVariable<2U>& ts ) const
 {
    if ( rows != 2U || rows != cols ) {
         std::cout <<"\nMatrix<>::ExportTo (tensor 2D): This matrix has wrong size.\n";
         return;
      } 
    ts(0,0) = data[0][0];
    ts(0,1) = data[0][1];
    ts(1,0) = data[1][0];
    ts(1,1) = data[1][1];
 }



void Matrix::ExportTo( TensorVariable<3U>& ts ) const
 {
    if ( rows != 3U || rows != cols ) {
         std::cout <<"\nMatrix<>::ExportTo (tensor 3D): This matrix has wrong size.\n";
         return;
      } 
    ts(0,0) = data[0][0];
    ts(0,1) = data[0][1];
    ts(0,2) = data[0][2];
    ts(1,0) = data[1][0];
    ts(1,1) = data[1][1];
    ts(1,2) = data[1][2];
    ts(2,0) = data[2][0];
    ts(2,1) = data[2][1];
    ts(2,2) = data[2][2];
 }

void Matrix::LUDecomposition(std::vector<size_t>& index, double64& d)
{
    const double64 TINY(1.0e-20);
    int i,imax(0),j,k;
    double64 big,dum,sum,temp;

    std::vector<double64> vv(rows);
    index.resize(rows);

    d=1.0;
    for (i=0;i<rows;i++) {
            big=0.0;
            for (j=0;j<rows;j++)
                if ((temp=std::fabs(data[i][j])) > big) big=temp;
            if (big == 0.0) {
                std::cout << "\nMatrix::LUDecomposition: Singular matrix" << std::endl;
                throw std::logic_error("Matrix::LUDecomposition: Singular matrix");
            }
            vv[i]=1.0/big;
    }
    for (j=0;j<rows;j++) {
            for (i=0;i<j;i++) {
                    sum=data[i][j];
                    for (k=0;k<i;k++) sum -= data[i][k]*data[k][j];
                    data[i][j]=sum;
            }
            big=0.0;
            for (i=j;i<rows;i++) {
                    sum=data[i][j];
                    for (k=0;k<j;k++) sum -= data[i][k]*data[k][j];
                    data[i][j]=sum;
                    if ((dum=vv[i]*std::fabs(sum)) >= big) {
                            big=dum;
                            imax=i;
                    }
            }
            if (j != imax) {
                    for (k=0;k<rows;k++) {
                            dum=data[imax][k];
                            data[imax][k]=data[j][k];
                            data[j][k]=dum;
                    }
                    d = -d;
                    vv[imax]=vv[j];
            }
            index[j]=imax;
            if (data[j][j] == 0.0) data[j][j]=TINY;
            if (j != rows-1) {
                    dum=1.0/(data[j][j]);
                    for (i=j+1;i<rows;i++) data[i][j] *= dum;
            }
    }


}

void Matrix::LUBackSubstitution(std::vector<size_t>& index, std::vector<double64> &b)
{

   if (index.size() != rows or b.size() != rows) {
        std::cout << "\nMatrix::LUBackSubstitution: Input vectors must be of length M.Rows" << std::endl;
        throw std::logic_error("Matrix::LUBackSubstitution: Input vectors must be of length M.Rows");
    }

  long i,ii=0,ip,j;
  double64 sum;
  for (i=0;i<rows;i++) {
      ip=index[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii != 0) for (j=ii-1;j<i;j++) sum -= data[i][j]*b[j];
      else if (sum != 0.0) ii=i+1;
      b[i]=sum;
    }
  for (i=rows-1;i>=0;i--) {
      sum=b[i];
      for (j=i+1;j<rows;j++) sum -= data[i][j]*b[j];
      b[i]=sum/data[i][i];
    }


}


// ReturnRow()
// ---------------------------------------------
std::vector<double64> Matrix::ReturnRow( size_t row ) const
 {
    std::vector<double64> row_vals(cols,0.);
    if ( row < rows ) {
        for ( size_t j=0U; j<cols; j++ ) row_vals[j] = data[row][j];
      }
    else {
       std::cout <<"\nMatrix<double>::ReturnRow: ";
       std::cout <<"Target rows does not exist: "<< row << std::endl;
      }
    return row_vals;
 }

// ReturnCol()
// ---------------------------------------------
std::vector<double64>  Matrix::ReturnCol( size_t col ) const
 {
    std::vector<double64> col_vals(rows,0.);
    if ( col < cols ) {
        for ( size_t i=0U; i<rows; i++ ) col_vals[i] = data[i][col];
      }
    else {
        std::cout <<"\nMatrix<double>::ReturnCol: ";
        std::cout <<"Target column does not exist: "<< col << std::endl;
      }
    return col_vals;
 }



void Matrix::AssignToRow( size_t row, std::vector<double64>& vec )
 {
    if ( vec.size() != cols ) {
      std::cout <<"\nMatrix<double>::AssignToRow: ";
      std::cout <<"Vector must have same length as number of cols: " << std::endl;
      return;
     }
    if ( row < rows ) for ( size_t i=0U; i<cols; i++ ) data[row][i] = vec[i];
    else {
      std::cout <<"\nMatrix<double>::AssignToRow: ";
      std::cout <<"Target column does not exist: "<< row << std::endl;
     }
 }



void Matrix::AssignToCol( size_t col, std::vector<double64>& vec )
{
  if ( vec.size() != rows ) {
    std::cout <<"\nMatrix<double>::AssignToCol: ";
    std::cout <<"Vector must have same length as number of rows: " << std::endl;
    return;
  }
  if ( col < cols ) for ( size_t i=0U; i<rows; i++ ) data[i][col] = vec[i];
  else {
    std::cout <<"\nMatrix<double>::AssignToCol: ";
    std::cout <<"Target column does not exist: "<< col << std::endl;
  }
}

/**
    @todo SKM whoever wrote this code should clarify its purpose
*/
Matrix Matrix::Minor(const size_t row, const size_t col)
{
  Matrix res;
  if (row <= rows-1 && col <= cols-1)
  {
    res = Matrix(rows - 1, cols - 1); // minor matrix
    // copy the content of the matrix to the minor, except the selected
    for (int r = 0; r <= (rows - (row >= (rows-1)))-1; r++)
    {
      for (int c = 0; c <= (cols - (col >= (cols-1)))-1; c++)
      {
        //printf("r=%i, c=%i, value=%f, rr=%i, cc=%i \n", r, c, p[r-1][c-1], r - (r > row), c - (c > col));
        res(r - (r > row), c - (c > col)) = data[r][c];
      }
    }
  }
  else
  {
    cout<<"Index out of range"<<endl;
  }
  return res;
}

} // csmp



