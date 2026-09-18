// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Matrix.h"

using namespace std;

namespace csmp {



Matrix& Matrix::operator=( double val )
 {
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) data(i,j) = val;
      
    return *this;
 }






// constructor (i,j, value)
// ---------------------------------------

Matrix::Matrix( size_t m, size_t n, double val )
 : data(m,n,val)
 {
 }




Matrix&  Matrix::operator*=( double val )
 {
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) data(i,j) *= val;
       
    return *this;
 }







// operator*=  matrix multiplication
// ---------------------------------------
/// @warning creates a temporary

Matrix&  Matrix::operator*=( const Matrix& mat )
 {
    if ( data.cols() != mat.data.rows() ) {
         std::cerr <<"\nMatrix<double>::operator*=: Matrices cannot be multiplied ";
         std::cerr <<"because of incompatible sizes (A *= B, see matrices below): "<< endl;
         Out();
         mat.Out();
         throw length_error("Matrix<mn_max>::operator*=");
      }

    //cout << "cols vs. rows."<<cols<<" "<<rows<<endl;
    // C++ named return value optimization
    //Matrix  temp( rows, mat.cols );
    Matrix  temp( data.cols(), mat.data.rows(),0.0);

    for ( size_t i{0U}; i<data.cols(); i++ )
    {
      for ( size_t j{0U}; j<mat.data.rows(); j++ ) {
        //temp.data(i,j) = static_cast<double>(0.0);
        for ( size_t k{0U}; k<mat.data.rows(); k++ )
        {
          temp.data(i,j) += this->data(i,k) * mat.data(k,j);
        }
      }
    }

    *this = temp;  
    //(*this).Out(3);
    return *this;
 }




// RowCondenseTo
// ---------------------------------------
/// vec = Mat * unity vector

void Matrix::RowCondenseTo( std::vector<double>& vec ) const
 {
    vec.resize(data.rows());
    vector<double>( vec ).swap( vec );
    for ( size_t i{0U}; i<data.rows(); i++ ) vec[i] = RowSum(i);
 }



// operator+=
// ---------------------------------------

Matrix& Matrix::operator+=( const Matrix& mat )
 {
    CheckSizes( mat, "Matrix::operator+=" );
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) data(i,j) += mat.data(i,j);
      
    return *this;
 }
 



// operator-=
// ---------------------------------------

Matrix& Matrix::operator-=( const Matrix& mat )
 {
    CheckSizes( mat, "Matrix::operator-=" );
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) data(i,j) -= mat.data(i,j);
      
    return *this;
 }





// operator+  
// ---------------------------------------

Matrix  Matrix::operator+( const Matrix& mat ) const
 {
    CheckSizes( mat, "Matrix::operator+" );
    Matrix  temp(data.cols(),data.rows());
    
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ )
        temp.data(i,j) = data(i,j) + mat.data(i,j);
      
    return temp;
 }

 
// operator-  
// ---------------------------------------

Matrix Matrix::operator-( const Matrix& mat ) const
 {
    CheckSizes( mat, "Matrix::operator-" );
    Matrix  temp(data.cols(),data.rows());
    
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ )
        temp.data(i,j) = data(i,j) - mat.data(i,j);
      
    return temp;
 }












// operator*  matrix multiplication
// ---------------------------------------

Matrix  Matrix::operator*( const Matrix& mat ) const
 {
//    cout << "cols vs. rows."<<cols<<" "<<rows<<endl;
    if ( data.cols() != mat.data.rows() ) {
         std::cout <<"\nMatrix<double>::operator*: Matrices cannot be multiplied ";
         std::cout <<"because of incompatible sizes (A *= B, see matrices below): "<< std::endl;
         Out(3);
         mat.Out(3);
         throw std::length_error("Matrix<mn_max>::operator*");
      }

    // use STANDARD return value optimization
    Matrix  temp( data.rows(), mat.data.cols() );
    //temp.Out(3);
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<mat.data.cols(); j++ ) {
           temp.data(i,j) = static_cast<double>(0.0);
           for ( size_t k{0U}; k<mat.data.rows(); k++ )
             temp.data(i,j) += data(i,k) * mat.data(k,j);
        }
    return temp;
 }




// operator*  matrix with STL vector<double> multiplication
// ------------------------------------------------

Matrix& Matrix::operator*=( const std::vector<double>& vecT )
 {
    if ( data.cols() != vecT.size() ) {
         std::cout <<"\nMatrix<double>::operator*=: ";
         std::cout <<"Matrix and vector<double> cannot be multiplied "; 
         std::cout <<"because of incompatible sizes (A(cols != vecT(size)): "<< vecT.size() << std::endl;
         throw std::length_error("Matrix<mn_max>::operator*=");
      }
    
    size_t  j;
    double  sum;
    
    for ( size_t i=0U; i<data.rows(); i++ ) {
         for ( sum=static_cast<double>(0.0), j=0ul; j<vecT.size(); j++ ) sum += data(i,j) * vecT[j];
         data(i,0U) = sum;
      }
      
    data.resize( data.rows(), 1 );
      
    return *this;
 }





/// RES = A B^T
// ---------------------------------------

void 
Matrix::MultiplyWithTransposedOf( const Matrix& B, Matrix& RES ) const
 {
    if ( data.cols() != B.data.cols() ) {
         cout <<"\nMatrix<double>::MultiplyWithTransposedOf: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out();
         B.Out();
         throw length_error("Matrix<mn_max>::MultiplyWithTransposedOf");
      }
    
    RES.Resize( data.rows(), B.data.rows() );
    double  sum;
    
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<B.data.rows(); j++ ) {
           sum = static_cast<double>(0.0);
           for ( size_t k{0U}; k<data.cols(); k++ )
             sum += data(i,k) * B.data(j,k);
           RES.data(i,j) = sum;
        }
 }




/// RES = A^T B
// ---------------------------------------

void 
Matrix::MultiplyTransposedOfWith( const Matrix& B, Matrix& RES ) const
 {
    if ( data.rows() != B.data.rows() ) {
         cout <<"\nMatrix<double>::MultiplyTransposedOfWith: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out();
         B.Out();
         throw length_error("Matrix<mn_max>::MultiplyTransposedOfWith");
      }
    
    RES.Resize( data.cols(), B.data.cols() );
    double  sum;
    
    for ( size_t i{0U}; i<data.cols(); i++ )
      for ( size_t j{0U}; j<B.data.cols(); j++ ) {
           sum = static_cast<double>(0.0);
           for ( size_t k{0U}; k<B.data.rows(); k++ )
             sum += data(k,i) * B.data(k,j);
           RES.data(i,j) = sum;
        }
 }




/// RES = A^T A
// ---------------------------------------

void 
Matrix::TransposedProduct( Matrix& RES ) const
 {
    RES.Resize( data.cols(), data.cols() );
    double  sum;
    
    for ( size_t i{0U}; i<data.cols(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) {
           sum = static_cast<double>(0.0);
           for ( size_t k{0U}; k<data.rows(); k++ )
             sum += data(k,i) * data(k,j);
           RES.data(j,i) = sum;
        }
 }




/// unscaled L_infinity matrix norm ||A||inf = max_i * sum_j |aij|

double Matrix::NormL_Infinity() const
 {
	if ( data.rows() == 0 || data.cols() == 0 ) {
         cerr <<"\nMatrix<double>::NormL_Infinity: empty matrix!\n";
         return std::numeric_limits<double>::signaling_NaN();
      }

	double  maxval(static_cast<double>(0.0)), sum;

	for ( size_t i=0U; i<data.rows(); i++ ) {
		 sum = static_cast<double>(0.0);
		 for ( size_t j{0U}; j<data.cols(); j++ )
		   sum += fabs(data(i,j));
		 maxval = max(maxval,sum);
	  }

	return maxval;
 }


/// unscaled L1 matrix norm ||A1|| = max_j * sum_i |aij|

double Matrix::NormL1() const
 {
	if ( data.rows() == 0 || data.cols() == 0 ) {
         cerr <<"\nMatrix<double>::NormL1: empty matrix!\n";
         return std::numeric_limits<double>::signaling_NaN();
      }

	double  maxval(static_cast<double>(0.0)), sum;

	for ( size_t j{0U}; j<data.cols(); j++ ) {
		 sum = static_cast<double>(0.0);
		 for ( size_t i{0U}; i<data.rows(); i++ )
		   sum += fabs(data(i,j));
		 maxval = max(maxval,sum);
	  }

	return maxval;
 }


// -------------------------------------------------------------------------------------------
// OPERATOR FUNCTIONS
// -------------------------------------------------------------------------------------------

/// vector^T = Matrix * vector^T
// ---------------------------------------

vector<double>  operator*( const Matrix& mat, const vector<double>& vec )
 {
    assert ( mat.Cols() == vec.size() );
    
    vector<double> temp(mat.Rows(),static_cast<double>(0.0));
    
    for ( size_t i=0U; i<mat.Rows(); i++ )
      for ( size_t j{0U}; j<mat.Cols(); j++ ) temp[i] += mat(i,j) * vec[j];
    
    return temp;
      
 } // end operator*



/// Matrix = vector^T * Matrix
// ---------------------------------------

Matrix  operator*( const vector<double>& vec, const Matrix& mat )
 {
    if ( vec.size() != mat.Rows() ) {
         cout <<"\noperator*: vector<double> cannot be multiplied with matrix"; 
         cout <<"because of incompatible sizes (v * M): "<< endl;
         throw length_error("Matrix<mn_max>::operator*");
      }
    Matrix  temp( vec.size(), mat.Cols() );
    
    for ( auto i{0}; i<vec.size(); i++ )
      for ( size_t j{0U}; j<mat.Cols(); j++ ) {
           temp(i,j) = static_cast<double>(0.0);
           for ( size_t k{0U}; k<mat.Rows(); k++ ) 
             temp(i,j) += vec[k] * mat(k,j);
        }
      
    return temp;
    
 } // end operator




// Resize()
// ---------------------------------------------

void Matrix::Resize( size_t m, size_t n )
 {
    data.resize( m, n );
 }


// Identity()
// ---------------------------------------------

void Matrix::Identity()
 {
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ )
        if ( i == j ) data(i,j) = static_cast<double>(1.0);
        else          data(i,j) = static_cast<double>(0.0);         
 }
 

// Zero()
// ---------------------------------------------

void Matrix::Zero()
 {
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ )
        data(i,j) = static_cast<double>(0.0);
 }
 

// ZeroRow()
// ---------------------------------------------

void Matrix::ZeroRow( size_t row )
 {
    if ( row < data.rows() ) for ( size_t j{0U}; j<data.cols(); j++ )
      data(row,j) = static_cast<double>(0.0);
      
    else {
        std::cout <<"\nMatrix<double>::ZeroRow: ";
        std::cout <<"Target rows does not exist: "<< row << std::endl;
        throw std::range_error("Matrix<double>::ZeroRow");
     }
 }


// ZeroCol()
// ---------------------------------------------

void Matrix::ZeroCol( size_t col )
 {
    if ( col < data.cols() ) for ( size_t i{0U}; i<data.rows(); i++ ) data(i,col) = static_cast<double>(0.0);
    else
    std::cout <<"\nMatrix<double>::ZeroCol: Target column does not exist: "<< col << std::endl;
 }



// Fill()
// ---------------------------------------------

void Matrix::Fill( double val )
 {
    for ( auto it=data.begin(); it!=data.end(); ++it ) (*it) = val;
 }


// FillRow()
// ---------------------------------------------

void Matrix::FillRow( size_t row, double val )
 {
    if ( row < data.rows() ) for ( size_t j{0U}; j<data.cols(); j++ ) data(row,j) = val;
    else {
      std::cout <<"\nMatrix<double>::FillRow: ";
      std::cout <<"Target rows does not exist: "<< row << std::endl;
     }
 }
 
 
// FillCol()
// ---------------------------------------------

void Matrix::FillCol( size_t col, double val )
 {
    if ( col < data.cols() ) for ( size_t i{0U}; i<data.rows(); i++ ) data(i,col) = val;
    else {
      std::cout <<"\nMatrix<double>::FillCol: ";
      std::cout <<"Target column does not exist: "<< col << std::endl;
     }
 }


// Transposed()
// ---------------------------------------------

void Matrix::Transposed( Matrix& M ) const
 {
    M.Resize(data.cols(),data.rows());
 
    for ( size_t i{0U}; i<data.rows(); i++ )
      for ( size_t j{0U}; j<data.cols(); j++ ) M.data(j,i) = data(i,j);
 }




// RowSum()
// ---------------------------------------------

double   Matrix::RowSum( size_t row ) const
 {
    double  sum(0.);
 
    if ( row < data.rows() ) {
         for ( size_t j{0U}; j<data.cols(); j++ ) sum += data(row,j);
         return sum;
      }

    std::cout <<"\nMatrix<double>::RowSum: ";
    std::cout <<"Target rows does not exist: "<< row << std::endl;
    throw std::length_error("Matrix<mn_max>::RowSum");

    return static_cast<double>(0.0);
 }
 

// ColSum()
// ---------------------------------------------

double   Matrix::ColSum( size_t col ) const
 {
    double  sum(0.0);
 
    if ( col < data.cols() ) {
         for ( size_t i{0ul}; i<data.rows(); i++ ) sum += data(i,col);
         return sum;
      }

    std::cout <<"\nMatrix<double>::ColSum: ";
    std::cout <<"Target column does not exist: "<< col << std::endl;
    throw std::length_error("Matrix<mn_max>::ColSum");

    return static_cast<double>(0.0);
 }



// In()
// ---------------------------------------

void Matrix::In()
 {
    size_t rows{0}, cols{0};
    cout <<"\nMatrix<double>: Please enter matrix entries as prompted for: "<< endl;
    cout <<"Enter number of rows and columns: ";
    cin >> rows >> cols;
    
    data.resize( rows, cols );
    
    for ( size_t i=0U; i<rows; i++ )
      {
         cout <<"\nEnter entries of row "<< i+1 <<": ";
         for ( size_t j{0U}; j<cols; j++ ) cin >> data(i,j);
      }

    cout <<"\nThank you."<< endl;    
 }



// Out( digits )
// ---------------------------------------

void Matrix::Out( long digits ) const
 {
    long  prec(cout.precision(digits));
    cout <<"\nMatrix<double>::Out(): m="<< data.rows() <<", n="<< data.cols() << endl;
    if ( digits != 0U ) cout.setf(ios::scientific);
    size_t row_break, split_adoubleer(10U);
     
    for ( size_t i=0U; i<data.rows(); i++ )
      {
         row_break = 1;
         for ( size_t j{0U}; j<data.cols(); j++, row_break++ )
           {
              if ( data(i,j) >= 0. ) cout <<" ";
              cout << data(i,j) <<" ";
              if ( row_break == split_adoubleer )
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





/// constructor (i,j)
Matrix::Matrix( size_t m, size_t n )
 : data(m,n)
 {
 }





/// operator (i,j)
double& Matrix::operator()( size_t m, size_t n )
 {
#ifndef NDEBUG
    CheckRange( m, n, "Matrix::operator()");
#endif
    return data(m,n);
 }



/// operator (i,j) const
const double& Matrix::operator()( size_t m, size_t n ) const
 {
#ifndef NDEBUG
    CheckRange( m, n, "Matrix::operator()");
#endif
    return data(m,n);
 }





/// returns number of rows
size_t Matrix::Rows() const { return data.rows(); }


/// returns number of columns
size_t Matrix::Cols() const { return data.cols(); }




#ifndef NDEBUG
/// checks whether row and column numbers are compatible with the size of the matrix
bool Matrix::CheckRange( size_t m, size_t n,
                         const char* originator ) const
 {
    if ( m >= data.rows() ) {
         std::cerr <<"\n"<< originator <<" row index violation, index="<< m;
         std::cerr <<" versus, row-max=" << data.rows() << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckRange");
         return false;
      }
    if ( n >= data.cols() ) {
         std::cerr <<"\n"<< originator <<" column index violation, index="<< n;
         std::cerr <<" versus, column-max=" << data.cols() << std::endl;
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
    if ( data.rows() != mat.data.rows() ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; rows1="<< data.rows();
         std::cerr <<" versus, rows2=" << mat.data.rows() << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckSizes");
         return false;
      }
    if ( data.cols() != mat.data.cols() ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; columns1="<< data.cols();
         std::cerr <<" versus, columns2=" << mat.data.cols() << std::endl;
         throw std::length_error("Matrix<mn_max>::CheckSizes");
         return false;
      }
    return true;
 }



void Matrix::LUDecomposition( std::vector<uint32_t>& index, double& d )
{
    const double TINY(1.0e-20);
    uint32_t i,imax(0),j,k;
    double big,dum,sum,temp;

    std::vector<double> vv(data.rows());
    index.resize(data.rows());

    d=1.0;
    for (i=0U;i<data.rows();i++) {
            big=0.0;
            for (j=0;j<data.rows();j++)
                if ((temp=std::fabs(data(i,j))) > big) big=temp;
            if (big == 0.0) {
                std::cout << "\nMatrix::LUDecomposition: Singular matrix" << std::endl;
                throw std::logic_error("Matrix::LUDecomposition: Singular matrix");
            }
            vv[i]=1.0/big;
    }
    for (j=0U;j<data.rows();j++) {
            for (i=0;i<j;i++) {
                    sum=data(i,j);
                    for (k=0;k<i;k++) sum -= data(i,k) * data(k,j);
                    data(i,j)=sum;
            }
            big=0.0;
            for (i=j;i<data.rows();i++) {
                    sum=data(i,j);
                    for (k=0;k<j;k++) sum -= data(i,k) * data(k,j);
                    data(i,j)=sum;
                    if ((dum=vv[i]*std::fabs(sum)) >= big) {
                            big=dum;
                            imax=i;
                    }
            }
            if (j != imax) {
                    for (k=0;k<data.rows();k++) {
                            dum=data(imax,k);
                            data(imax,k) = data(j,k);
                            data(j,k)=dum;
                    }
                    d = -d;
                    vv[imax]=vv[j];
            }
            index[j] = imax;
            if (data(j,j) == 0.0) data(j,j) = TINY;
            if (j != data.rows()-1) {
                    dum=1.0/(data(j,j));
                    for (i=j+1;i<data.rows();i++) data(i,j) *= dum;
            }
    }


}

void Matrix::LUBackSubstitution(std::vector<uint32_t>& index, std::vector<double> &b)
{

   if (index.size() != data.rows() or b.size() != data.rows()) {
        std::cout << "\nMatrix::LUBackSubstitution: Input vectors must be of length M.Rows" << std::endl;
        throw std::logic_error("Matrix::LUBackSubstitution: Input vectors must be of length M.Rows");
    }

  long i,ii=0,ip,j;
  double sum;
  for (i=0U;i<data.rows();i++) {
      ip=index[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii != 0) for (j=ii-1;j<i;j++) sum -= data(i,j)*b[j];
      else if (sum != 0.0) ii=i+1;
      b[i]=sum;
    }
  for (i=data.rows()-1;i>=0;i--) {
      sum=b[i];
      for (j=i+1;j<data.rows();j++) sum -= data(i,j)*b[j];
      b[i]=sum/data(i,i);
    }


}


// ReturnRow()
// ---------------------------------------------
std::vector<double> Matrix::ReturnRow( size_t row ) const
 {
    std::vector<double> row_vals(data.cols(),0.);
    if ( row < data.rows() ) {
        for ( size_t j{0U}; j<data.cols(); j++ ) row_vals[j] = data(row,j);
      }
    else {
       std::cout <<"\nMatrix<double>::ReturnRow: ";
       std::cout <<"Target rows does not exist: "<< row << std::endl;
      }
    return row_vals;
 }

// ReturnCol()
// ---------------------------------------------
std::vector<double>  Matrix::ReturnCol( size_t col ) const
 {
    std::vector<double> col_vals(data.rows(),0.);
    if ( col < data.cols() ) {
        for ( size_t i{0U}; i<data.rows(); i++ ) col_vals[i] = data(i,col);
      }
    else {
        std::cout <<"\nMatrix<double>::ReturnCol: ";
        std::cout <<"Target column does not exist: "<< col << std::endl;
      }
    return col_vals;
 }



void Matrix::AssignToRow( size_t row, std::vector<double>& vec )
 {
    if ( vec.size() != data.cols() ) {
      std::cout <<"\nMatrix<double>::AssignToRow: ";
      std::cout <<"Vector must have same length as number of cols: " << std::endl;
      return;
     }
    if ( row < data.rows() ) for ( size_t i{0U}; i<data.cols(); i++ ) data(row,i) = vec[i];
    else {
      std::cout <<"\nMatrix<double>::AssignToRow: ";
      std::cout <<"Target column does not exist: "<< row << std::endl;
     }
 }



void Matrix::AssignToCol( size_t col, std::vector<double>& vec )
{
  if ( vec.size() != data.rows() ) {
    std::cout <<"\nMatrix<double>::AssignToCol: ";
    std::cout <<"Vector must have same length as number of rows: " << std::endl;
    return;
  }
  if ( col < data.cols() ) for ( size_t i{0U}; i<data.rows(); i++ ) data(i,col) = vec[i];
  else {
    std::cout <<"\nMatrix<double>::AssignToCol: ";
    std::cout <<"Target column does not exist: "<< col << std::endl;
  }
}

/**
    Returns a submatrix of the matrix with the dimensions row x col
*/
Matrix Matrix::Minor(const size_t row, const size_t col)
{
   if ( row >= data.rows() || col >= data.cols() )
      throw out_of_range("Matrix::Minor");

    Matrix res = Matrix(data.rows() - 1, data.cols() - 1); // minor matrix
    // copy the content of the matrix to the minor, except the selected
    for ( long r = 0; r <= (data.rows() - (row >= (data.rows()-1)))-1; r++)
      {
        for ( long c = 0; c <= (data.cols() - (col >= (data.cols()-1)))-1; c++)
          {
            //printf("r=%i, c=%i, value=%f, rr=%i, cc=%i \n", r, c, p[r-1][c-1], r - (r > row), c - (c > col));
            res(r - (r > row), c - (c > col)) = data(r,c);
          }
      }

  return res;
}

} // csmp



