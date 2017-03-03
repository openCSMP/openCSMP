#include "DenseMatrix.h"
#include <algorithm>

// #define CSMP_DENSE_MATRIX_DEBUG // uncomment this to invoke debugging

using namespace std;

namespace csmp {

/// default constructor
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>::DenseMatrix()
 : rows(mn_max), cols(mn_max)
 {
//    cout <<"\nDenseMatrix<"<< mn_max <<","<< mn_max <<">: called constructor."<< endl;
 }


/// copy constructor
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( const DenseMatrix<mn_max>& mat )
 : data(mat.data),
   rows(mat.rows),
   cols(mat.cols)
 {
 }


<<<<<<< HEAD

=======
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
// constructor (i,j, value)
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( size_t m, size_t n, double64 val )
 : rows(m), cols(n)
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] = val;
 }


// destructor
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>::~DenseMatrix()
 {
 }


// operator=( DM )
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator=( const DenseMatrix<mn_max>& mat )
 {
    if ( &mat != this ) {
         rows  = mat.rows;
         cols  = mat.cols;
         data  = mat.data;
      }
    return *this;
 }



// Identity()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::Identity()
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        if ( i == j ) data[i][j] = static_cast<double64>(1.0);
        else          data[i][j] = static_cast<double64>(0.0);         
 }
 
 
/// constructor (i,j)
template<size_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( size_t m, size_t n )
 : rows(m), cols(n)
 {
 }



/// operator (i,j)
template<size_t mn_max>
double64& DenseMatrix<mn_max>::operator()( size_t m, size_t n )
 {
    CheckRange( m, n, "DenseMatrix<mn_max>::operator()");
    return data[m][n];
 }



/// operator (i,j) const
template<size_t mn_max>
const double64& DenseMatrix<mn_max>::operator()( size_t m, size_t n ) const
 {
    CheckRange( m, n, "DenseMatrix<mn_max>::operator()");
    return data[m][n];
 }


 
// Rows()
template<size_t mn_max>
size_t DenseMatrix<mn_max>::Rows() const { return rows; }


// Cols()
template<size_t mn_max>
size_t DenseMatrix<mn_max>::Cols() const { return cols; }


/**  
    Resizes DenseMatrix without allocation of new memory.
    If the capacitiy is exceeded an exception is thrown,
    but only if the code is compiled in debug mode.
*/
template<size_t mn_max>
void DenseMatrix<mn_max>::Resize( size_t m, size_t n )
 {
#ifndef NDEBUG 
    if ( m > mn_max ) {
         std::cerr <<"\nDenseMatrix<"<< mn_max;
         std::cerr <<">::Resize: Requested m-rows exceed matric capacity (";
         std::cerr << m <<" versus "<< rows <<")."<< std::endl;
         throw std::length_error("DenseMatrix<mn_max>::Resize");
      }
    if ( n > mn_max ) {
         std::cerr <<"\nDenseMatrix<"<<  mn_max;
         std::cerr <<">::Resize: Requested n-columns exceed matric capacity (";
         std::cerr << n <<" versus "<< cols <<")."<< std::endl;
         throw std::length_error("DenseMatrix<mn_max>::Resize");
      }
#endif
    rows = m;
    cols = n;
 }


#ifndef NDEBUG 
/// indices checking but only in the debug version
template<size_t mn_max>
bool DenseMatrix<mn_max>::CheckRange( size_t m, size_t n, 
                                             const char* originator ) const
 {
    if ( m >= rows ) {
         std::cerr <<"\n"<< originator <<" row index violation, index="<< m;
         std::cerr <<" versus, row-max=" << rows << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckRange");
         return false;
      }
    if ( n >= cols ) {
         std::cerr <<"\n"<< originator <<" column index violation, index="<< n;
         std::cerr <<" versus, column-max=" << cols << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckRange");
         return false;
      }
    return true;
 }
#else
template<size_t mn_max>
bool DenseMatrix<mn_max>::CheckRange( size_t, size_t, 
                                             const char* ) const
 {
    return true;
 }
#endif




#ifndef NDEBUG 
template<size_t mn_max>
/// checks (in DEBUG mode) whether the sizes of the matrices on either side of the expression match
bool DenseMatrix<mn_max>::CheckSizes( const DenseMatrix& mat, 
                                                const char* originator ) const
#else
template<size_t mn_max>
bool DenseMatrix<mn_max>::CheckSizes( const DenseMatrix&, 
                                                const char* ) const
#endif
 {
#ifndef NDEBUG 
    if ( rows != mat.rows ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; rows1="<< rows;
         std::cerr <<" versus, rows2=" << mat.rows << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckSizes");
         return false;
      }
    if ( cols != mat.cols ) {
         std::cerr <<"\n"<< originator <<" matrices have different sizes; columns1="<< cols;
         std::cerr <<" versus, columns2=" << mat.cols << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::CheckSizes");
         return false;
      }
#endif
    return true;
 }
 
 
 

// Zero()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::Zero()
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        data[i][j] = static_cast<double64>(0.0);
 }
 


// ZeroRow()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::ZeroRow( size_t row )
 {
    if ( row < rows ) for ( size_t j=0; j<cols; j++ ) 
      data[row][j] = static_cast<double64>(0.0);
      
    else {
      std::cout <<"\nDenseMatrix<"<< mn_max <<">::ZeroRow: ";
      std::cout <<"Target rows does not exist: "<< row << std::endl;
    }
 }


// ZeroCol()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::ZeroCol( size_t col )
 {
    if ( col < cols ) for ( size_t i=0; i<rows; i++ ) data[i][col] = static_cast<double64>(0.0);
    else
    std::cout <<"\nDenseMatrix<"<< mn_max <<">::ZeroCol: Target column does not exist: "<< col << std::endl;
 }


// Fill()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::Fill( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] = val;
 }


// FillRow()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::FillRow( size_t row, double64 val )
 {
    if ( row < rows ) for ( size_t j=0U; j<cols; j++ ) data[row][j] = val;
    else {
      std::cerr <<"\nDenseMatrix<"<< mn_max <<">::FillRow: ";
      std::cerr <<"Target rows does not exist: "<< row << std::endl;
    }
 }
 
 
// FillCol()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::FillCol( size_t col, double64 val )
 {
    if ( col < cols ) for ( size_t i=0U; i<rows; i++ ) data[i][col] = val;
    else {
      std::cerr <<"\nDenseMatrix<"<< mn_max <<">::FillCol: ";
      std::cerr <<"Target column does not exist: "<< col << std::endl;
     }
 }


// Transposed()
// ---------------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::Transposed( DenseMatrix& M ) const
 {
    M.Resize(cols,rows);
 
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) M.data[j][i] = data[i][j];
 }




// RowSum()
// ---------------------------------------------
template<size_t mn_max>
double64   DenseMatrix<mn_max>::RowSum( size_t row ) const
 {
    double64  sum(0.0);
 
    if ( row < rows ) {
         for ( size_t j=0U; j<cols; j++ ) sum += data[row][j];
         return sum;
      }

    std::cerr <<"\nDenseMatrix<"<< mn_max <<">::RowSum: ";
    std::cerr <<"Target rows does not exist: "<< row << std::endl;
    throw std::length_error("DenseMatrix<mn_max>::RowSum");

    return static_cast<double64>(0.0);
 }



// ColSum()
// ---------------------------------------------
template<size_t mn_max>
double64   DenseMatrix<mn_max>::ColSum( size_t col ) const
 {
    double64  sum(0.0);
 
    if ( col < cols ) {
         for ( size_t i=0U; i<rows; i++ ) sum += data[i][col];
         return sum;
      }

    std::cerr <<"\nDenseMatrix<"<< mn_max <<">::ColSum: ";
    std::cerr <<"Target column does not exist: "<< col << std::endl;
    throw std::length_error("DenseMatrix<mn_max>::ColSum");

    return static_cast<double64>(0.0);
 }


// CheckRange(i,j, message )









// POINT

#ifdef USED_TOGETHER_WITH_CSMP

// AssignRow
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const Point<1U>& pt )
 {
    cols = 1U;
    data[i][0] = pt[0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const Point<2U>& pt )
 {
    cols = 2U;
    data[i][0] = pt[0];
    data[i][1] = pt[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const Point<3U>& pt )
 {
    cols = 3U;
    data[i][0] = pt[0];
    data[i][1] = pt[1];
    data[i][2] = pt[2];
 }

// AssignCol
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const Point<1U>& pt )
 {
    rows = 1U;
    data[0][j] = pt[0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const Point<2U>& pt )
 {
    rows = 2U;
    data[0][j] = pt[0];
    data[1][j] = pt[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const Point<3U>& pt )
 {
    rows = 3U;
    data[0][j] = pt[0];
    data[1][j] = pt[1];
    data[2][j] = pt[2];
 }

// AssignToDiagonal
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const Point<1U>& pt )
 {
    Resize( 1U, 1U );
    data[0][0] = pt[0];
 }


template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const Point<2U>& pt )
 {
    Resize( 2U, 2U );
    data[0][0] = pt[0];
    data[0][1] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = pt[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const Point<3U>& pt )
 {
    Resize( 3U, 3U );
    data[0][0] = pt[0];
    data[0][1] = static_cast<double64>(0.0);
    data[0][2] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = pt[1];
    data[1][2] = static_cast<double64>(0.0);
    data[2][0] = static_cast<double64>(0.0);
    data[2][1] = static_cast<double64>(0.0);
    data[2][2] = pt[2];
 }

/// 1D version
template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const Point<1U>& pt )
 {
    assert ( cols == 1U );

    // when loop unrolling is not possible
    if ( rows != cols ) {
         for ( size_t i=0; i<rows; i++ )
            data[i][0] *= pt[0];
         cols = 1U;
         return *this;
      }

    data[0][0] *= pt[0];
    cols = 1U;

    return *this;
 }


/// 2D version
template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const Point<2U>& pt )
 {
    assert ( cols == 2U );
    double64 sum;

    // when loop unrolling is not possible
    if ( rows != cols ) {
        for ( size_t i=0; i<rows; i++ ) {
               sum = static_cast<double64>(0.0);
               for ( size_t j=0; j<cols; j++ )
                 sum +=  data[i][j] * pt[j];
               data[i][0] = sum;
            }
         cols = 1U;
         return *this;
      }

    sum = data[0][0] * pt[0] + data[0][1] * pt[1];
    data[0][0] = sum;
    sum = data[1][0] * pt[0] + data[1][1] * pt[1];
    data[1][0] = sum;

    cols = 1U;

    return *this;
 }


/// 3D version
template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const Point<3U>& pt )
 {
    assert ( cols == 3U );
    double64 sum;

    // when loop unrolling is not possible
    if ( rows != cols ) {
        for ( size_t i=0; i<rows; i++ ) {
               sum = static_cast<double64>(0.0);
               for ( size_t j=0; j<cols; j++ )
                 sum +=  data[i][j] * pt[j];
               data[i][0] = sum;
            }
         cols = 1U;
         return *this;
      }

    sum = data[0][0] * pt[0] + data[0][1] * pt[1] + data[0][2] * pt[2];
    data[0][0] = sum;
    sum = data[1][0] * pt[0] + data[1][1] * pt[1] + data[1][2] * pt[2];
    data[1][0] = sum;
    sum = data[2][0] * pt[0] + data[2][1] * pt[1] + data[2][2] * pt[2];
    data[2][0] = sum;

    cols = 1U;

    return *this;
 }





// SCALAR VARIABLE

// AssignToDiagonal ( ScalarVariable)
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( size_t diag_elmts,
                                            const ScalarVariable& sc )
 {
    Resize( diag_elmts, diag_elmts );
    for ( size_t i=0U; i<diag_elmts; i++ )
      for ( size_t j=0U; j<diag_elmts; j++ )
        data[i][j] = (i==j) ? sc() : static_cast<double64>(0.);
 }

template<size_t mn_max>
DenseMatrix<mn_max>&
  DenseMatrix<mn_max>::operator*=( const ScalarVariable& sc )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ )
          data[i][j] *= sc();
    return *this;
 }





// VECTOR VARIABLE

// AssignRow
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const VectorVariable<1U>& vc )
 {
    cols = 1U;
    data[i][0] = vc[0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const VectorVariable<2U>& vc )
 {
    cols = 2U;
    data[i][0] = vc[0];
    data[i][1] = vc[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const VectorVariable<3U>& vc )
 {
    cols = 3U;
    data[i][0] = vc[0];
    data[i][1] = vc[1];
    data[i][2] = vc[2];
 }

// AssignCol
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const VectorVariable<1U>& vc )
 {
    rows = 1U;
    data[0][j] = vc[0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const VectorVariable<2U>& vc )
 {
    rows = 2U;
    data[0][j] = vc[0];
    data[1][j] = vc[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const VectorVariable<3U>& vc )
 {
    rows = 3U;
    data[0][j] = vc[0];
    data[1][j] = vc[1];
    data[2][j] = vc[2];
 }

// AssignToDiagonal
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const VectorVariable<1U>& vc )
 {
    Resize( 1U, 1U );
    data[0][0] = vc[0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const VectorVariable<2U>& vc )
 {
    Resize( 2U, 2U );
    data[0][0] = vc[0];
    data[0][1] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = vc[1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const VectorVariable<3U>& vc )
 {
    Resize( 3U, 3U );
    data[0][0] = vc[0];
    data[0][1] = static_cast<double64>(0.0);
    data[0][2] = static_cast<double64>(0.0);
    data[1][0] = static_cast<double64>(0.0);
    data[1][1] = vc[1];
    data[1][2] = static_cast<double64>(0.0);
    data[2][0] = static_cast<double64>(0.0);
    data[2][1] = static_cast<double64>(0.0);
    data[2][2] = vc[2];
 }


template<size_t mn_max>
DenseMatrix<mn_max>&
  DenseMatrix<mn_max>::operator*=( const VectorVariable<1U>& vc )
 {
    assert ( cols == 1U );

    // when loop unrolling is not possible
    if ( rows != cols ) {
         for ( size_t i=0; i<rows; i++ )
            data[i][0] *= vc[0];
         cols = 1U;
         return *this;
      }

    data[0][0] *= vc[0];
    cols = 1U;

    return *this;
 }

template<size_t mn_max>
DenseMatrix<mn_max>&
  DenseMatrix<mn_max>::operator*=( const VectorVariable<2U>& vc )
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


template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const VectorVariable<3U>& vc )
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




// ARRAY VARIABLE

// AssignRow
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const ArrayVariable& ar )
 {
    cols = ar.Size();
    for ( size_t j=0U; j<cols; j++ )
        data[i][j] = ar[j];
 }

// AssignCol
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const ArrayVariable& ar )
 {
    rows = ar.Size();
    for ( size_t i=0U; i<rows; i++ )
        data[i][j] = ar[i];
 }

// AssignToDiagonal
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const ArrayVariable& ar )
 {
    Resize( ar.Size(), ar.Size() );
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ )
        data[i][j] = (i==j) ? ar[i] : static_cast<double64>(0.);
 }

template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const ArrayVariable& ar )
 {
    if ( cols != ar.Size() ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=(ArrayVariable): matrix.cols != array.size\n";
         throw std::range_error("DenseMatrix<>::operator*=(ArrayVariable): matrix.cols != array.size");
      }
    double64 sum;
    for ( size_t i=0; i<rows; i++ ) {
           sum = static_cast<double64>(0.0);
           for ( size_t j=0; j<cols; j++ )
             sum +=  data[i][j] * ar[j];
           data[i][0] = sum;
        }
    cols = 1U;
    return *this;
 }




// FLAGGED ARRAY VARIABLE

// AssignRow
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignRow( size_t i, const FlaggedArrayVariable& fr )
 {
    cols = fr.Size();
    for ( size_t j=0U; j<cols; j++ )
        data[i][j] = fr[j];
 }

// AssignCol
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignCol( size_t j, const FlaggedArrayVariable& fr )
 {
    rows = fr.Size();
    for ( size_t i=0U; i<rows; i++ )
        data[i][j] = fr[i];
 }

// AssignToDiagonal
// ---------------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const FlaggedArrayVariable& fr )
 {
    Resize( fr.Size(), fr.Size() );
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ )
        data[i][j] = (i==j) ? fr[i] : static_cast<double64>(0.);
 }

template<size_t mn_max>
DenseMatrix<mn_max>&
DenseMatrix<mn_max>::operator*=( const FlaggedArrayVariable& fr )
 {
    if ( cols != fr.Size() ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=(FlaggedArrayVariable): matrix.cols != array.size\n";
         throw std::range_error("DenseMatrix<>::operator*=(FlaggedArrayVariable): matrix.cols != array.size");
      }
    double64 sum;
    for ( size_t i=0; i<rows; i++ ) {
           sum = static_cast<double64>(0.0);
           for ( size_t j=0; j<cols; j++ )
             sum +=  data[i][j] * fr[j];
           data[i][0] = sum;
        }
    cols = 1U;
    return *this;
 }




// TENSOR VARIABLE

// operator= ( TensorVariable 1D )
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& 
  DenseMatrix<mn_max>::operator=( const TensorVariable<1U>& ts )
 {
    rows = cols = 1U;
    data[0][0] = ts(0,0);
      
    return *this; 
 }


// operator= ( TensorVariable 2D )
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& 
  DenseMatrix<mn_max>::operator=( const TensorVariable<2U>& ts )
 {
    rows = cols = 2U;
    data[0][0] = ts(0,0);
    data[0][1] = ts(0,1);
    data[1][0] = ts(1,0);
    data[1][1] = ts(1,1);
      
    return *this; 
 }



// operator= ( TensorVariable 3D )
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& 
  DenseMatrix<mn_max>::operator=( const TensorVariable<3U>& ts )
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

// assign matrix to tensor variable
// --------------------------------

template<size_t mn_max>
void DenseMatrix<mn_max>::ExportTo( TensorVariable<1U>& ts ) const
 {
    if ( rows != 1U || rows != cols ) {
         std::cerr <<"\nDenseMatrix<>::ExportTo (tensor 1D): This matrix has wrong size.\n";
         throw std::range_error("DenseMatrix<>::ExportTo(TensorVar.1D): matrix.cols != tensor.dim");
         return;
      }
    ts(0,0) = data[0][0];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::ExportTo( TensorVariable<2U>& ts ) const
 {
    if ( rows != 2U || rows != cols ) {
         std::cerr <<"\nDenseMatrix<>::ExportTo (tensor 2D): This matrix has wrong size.\n";
         throw std::range_error("DenseMatrix<>::ExportTo(TensorVar.2D): matrix.cols != tensor.dim");
         return;
      }
    ts(0,0) = data[0][0];
    ts(0,1) = data[0][1];
    ts(1,0) = data[1][0];
    ts(1,1) = data[1][1];
 }

template<size_t mn_max>
void DenseMatrix<mn_max>::ExportTo( TensorVariable<3U>& ts ) const
 {
    if ( rows != 3U || rows != cols ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max;
         std::cerr <<">::ExportTo (tensor 3D): This matrix has wrong size.\n";
         throw std::range_error("DenseMatrix<>::ExportTo(TensorVar.3D): matrix.cols != tensor.dim");
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

// matrix - tensor multiplication
// ------------------------------

template<size_t mn_max>
DenseMatrix<mn_max>&
  DenseMatrix<mn_max>::operator*=( const TensorVariable<1U>& ts )
 {
    if ( cols != 1U ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=(TensorVar.1D): matrix.cols != tensor.dim\n";
         throw std::range_error("DenseMatrix<>::operator*=(TensorVar.1D): matrix.cols != tensor.dim");
      }

    // when loop unrolling is not possible
    if ( rows != cols ) {
         for ( size_t i=0; i<rows; i++ )
            data[i][0] *= ts(0,0);
         cols = 1U;
         return *this;
      }

    data[0][0] *= ts(0,0);
    cols = 1U;

    return *this;
 }

template<size_t mn_max>
DenseMatrix<mn_max>& 
  DenseMatrix<mn_max>::operator*=( const TensorVariable<2U>& ts )
 {
    if ( cols != 2U ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=(TensorVar.2D): matrix.cols != tensor.dim\n";
         throw std::range_error("DenseMatrix<>::operator*=(TensorVar.2D): matrix.cols != tensor.dim");
      }
      
    // when loop unrolling is not possible
    if ( rows != cols ) {
         // C++ named return value optimization
         DenseMatrix<mn_max>  temp( rows, 2U );
    
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


template<size_t mn_max>
DenseMatrix<mn_max>&
  DenseMatrix<mn_max>::operator*=( const TensorVariable<3U>& ts )
 {
    if ( cols != 3U ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=(TensorVar.3D): matrix.cols != tensor.dim\n";
         throw(range_error("DenseMatrix<>::operator*=(TensorVar.3D): matrix.cols != tensor.dim"));
      }

    // when loop unrolling is not possible
    if ( rows != cols ) {
         // C++ named return value optimization
         DenseMatrix<mn_max>  temp( rows, 3U );

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

    return *this = std::move(temp);
 }

#endif // USED_TOGETHER WITH CSMP












// AssignToDiagonal (double64)
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( size_t diag_elmts,
                                            const double64& sc )
 {
    Resize( diag_elmts, diag_elmts );
    for ( size_t i=0U; i<diag_elmts; i++ )
      for ( size_t j=0U; j<diag_elmts; j++ )
        data[i][j] = (i==j) ? sc : static_cast<double64>(0.);
 }

// RowCondenseTo
// ---------------------------------------
/// vec = Mat * unity vector
template<size_t mn_max>
void DenseMatrix<mn_max>::RowCondenseTo( std::vector<double64>& vec ) const
 {
    vec.resize(rows);
    vector<double64>( vec ).swap( vec );
    for ( size_t i=0U; i<rows; i++ ) vec[i] = RowSum(i);
 } 



// operator+=
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator+=( const DenseMatrix<mn_max>& mat )
 {
    CheckSizes( mat, "DenseMatrix<mn_max>::operator+=" );

    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] += mat.data[i][j]; 
      
    return *this;
 }
 

// operator-=
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator-=( const DenseMatrix<mn_max>& mat )
 {
    CheckSizes( mat, "DenseMatrix<mn_max>::operator-=" );

    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] -= mat.data[i][j]; 
      
    return *this;
 }





// operator+  
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>  DenseMatrix<mn_max>::operator+( 
                                             const DenseMatrix<mn_max>& mat ) const
 {
    CheckSizes( mat, "DenseMatrix<mn_max>::operator+" );

    DenseMatrix<mn_max>  temp(cols,rows);
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        temp.data[i][j] = data[i][j] + mat.data[i][j]; 
      
    return temp;
 }

 
// operator-  
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max> DenseMatrix<mn_max>::operator-( 
                                            const DenseMatrix<mn_max>& mat ) const
 {
    CheckSizes( mat, "DenseMatrix<mn_max>::operator-" );

    DenseMatrix<mn_max>  temp(cols,rows);
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) 
        temp.data[i][j] = data[i][j] - mat.data[i][j]; 
      
    return temp;
 }




// operator*  matrix multiplication
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>  DenseMatrix<mn_max>::operator*(
                                             const DenseMatrix<mn_max>& mat ) const
 {
    if ( cols != mat.rows ) {
         cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*: Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B, see matrices below): "<< std::endl;
         Out(cout, 3);
         mat.Out(cout, 3);
         throw std::length_error("DenseMatrix<mn_max>::operator*");
      }

    // use STANDARD return value optimization
    DenseMatrix<mn_max>  temp( rows, mat.cols );
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<mat.cols; j++ ) {
           temp.data[i][j] = static_cast<double64>(0.0);
           for ( size_t k=0U; k<mat.rows; k++ ) 
             temp.data[i][j] += data[i][k]*mat.data[k][j];
        }  
    return temp;
 }




// operator*  matrix with STL vector multiplication
// ------------------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const std::vector<double64>& vecT )
 {
    if ( cols != vecT.size() ) {
         std::cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=: ";
         std::cout <<"Matrix and vector cannot be multiplied "; 
         std::cout <<"because of incompatible sizes (A(cols != vecT(size)): "<< vecT.size() << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::operator*=");
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
template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const double64* vecT ) 
 {
    if ( vecT == NULL ) {
         std::cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::operator*=: Matrix and vector cannot be multiplied "; 
         std::cout <<"because vecT C-array is not initialized." << std::endl;
         throw std::length_error("DenseMatrix<mn_max>::operator*=");
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





template<size_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator=( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] = val; 
      
    return *this;
 }






template<size_t mn_max>
DenseMatrix<mn_max>& 
  DenseMatrix<mn_max>::operator*=( double64 val )
 {
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<cols; j++ ) data[i][j] *= val;
       
    return *this;
 }



/// operator*=  matrix multiplication
// ---------------------------------------
// creates a temporary   
template<size_t mn_max>
DenseMatrix<mn_max>&  DenseMatrix<mn_max>::operator*=( const DenseMatrix<mn_max>& mat )
 {
#ifdef CSMP_DENSE_MATRIX_DEBUG
    if ( cols != mat.rows ) {
         cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max;
         cout <<">::operator*=: Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B, see matrices below): "<< endl;
         Out(cout);
         mat.Out(cout);
         throw length_error("DenseMatrix<mn_max>::operator*=");
      }
#endif
    
    // C++ named return value optimization
    DenseMatrix<mn_max>  temp( rows, mat.cols );
    
    for ( size_t i=0U; i<rows; i++ )
      for ( size_t j=0U; j<mat.cols; j++ ) {
           temp.data[i][j] = static_cast<double64>(0.0);
           for ( size_t k=0U; k<mat.rows; k++ ) 
             temp.data[i][j] += data[i][k] * mat.data[k][j];
        }
        
    *this = temp;  
      
    return *this;
 }



/// RES = A B^T
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::MultiplyWithTransposedOf( const DenseMatrix<mn_max>& B,
                                                    DenseMatrix<mn_max>& RES ) const
 {
    if ( cols != B.cols ) {
         cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max;
         cout <<">::MultiplyWithTransposedOf: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out(cout);
         B.Out(cout);
         throw length_error("DenseMatrix<mn_max>::MultiplyWithTransposedOf");
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
template<size_t mn_max>
void 
DenseMatrix<mn_max>::MultiplyTransposedOfWith( const DenseMatrix<mn_max>& B, 
                                                      DenseMatrix<mn_max>& RES ) const
 {
    if ( rows != B.rows ) {
         cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::MultiplyTransposedOfWith: ";
         cout <<"Matrices cannot be multiplied "; 
         cout <<"because of incompatible sizes (A *= B^T, see matrices below): "<< endl;
         Out(cout);
         B.Out(cout);
         throw length_error("DenseMatrix<mn_max>::MultiplyTransposedOfWith");
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
template<size_t mn_max>
void 
DenseMatrix<mn_max>::TransposedProduct( DenseMatrix<mn_max>& RES ) const
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
template<size_t mn_max>
double64 DenseMatrix<mn_max>::NormL_Infinity() const
 {
	if ( rows == 0 || cols == 0 ) {
         cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::NormL_Infinity: empty matrix!\n";
         return std::numeric_limits<double>::quiet_NaN();
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



/// operator +
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>  operator+( const DenseMatrix<mn_max>& a, 
                                   const DenseMatrix<mn_max>& b )
 {
    DenseMatrix<mn_max> temp(a); 
    temp += b;
    return temp; 
 }


/// operator -
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>  operator-( const DenseMatrix<mn_max>& a, 
                                   const DenseMatrix<mn_max>& b )
 {
    DenseMatrix<mn_max> temp(a);  
    temp -= b;
    return temp; 
 }


/// matrix multiplication operator* that returns temporary
template<size_t mn_max>
DenseMatrix<mn_max>  operator*( const DenseMatrix<mn_max>& a,
                                const DenseMatrix<mn_max>& b )
 {
#ifdef CSMP_DENSE_MATRIX_DEBUG
    if ( a.Cols() != b.Rows() ) {
         std::cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max;
         std::cerr <<"> operator*: Matrices cannot be multiplied ";
         std::cerr <<"because of incompatible sizes (A * B, see matrices below): "<< endl;
         a.Out(std::cerr, 3L);
         b.Out(std::cerr, 3L);
         throw length_error("DenseMatrix<mn_max>::operator*");
      }
#endif
    DenseMatrix<mn_max>  temp( a.Rows(), b.Cols() );
    
    for ( size_t i=0U; i<a.Rows(); i++ )
      for ( size_t j=0U; j<b.Cols(); j++ ) {
           temp(i,j) = static_cast<double64>(0.0);
           for ( size_t k=0U; k<b.Rows(); k++ ) 
             temp(i,j) += a(i,k) * b(k,j);
        }
      
    return temp;
 }




/// computes unscaled L1 matrix norm ||A1|| = max_j * sum_i |aij|
template<size_t mn_max>
double64 DenseMatrix<mn_max>::NormL1() const
 {
	if ( rows == 0 || cols == 0 ) {
         cerr <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">::NormL1: empty matrix!\n";
         return std::numeric_limits<double>::quiet_NaN();
      }

	double64  maxval(static_cast<double64>(0.0)), sum;

	for ( size_t j=0; j<cols; j++ ) {
		 sum = static_cast<double64>(0.0);
		 for ( size_t i=0; i<rows; i++ )
		   sum += std::fabs(data[i][j]);
		 maxval = std::max(maxval,sum);
	  }

	return maxval;
 }
 






// -------------------------------------------------------------------------------------------
// OPERATOR FUNCTIONS
// -------------------------------------------------------------------------------------------

/// vector^T = Matrix * vector^T
// ---------------------------------------
template<size_t mn_max>
std::vector<double64>  operator*( const DenseMatrix<mn_max>& mat, const std::vector<double64>& vec )
 {
    assert ( mat.Cols() == vec.size() );
    
    std::vector<double64> temp(mat.Rows(),static_cast<double64>(0.0));
    
    for ( size_t i=0; i<mat.Rows(); i++ )
      for ( size_t j=0; j<mat.Cols(); j++ ) temp[i] += mat(i,j) * vec[j];
    
    return temp;
      
 } // end operator*



/// Matrix = vector^T * Matrix
// ---------------------------------------
template<size_t mn_max>
DenseMatrix<mn_max>  operator*( const std::vector<double64>& vec, const DenseMatrix<mn_max>& mat )
 {
    if ( vec.size() != mat.Rows() ) {
         std::cout <<"\noperator*: vector cannot be multiplied with matrix"; 
         std::cout <<"because of incompatible sizes (v * M): "<< std::endl;
         throw std::length_error("DenseMatrix<mn_max>::operator*");
      }

    DenseMatrix<mn_max>  temp( vec.size(), mat.Cols() );
    
    for ( size_t i=0U; i<vec.size(); i++ )
      for ( size_t j=0U; j<mat.Cols(); j++ ) {
           temp(i,j) = static_cast<double64>(0.0);
           for ( size_t k=0U; k<mat.Rows(); k++ ) 
             temp(i,j) += vec[k] * mat(k,j);
        }
      
    return temp;
    
 } // end operator


// In()
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::In()
 {
    cout <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max;
    cout <<">: Please enter matrix entries as prompted for: "<< endl;
    cout <<"Enter number of rows and columns: ";
    cin >> rows >> cols;
    
    if ( rows >= mn_max || cols >= mn_max ) {
         cout <<"\nTarget matrix cannot accomodate this size; rebuild it please."<< endl;
         throw length_error("DenseMatrix<mn_max>In");
      }
    
    for ( size_t i=0; i<rows; i++ )
      {
         cout <<"\nEnter entries of row "<< i+1 <<": ";
         for ( size_t j=0; j<cols; j++ ) cin >> data[i][j];
      }

    cout <<"\nThank you."<< endl;    
 }



// Out( digits )
// ---------------------------------------
template<size_t mn_max>
void DenseMatrix<mn_max>::Out( std::ostream& os, long digits ) const
 {
    long  prec(2);
    os <<"\nDenseMatrix<"<< typeid(double64).name() <<","<< mn_max <<">: m=n="<< mn_max << endl;
    if ( digits != 0U ) { 
         os.setf(ios::scientific);
         prec = os.precision(digits);
      }
    size_t row_break, split_after(10U);
     
    for ( size_t i=0; i<rows; i++ )
      {
         row_break = 1;
         for ( size_t j=0; j<cols; j++, row_break++ )
           {
              if ( data[i][j] >= 0. ) os <<" ";
              os << data[i][j] <<" ";
              if ( row_break == split_after )
                {
                   os << endl;
                   row_break = 0U;
                }
           }
         os << endl;
      }
      
    if ( digits != 0U ) {
         os.unsetf( ios::scientific );
         os.precision(prec);
      }
      
    os << endl;    
    
 } // end Out()


// ------------------------------------------------------------------------------
// template instantiations
// ------------------------------------------------------------------------------

#ifdef CSMP_NUMBER_TYPES_H  // if matrix is used with CSMP


template class DenseMatrix<DM3>;

template DenseMatrix<DM3>  operator+( 
                      const DenseMatrix<DM3>& a, 
                      const DenseMatrix<DM3>& b );

template DenseMatrix<DM3>  operator-(
                      const DenseMatrix<DM3>& a, 
                      const DenseMatrix<DM3>& b );

template DenseMatrix<DM3>  operator*(
                      const DenseMatrix<DM3>& a, 
                      const DenseMatrix<DM3>& b );
template
vector<double64>  operator*( const DenseMatrix<DM3>& mat, 
                             const vector<double64>& vec );



                      
template class DenseMatrix<DM4>;

template DenseMatrix<DM4>  operator+( 
                      const DenseMatrix<DM4>& a, 
                      const DenseMatrix<DM4>& b );

template DenseMatrix<DM4>  operator-(
                      const DenseMatrix<DM4>& a, 
                      const DenseMatrix<DM4>& b );

template DenseMatrix<DM4>  operator*(
                      const DenseMatrix<DM4>& a, 
                      const DenseMatrix<DM4>& b );
template
vector<double64>  operator*( const DenseMatrix<DM4>& mat, 
                              const vector<double64>& vec );


                      
                      
template class DenseMatrix<DM12>;

template DenseMatrix<DM12>  operator+( 
                      const DenseMatrix<DM12>& a, 
                      const DenseMatrix<DM12>& b );

template DenseMatrix<DM12>  operator-(
                      const DenseMatrix<DM12>& a, 
                      const DenseMatrix<DM12>& b );

template DenseMatrix<DM12>  operator*(
                      const DenseMatrix<DM12>& a, 
                      const DenseMatrix<DM12>& b );
template
vector<double64>  operator*( const DenseMatrix<DM12>& mat, 
                              const vector<double64>& vec );




// big enough for TRI_3 (2D, 2DOF)
template class DenseMatrix<DM_MIN>;

template DenseMatrix<DM_MIN>  operator+( 
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );

template DenseMatrix<DM_MIN>  operator-(
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );

template DenseMatrix<DM_MIN>  operator*(
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );
template
vector<double64>  operator*( const DenseMatrix<DM_MIN>& mat, 
                              const vector<double64>& vec );




// big enough for HEXA_27 (3D, 3DOF) -> 52 kbyte
template class DenseMatrix<DM_MAX>;

template DenseMatrix<DM_MAX>  operator+( 
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );

template DenseMatrix<DM_MAX>  operator-(
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );

template DenseMatrix<DM_MAX>  operator*(
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );
template
vector<double64>  operator*( const DenseMatrix<DM_MAX>& mat, 
                              const vector<double64>& vec );






#else // if CSMP is not known

// big enough for TRI_3 (2D, 2DOF)
template class DenseMatrix<DM_MIN>;

template DenseMatrix<DM_MIN>  operator+( 
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );

template DenseMatrix<DM_MIN>  operator-(
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );

template DenseMatrix<DM_MIN>  operator*(
                      const DenseMatrix<DM_MIN>& a, 
                      const DenseMatrix<DM_MIN>& b );

// big enough for HEXA_27 (3D, 3DOF) -> 52 kbyte
template class DenseMatrix<DM_MAX>;

template DenseMatrix<DM_MAX>  operator+( 
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );

template DenseMatrix<DM_MAX>  operator-(
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );

template DenseMatrix<DM_MAX>  operator*(
                      const DenseMatrix<DM_MAX>& a, 
                      const DenseMatrix<DM_MAX>& b );

template
vector<double64>  operator*( const DenseMatrix<DM_MIN>& mat, 
                           const vector<double64>& vec );
template
vector<double64>  operator*( const DenseMatrix<DM_MAX>& mat, 
                           const vector<double64>& vec );
#endif


// extra operators

template DenseMatrix<DM_MAX>  
operator*( const vector<double64>& vec, const DenseMatrix<DM_MAX>& mat );

} // end namespace csmp
