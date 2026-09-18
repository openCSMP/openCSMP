// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "DenseMatrix.h"

#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#endif

using namespace std;

namespace csmp {

// =============================================================================
// Constructors
// =============================================================================

template<uint32_t mn_max>
DenseMatrix<mn_max>::DenseMatrix()
  : rows(mn_max), cols(mn_max)
 {}


template<uint32_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( const std::initializer_list<std::initializer_list<double>>& init )
 {
    if ( (rows = static_cast<uint32_t>(init.size())) > mn_max )
      throw out_of_range("DenseMatrix(ctor): initialiser list has more rows than storage capacity");
    const auto it = init.begin();
    if ( (cols = static_cast<uint32_t>((*it).size())) > mn_max )
      throw out_of_range("DenseMatrix(ctor): first row in initialiser list has more columns than storage capacity");

    uint32_t i{0U};
    for ( const auto& row : init )
      {
        if ( cols != row.size() )
          throw out_of_range("DenseMatrix(ctor): inconsistent row length in initialiser list");
        uint32_t j{0U};
        for ( const auto& val : row )
          data[i][j++] = val;
        ++i;
      }
 }


template<uint32_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( uint32_t m, uint32_t n )
  : rows(m), cols(n)
 {
    static_assert( mn_max > 0U, "DenseMatrix: mn_max must be > 0" );
 }


template<uint32_t mn_max>
DenseMatrix<mn_max>::DenseMatrix( uint32_t m, uint32_t n, double val )
  : rows(m), cols(n)
 {
    static_assert( mn_max > 0U, "DenseMatrix: mn_max must be > 0" );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = val;
 }


// =============================================================================
// Resize and range checking
// =============================================================================

template<uint32_t mn_max>
void DenseMatrix<mn_max>::Resize( uint32_t m, uint32_t n )
 {
#ifndef NDEBUG
    if ( m > mn_max )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::Resize: requested rows "<< m
             <<" exceed capacity "<< mn_max << endl;
        throw length_error("DenseMatrix::Resize: rows exceed capacity");
      }
    if ( n > mn_max )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::Resize: requested cols "<< n
             <<" exceed capacity "<< mn_max << endl;
        throw length_error("DenseMatrix::Resize: cols exceed capacity");
      }
#endif
    rows = m;
    cols = n;
 }




// =============================================================================
// Assignment helpers
// =============================================================================

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonalAndZeroOffDiagonal( uint32_t diag_elmts, double sc )
 {
    Resize( diag_elmts, diag_elmts );
    for ( uint32_t i{0U}; i < diag_elmts; ++i )
      for ( uint32_t j{0U}; j < diag_elmts; ++j )
        data[i][j] = ( i == j ) ? sc : 0.0;
 }


template<uint32_t mn_max>
void DenseMatrix<mn_max>::RowCondenseTo( vector<double>& vec ) const
 {
    vec.resize( rows );
    for ( uint32_t i{0U}; i < rows; ++i )
      vec[i] = RowSum(i);
 }


// =============================================================================
// Row / column queries
// =============================================================================

template<uint32_t mn_max>
double DenseMatrix<mn_max>::RowSum( uint32_t row ) const
 {
    if ( row >= rows )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::RowSum: row "<< row <<" does not exist" << endl;
        throw length_error("DenseMatrix::RowSum: row index out of range");
      }
    double sum{0.0};
    for ( uint32_t j{0U}; j < cols; ++j )
      sum += data[row][j];
    return sum;
 }


template<uint32_t mn_max>
double DenseMatrix<mn_max>::ColSum( uint32_t col ) const
 {
    if ( col >= cols )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::ColSum: col "<< col <<" does not exist" << endl;
        throw length_error("DenseMatrix::ColSum: col index out of range");
      }
    double sum{0.0};
    for ( uint32_t i{0U}; i < rows; ++i )
      sum += data[i][col];
    return sum;
 }


// =============================================================================
// Norms
// =============================================================================

template<uint32_t mn_max>
double DenseMatrix<mn_max>::NormL1() const noexcept
 {
    if ( rows == 0U || cols == 0U )
      return numeric_limits<double>::signaling_NaN();

    double maxval{0.0};
    for ( uint32_t j{0U}; j < cols; ++j )
      {
        double sum{0.0};
        for ( uint32_t i{0U}; i < rows; ++i )
          sum += std::fabs( data[i][j] );
        maxval = std::max( maxval, sum );
      }
    return maxval;
 }


template<uint32_t mn_max>
double DenseMatrix<mn_max>::NormL_Infinity() const noexcept
 {
    if ( rows == 0U || cols == 0U )
      return numeric_limits<double>::signaling_NaN();

    double maxval{0.0};
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        double sum{0.0};
        for ( uint32_t j{0U}; j < cols; ++j )
          sum += std::fabs( data[i][j] );
        maxval = std::max( maxval, sum );
      }
    return maxval;
 }




template<uint32_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const vector<double>& vec )
 {
    if ( cols != vec.size() )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::operator*=(vector): incompatible sizes" << endl;
        throw length_error("DenseMatrix::operator*=(vector): incompatible sizes");
      }
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        double sum{0.0};
        for ( uint32_t j{0U}; j < cols; ++j )
          sum += data[i][j] * vec[j];
        data[i][0U] = sum;
      }
    cols = 1U;
    return *this;
 }


template<uint32_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const double* vec )
 {
    if ( vec == nullptr )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::operator*=(double*): null pointer" << endl;
        throw length_error("DenseMatrix::operator*=(double*): null pointer");
      }
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        double sum{0.0};
        for ( uint32_t j{0U}; j < cols; ++j )
          sum += data[i][j] * vec[j];
        data[i][0U] = sum;
      }
    cols = 1U;
    return *this;
 }


// =============================================================================
// Transposed products
// =============================================================================

template<uint32_t mn_max>
template<uint32_t other_max1, uint32_t other_max2>
void DenseMatrix<mn_max>::MultiplyWithTransposedOf( const DenseMatrix<other_max1>& A,
                                                    DenseMatrix<other_max2>& res ) const
 {
    if ( cols != A.cols )
      {
        std::cerr <<"\nDenseMatrix<"<< mn_max <<">::MultiplyWithTransposedOf: incompatible sizes" << std::endl;
        throw std::length_error("DenseMatrix::MultiplyWithTransposedOf: incompatible sizes");
      }
    res.Resize( rows, A.rows );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < A.rows; ++j )
        {
          double sum{0.0};
          for ( uint32_t k{0U}; k < cols; ++k )
            sum += data[i][k] * A.data[j][k];
          res.data[i][j] = sum;
        }
 }

template<uint32_t mn_max>
template<uint32_t other_max1, uint32_t other_max2>
void DenseMatrix<mn_max>::MultiplyTransposedOfWith( const DenseMatrix<other_max1>& A,
                                                    DenseMatrix<other_max2>& res ) const
 {
    if ( rows != A.rows )
      {
        std::cerr <<"\nDenseMatrix<"<< mn_max <<">::MultiplyTransposedOfWith: incompatible sizes" << std::endl;
        throw std::length_error("DenseMatrix::MultiplyTransposedOfWith: incompatible sizes");
      }
    res.Resize( cols, A.cols );
    for ( uint32_t i{0U}; i < cols; ++i )
      for ( uint32_t j{0U}; j < A.cols; ++j )
        {
          double sum{0.0};
          for ( uint32_t k{0U}; k < A.rows; ++k )
            sum += data[k][i] * A.data[k][j];
          res.data[i][j] = sum;
        }
 }

// =============================================================================
// CSMP variable operations
// =============================================================================

#ifdef DENSE_MATRIX_USED_TOGETHER_WITH_CSMP

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( uint32_t diag_elmts, const ScalarVariable& sc )
 {
    Resize( diag_elmts, diag_elmts );
    for ( uint32_t i{0U}; i < diag_elmts; ++i )
      for ( uint32_t j{0U}; j < diag_elmts; ++j )
        data[i][j] = ( i == j ) ? sc() : 0.0;
 }

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignRow( uint32_t i, const ArrayVariable& ar )
 {
    assert( ar.Size() < mn_max );
    cols = ar.Size();
    for ( uint32_t j{0U}; j < cols; ++j ) data[i][j] = ar[j];
 }

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignCol( uint32_t j, const ArrayVariable& ar )
 {
    assert( ar.Size() < mn_max );
    rows = ar.Size();
    for ( uint32_t i{0U}; i < rows; ++i ) data[i][j] = ar[i];
 }

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const ArrayVariable& ar )
 {
    assert( ar.Size() < mn_max );
    Resize( ar.Size(), ar.Size() );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = ( i == j ) ? ar[i] : 0.0;
 }

template<uint32_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const ArrayVariable& ar )
 {
    assert( ar.Size() < mn_max );
    if ( cols != ar.Size() )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::operator*=(ArrayVariable): size mismatch" << endl;
        throw range_error("DenseMatrix::operator*=(ArrayVariable): size mismatch");
      }
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        double sum{0.0};
        for ( uint32_t j{0U}; j < cols; ++j ) sum += data[i][j] * ar[j];
        data[i][0] = sum;
      }
    cols = 1U;
    return *this;
 }


template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignRow( uint32_t i, const FlaggedArrayVariable& fr )
 {
    assert( fr.Size() < mn_max );
    cols = fr.Size();
    for ( uint32_t j{0U}; j < cols; ++j ) data[i][j] = fr[j];
 }

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignCol( uint32_t j, const FlaggedArrayVariable& fr )
 {
    assert( fr.Size() < mn_max );
    rows = fr.Size();
    for ( uint32_t i{0U}; i < rows; ++i ) data[i][j] = fr[i];
 }

template<uint32_t mn_max>
void DenseMatrix<mn_max>::AssignToDiagonal( const FlaggedArrayVariable& fr )
 {
    assert( fr.Size() < mn_max );
    Resize( fr.Size(), fr.Size() );
    for ( uint32_t i{0U}; i < rows; ++i )
      for ( uint32_t j{0U}; j < cols; ++j )
        data[i][j] = ( i == j ) ? fr[i] : 0.0;
 }

template<uint32_t mn_max>
DenseMatrix<mn_max>& DenseMatrix<mn_max>::operator*=( const FlaggedArrayVariable& fr )
 {
    assert( fr.Size() < mn_max );
    if ( cols != fr.Size() )
      {
        cerr <<"\nDenseMatrix<"<< mn_max <<">::operator*=(FlaggedArrayVariable): size mismatch" << endl;
        throw range_error("DenseMatrix::operator*=(FlaggedArrayVariable): size mismatch");
      }
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        double sum{0.0};
        for ( uint32_t j{0U}; j < cols; ++j ) sum += data[i][j] * fr[j];
        data[i][0] = sum;
      }
    cols = 1U;
    return *this;
 }



#endif // DENSE_MATRIX_USED_TOGETHER_WITH_CSMP


// =============================================================================
// I/O
// =============================================================================

template<uint32_t mn_max>
void DenseMatrix<mn_max>::In()
 {
    cout <<"\nDenseMatrix<"<< mn_max <<">: enter rows and columns: ";
    cin >> rows >> cols;
    if ( rows > mn_max || cols > mn_max )
      {
        cout <<"\nCapacity exceeded — rebuild the matrix." << endl;
        throw length_error("DenseMatrix::In: capacity exceeded");
      }
    for ( uint32_t i{0U}; i < rows; ++i )
      {
        cout <<"\nRow "<< i+1 <<": ";
        for ( uint32_t j{0U}; j < cols; ++j ) cin >> data[i][j];
      }
    cout <<"\nThank you." << endl;
 }




template<uint32_t mn_max>
void DenseMatrix<mn_max>::Out( long digits ) const
 {
    cout <<"\nDenseMatrix<"<< mn_max <<">: rows="<< rows <<", cols="<< cols << endl;
    cout.setf( ios::scientific );
    const long     prec        = cout.precision( digits );
    const uint32_t stride{     3U };
    const uint32_t split_after = ( digits != 0L ) ? 10U : 18U;

    for ( uint32_t i{0U}; i < rows; ++i )
      {
        uint32_t row_break{1U};
        for ( uint32_t j{0U}; j < cols; ++j, ++row_break )
          {
            if ( std::isnan(data[i][j]) )
              {
                cout << " NaN ";
              }
            else
              {
                const auto len = to_string( std::abs( std::lround(data[i][j]) ) ).length();
                // only pad if the integer part is shorter than stride — avoids unsigned underflow
                if ( len < stride )
                  for ( uint32_t k{0U}; k < stride - len; ++k ) cout << " ";
                if ( data[i][j] >= 0.0 )
                  cout << " " << data[i][j] << " ";
                else
                  cout       << data[i][j]  << " ";
              }

            if ( row_break == split_after ) { cout << endl; row_break = 0U; }
          }
        cout << endl;
      }

    cout.unsetf( ios::scientific );
    cout.precision( prec );
    cout << endl;
 }





// =============================================================================
// Explicit template instantiations
// =============================================================================

template class DenseMatrix<DM1>;
template class DenseMatrix<DM2>;
template class DenseMatrix<DM3>;

template DenseMatrix<DM3> operator+( const DenseMatrix<DM3>&, const DenseMatrix<DM3>& );
template DenseMatrix<DM3> operator-( const DenseMatrix<DM3>&, const DenseMatrix<DM3>& );
template DenseMatrix<DM3> operator*( const DenseMatrix<DM3>&, const DenseMatrix<DM3>& );
template vector<double>   operator*( const DenseMatrix<DM3>&, const vector<double>&   );

template class DenseMatrix<DM4>;

template DenseMatrix<DM4> operator+( const DenseMatrix<DM4>&, const DenseMatrix<DM4>& );
template DenseMatrix<DM4> operator-( const DenseMatrix<DM4>&, const DenseMatrix<DM4>& );
template DenseMatrix<DM4> operator*( const DenseMatrix<DM4>&, const DenseMatrix<DM4>& );
template vector<double>   operator*( const DenseMatrix<DM4>&, const vector<double>&   );

template class DenseMatrix<DM6>;

template DenseMatrix<DM6> operator+( const DenseMatrix<DM6>&, const DenseMatrix<DM6>& );
template DenseMatrix<DM6> operator-( const DenseMatrix<DM6>&, const DenseMatrix<DM6>& );
template DenseMatrix<DM6> operator*( const DenseMatrix<DM6>&, const DenseMatrix<DM6>& );
template vector<double>   operator*( const DenseMatrix<DM6>&, const vector<double>&   );

template class DenseMatrix<DM12>;

template DenseMatrix<DM12> operator+( const DenseMatrix<DM12>&, const DenseMatrix<DM12>& );
template DenseMatrix<DM12> operator-( const DenseMatrix<DM12>&, const DenseMatrix<DM12>& );
template DenseMatrix<DM12> operator*( const DenseMatrix<DM12>&, const DenseMatrix<DM12>& );
template vector<double>    operator*( const DenseMatrix<DM12>&, const vector<double>&    );

template class DenseMatrix<DM_MIN>;

template DenseMatrix<DM_MIN> operator+( const DenseMatrix<DM_MIN>&, const DenseMatrix<DM_MIN>& );
template DenseMatrix<DM_MIN> operator-( const DenseMatrix<DM_MIN>&, const DenseMatrix<DM_MIN>& );
template DenseMatrix<DM_MIN> operator*( const DenseMatrix<DM_MIN>&, const DenseMatrix<DM_MIN>& );
template vector<double>      operator*( const DenseMatrix<DM_MIN>&, const vector<double>&      );

template class DenseMatrix<DM_MAX>;

template DenseMatrix<DM_MAX> operator+( const DenseMatrix<DM_MAX>&, const DenseMatrix<DM_MAX>& );
template DenseMatrix<DM_MAX> operator-( const DenseMatrix<DM_MAX>&, const DenseMatrix<DM_MAX>& );
template DenseMatrix<DM_MAX> operator*( const DenseMatrix<DM_MAX>&, const DenseMatrix<DM_MAX>& );
template vector<double>      operator*( const DenseMatrix<DM_MAX>&, const vector<double>&  );


template void DenseMatrix<DM4>::MultiplyWithTransposedOf( const DenseMatrix<DM4>&, DenseMatrix<DM4>& ) const;
template void DenseMatrix<DM4>::MultiplyTransposedOfWith( const DenseMatrix<DM4>&, DenseMatrix<DM4>& ) const;


} // namespace csmp

