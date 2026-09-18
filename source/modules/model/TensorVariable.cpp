// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "TensorVariable.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/// fastest way to insert a tensor into an STL container; tensor only has flags for diagonal elements
TensorVariable<2U> makeTensor( VARIABLE_FLAG f1, VARIABLE_FLAG f2,
                               double v11, double v12,
                               double v21, double v22 ) noexcept
{
  return TensorVariable<2U>( f1, f2, v11, v12, v21, v22 );
}



// ============================================================================
//  Interactive I/O
// ============================================================================

template<uint32_t dim>
ostream&  operator<<( ostream& stream, const TensorVariable<dim>& o )
{
  for ( uint32_t i{0U}; i<dim; i++ )
  {
    stream << "flag" << i + 1U << ": " << parseStatus( o.Flag( i ) );
    stream << ", row" << i + 1U << ":";
    for ( uint32_t j{0U}; j<dim; j++ )
      stream <<" "<< o( i, j );
    stream <<"; ";

  }

  return stream;
}



void  TensorVariable<3U>::In()
{
  size_t  i;
  string  status;

  cout.flush();
  cout << "\nEnter [" << 3U << "] tensor variable status: ";
  cout.flush();
  cin >> status;
  for ( i = 0; i<3U; i++ )
    flag[i] = parseStatus( status.c_str() );

  cout << "\nEnter first row of elements : ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[0][i];
  cout << "Enter second row of elements: ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[1][i];
  cout << "Enter third row of elements : ";
  cout.flush();
  for ( i = 0; i<3U; i++ ) cin >> data[2][i];

} // end In





  /// @test tested: O.K.
void  TensorVariable<3U>::Out() const
{
  size_t   i, j;

  cout << "\nStatus: " << endl;
  for ( i = 0; i<3U; i++ )
    cout << parseStatus( flag[i] ) << "  ";

  cout << endl;

  cout << "\nValues: " << endl;
  for ( i = 0; i<3U; i++ )
  {
    for ( j = 0; j<3U; j++ ) cout << data[i][j] << "\t\t";
    cout << endl;
  }
} // end Out



// ============================================================================
//  Binary I/O
// ============================================================================

bool TensorVariable<3U>::Out( fstream& fp ) const
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    const int32_t f0( flag[0] ), f1( flag[1] ), f2( flag[2] );
    fp.write( reinterpret_cast<const char*>( &f0 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f1 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f2 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &data[0][0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1][0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[2][0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[0][1] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1][1] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[2][1] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[0][2] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1][2] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[2][2] ), data_size );
    return true;
}

bool TensorVariable<3U>::In( fstream& fp )
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    fp.read( reinterpret_cast<char*>( &flag[0] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[1] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[2] ), flag_size );
    fp.read( reinterpret_cast<char*>( &data[0][0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1][0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[2][0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[0][1] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1][1] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[2][1] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[0][2] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1][2] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[2][2] ), data_size );
    return true;
}




  /**

  Computes eigenvalues and eigenvectors assuming that the tensor variable
  is symmetric. If not symmetric, off-diagonal elements are averaged.
  No check of symmetry is performed.

  @attention the Eigenvectors are returned row by row into the evecs
  tensor.

  By default the length of the eigenvectors is equivalent to the eigenvalues
  but it can be normalized to one.

  @section arguments Input Arguments

  Flag to normalize length of eigenvectors to 1.

  @return Eigenvalues and vectors are returned into the supplied Vector and
  TensorVariables, respectively. Method returns false if the rank of the matrix is zero.

  @section implementation Implementation

  Cubic root finding algorithm using Cardano's formula and the eigenvectors
  are calculated via LU-backsubstitution.

  @test SKM - observed, by comparison with Maple, that all the values are
  the same, but that the first element of the
  second Eigenvector always has the opposite sign (does this matter?)

  */
bool TensorVariable<3U>::Eigen( VectorVariable<3U>& evals, TensorVariable<3U>& evecs, bool bNormalize ) const
{
  //calculate eigenvalues
  //set eigenvalues
  if ( !EigenValues( evals ) ) return false;

  //suppose tensor is symmetric, by averaging the non diagonal elements
  const double f_0_1( (data[1][0] + data[0][1]) / 2. ); //=1_0
  const double f_0_2( (data[2][0] + data[0][2]) / 2. ); //=2_0
  const double f_1_2( (data[1][2] + data[2][1]) / 2. ); //=1_2

                                                          //EigenValueCalculation(Ev0_, Ev1_, Ev2_);
  TensorVariable<3U> Ev0_I( ANY, evals( 0 ) );
  TensorVariable<3U> Ev1_I( ANY, evals( 1 ) );
  TensorVariable<3U> Ev2_I( ANY, evals( 2 ) );

  TensorVariable<3U> A( ANY, data[0][0], f_0_1, f_0_2, f_0_1, data[1][1], f_1_2, f_0_2, f_1_2, data[2][2] );
  //A.Out();

  TensorVariable<3U> A_Ev0_I = A - Ev0_I;
  TensorVariable<3U> A_Ev1_I = A - Ev1_I;
  TensorVariable<3U> A_Ev2_I = A - Ev2_I;

  TensorVariable<3U> AEv0byAEv1 = A_Ev0_I * A_Ev1_I;
  TensorVariable<3U> AEv0byAEv2 = A_Ev0_I * A_Ev2_I;
  TensorVariable<3U> AEv1byAEv2 = A_Ev1_I * A_Ev2_I;

  VectorVariable<3U> LEvec( ANY, ANY, ANY,
                            sqrt( AEv1byAEv2( 0, 0 )*AEv1byAEv2( 0, 0 ) + AEv1byAEv2( 1, 0 )*AEv1byAEv2( 1, 0 ) + AEv1byAEv2( 2, 0 )*AEv1byAEv2( 2, 0 ) ),
                            sqrt( AEv0byAEv2( 0, 1 )*AEv0byAEv2( 0, 1 ) + AEv0byAEv2( 1, 1 )*AEv0byAEv2( 1, 1 ) + AEv0byAEv2( 2, 1 )*AEv0byAEv2( 2, 1 ) ),
                            sqrt( AEv0byAEv1( 0, 2 )*AEv0byAEv1( 0, 2 ) + AEv0byAEv1( 1, 2 )*AEv0byAEv1( 1, 2 ) + AEv0byAEv1( 2, 2 )*AEv0byAEv1( 2, 2 ) ) );

  uint32_t EvecLengthColID[3] = { 0, 1, 2 };

  for ( uint32_t i( 0 ); i < 3; i++ ) {
    uint32_t j = i;
    if ( LEvec( i ) < numeric_limits<double>::epsilon() ) {
      j++;
      if ( j > 2 ) { j = 0; }
      LEvec( i ) = sqrt( AEv1byAEv2( 0, j )*AEv1byAEv2( 0, j ) + AEv1byAEv2( 1, j )*AEv1byAEv2( 1, j ) + AEv1byAEv2( 2, j )*AEv1byAEv2( 2, j ) );
      if ( LEvec( i ) < numeric_limits<double>::epsilon() ) {
        j++;
        if ( j > 2 ) { j = 0; }
        LEvec( i ) = sqrt( AEv1byAEv2( 0, j )*AEv1byAEv2( 0, j ) + AEv1byAEv2( 1, j )*AEv1byAEv2( 1, j ) + AEv1byAEv2( 2, j )*AEv1byAEv2( 2, j ) );
      }
      EvecLengthColID[i] = j;
    }
  }

  // Normalised eigen vectors. Each column of the following matrix is the eigenvector for the corresponding eigenvalue.
  for ( uint32_t i( 0 ); i < 3; i++ ) {
    evecs( i, 0 ) = AEv1byAEv2( i, EvecLengthColID[0] ) / LEvec( EvecLengthColID[0] );
    evecs( i, 1 ) = AEv0byAEv2( i, EvecLengthColID[1] ) / LEvec( EvecLengthColID[1] );
    evecs( i, 2 ) = AEv0byAEv1( i, EvecLengthColID[2] ) / LEvec( EvecLengthColID[2] );
  }

  // Not Normalised eignvectors; i.e. vectors of principal stresses/strains.
  if ( !bNormalize ) {
    for ( uint32_t j( 0 ); j< 3; j++ ) {
      for ( uint32_t i( 0 ); i < 3; i++ ) {
        evecs( i, j ) *= evals( j );
      }
    }
  }

  // evecs inherits the flags from the original tensor variable
  evecs.flag[0] = flag[0];
  evecs.flag[1] = flag[1];
  evecs.flag[2] = flag[2];

  return true;

} // end Eigen





  /**

  Computes eigenvalues and eigenvectors for 3d TensorVariable (doesn't have to be symmetric)

  @section arguments Input Arguments

  eigenVectors, eigenValues

  @return Eigenvalues and vectors are returned into the supplied Vector and
  TensorVariables, respectively.

  @section implementation Implementation
  Uses:
  1) Symmetric Householder reduction to tridiagonal form,
  derived from the Algol procedures tred2 by
  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  Fortran subroutine in EISPACK.
  2) Symmetric tridiagonal QL algorithm,
  derived from the Algol procedures tql2, by
  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
  Fortran subroutine in EISPACK.

  @author Lukas Mosser?

  */
bool TensorVariable<3U>::EigenSymmetric( VectorVariable<3U>& eigenVals,
                                         TensorVariable<3U>& eigenVecs ) const
{
    // This method implements the Householder tridiagonalisation + symmetric
    // tridiagonal QL algorithm.  It requires the input tensor to be symmetric.
    // Non-diagonal elements are symmetrised by averaging before processing.
    //
    // Eigenvalues are returned in descending order (largest to smallest),
    // consistent with EigenValuesPositiveDefiniteSymmetricMatrix.
    //
    // Source: Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    // Auto. Comp., Vol.ii-Linear Algebra (EISPACK).

    constexpr int n = 3;
    double V[3][3], d[3], e[3];

    // Copy tensor into V, symmetrising off-diagonal elements.
    for ( uint32_t i = 0; i < n; ++i )
        for ( uint32_t j = 0; j < n; ++j )
            V[i][j] = ( (*this)(i,j) + (*this)(j,i) ) / 2.0;

    // ------------------------------------------------------------------
    //  Symmetric Householder reduction to tridiagonal form.
    // ------------------------------------------------------------------
    for ( int j = 0; j < n; ++j )
        d[j] = V[n-1][j];

    for ( int i = n-1; i > 0; --i )
    {
        double scale = 0.0;
        double h     = 0.0;

        for ( int k = 0; k < i; ++k )
            scale += std::fabs( d[k] );

        if ( scale == 0.0 )
        {
            e[i] = d[i-1];
            for ( int j = 0; j < i; ++j )
            {
                d[j]    = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        else
        {
            for ( int k = 0; k < i; ++k )
            {
                d[k] /= scale;
                h    += d[k] * d[k];
            }
            double f = d[i-1];
            double g = std::sqrt( h );
            if ( f > 0.0 ) g = -g;
            e[i]    = scale * g;
            h       = h - f * g;
            d[i-1]  = f - g;

            for ( int j = 0; j < i; ++j )
                e[j] = 0.0;

            for ( int j = 0; j < i; ++j )
            {
                f       = d[j];
                V[j][i] = f;
                g       = e[j] + V[j][j] * f;
                for ( int k = j+1; k <= i-1; ++k )
                {
                    g    += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }

            f = 0.0;
            for ( int j = 0; j < i; ++j )
            {
                e[j] /= h;
                f    += e[j] * d[j];
            }
            double hh = f / ( h + h );
            for ( int j = 0; j < i; ++j )
                e[j] -= hh * d[j];

            for ( int j = 0; j < i; ++j )
            {
                f = d[j];
                g = e[j];
                for ( int k = j; k <= i-1; ++k )
                    V[k][j] -= ( f * e[k] + g * d[k] );
                d[j]    = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // ------------------------------------------------------------------
    //  Accumulate transformations.
    // ------------------------------------------------------------------
    for ( int i = 0; i < n-1; ++i )
    {
        V[n-1][i] = V[i][i];
        V[i][i]   = 1.0;
        double h  = d[i+1];
        if ( h != 0.0 )
        {
            for ( int k = 0; k <= i; ++k )
                d[k] = V[k][i+1] / h;
            for ( int j = 0; j <= i; ++j )
            {
                double g = 0.0;
                for ( int k = 0; k <= i; ++k )
                    g += V[k][i+1] * V[k][j];
                for ( int k = 0; k <= i; ++k )
                    V[k][j] -= g * d[k];
            }
        }
        for ( int k = 0; k <= i; ++k )
            V[k][i+1] = 0.0;
    }
    for ( int j = 0; j < n; ++j )
    {
        d[j]        = V[n-1][j];
        V[n-1][j]   = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0]        = 0.0;

    // ------------------------------------------------------------------
    //  Symmetric tridiagonal QL algorithm.
    // ------------------------------------------------------------------
    for ( int i = 1; i < n; ++i )
        e[i-1] = e[i];
    e[n-1] = 0.0;

    double f    = 0.0;
    double tst1 = 0.0;
    const double eps = std::numeric_limits<double>::epsilon();

    for ( int l = 0; l < n; ++l )
    {
        tst1 = std::max( tst1, std::fabs( d[l] ) + std::fabs( e[l] ) );

        // Find small subdiagonal element.
        int m = l;
        while ( m < n )
        {
            if ( std::fabs( e[m] ) <= eps * tst1 ) break;
            ++m;
        }

        if ( m > l )
        {
            int iter = 0;
            do
            {
                ++iter;

                // Compute implicit shift.
                double g  = d[l];
                double p  = ( d[l+1] - g ) / ( 2.0 * e[l] );
                double r  = std::sqrt( p*p + 1.0 );
                if ( p < 0.0 ) r = -r;
                d[l]      = e[l] / ( p + r );
                d[l+1]    = e[l] * ( p + r );
                double dl1 = d[l+1];
                double h   = g - d[l];
                for ( int i = l+2; i < n; ++i )
                    d[i] -= h;
                f += h;

                // Implicit QL transformation.
                p          = d[m];
                double c   = 1.0;
                double c2  = c;
                double c3  = c;
                double el1 = e[l+1];
                double s   = 0.0;
                double s2  = 0.0;

                for ( int i = m-1; i >= l; --i )
                {
                    c3          = c2;
                    c2          = c;
                    s2          = s;
                    g           = c * e[i];
                    h           = c * p;
                    r           = std::sqrt( p*p + e[i]*e[i] );
                    e[i+1]      = s * r;
                    s           = e[i] / r;
                    c           = p / r;
                    p           = c * d[i] - s * g;
                    d[i+1]      = h + s * ( c * g + s * d[i] );

                    for ( int k = 0; k < n; ++k )
                    {
                        h           = V[k][i+1];
                        V[k][i+1]   = s * V[k][i] + c * h;
                        V[k][i]     = c * V[k][i] - s * h;
                    }
                }
                p    = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

            } while ( std::fabs( e[l] ) > eps * tst1 );
        }
        d[l] += f;
        e[l]  = 0.0;
    }

    // ------------------------------------------------------------------
    //  Sort eigenvalues and eigenvectors in descending order.
    //  (Consistent with EigenValuesPositiveDefiniteSymmetricMatrix.)
    // ------------------------------------------------------------------
    for ( int i = 0; i < n-1; ++i )
    {
        int    k = i;
        double p = d[i];
        for ( int j = i+1; j < n; ++j )
        {
            if ( d[j] > p )   // > for descending (was < in original)
            {
                k = j;
                p = d[j];
            }
        }
        if ( k != i )
        {
            d[k] = d[i];
            d[i] = p;
            for ( int j = 0; j < n; ++j )
            {
                p       = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }

    // ------------------------------------------------------------------
    //  Copy results into output variables.
    // ------------------------------------------------------------------
    for ( uint32_t i = 0; i < n; ++i )
    {
        eigenVals(i) = d[i];
        for ( uint32_t j = 0; j < n; ++j )
            eigenVecs(i,j) = V[i][j];
    }

    return true;
}




/**
Computes eigenvalues and eigenvectors for a weakly non-symmetric 3D tensor.

A weakly non-symmetric tensor is one whose skew-symmetric part is small
relative to its symmetric part — for example a permeability or stress tensor
that has accumulated small numerical asymmetries during computation.

@section method Method

The tensor A is decomposed into its symmetric part S and skew-symmetric
part W:

    S = (A + A^T) / 2
    W = (A - A^T) / 2

Eigenvalues and eigenvectors are computed from S using the Householder
tridiagonalisation + symmetric tridiagonal QL algorithm (EigenSymmetric).

The Frobenius norm of W relative to S is computed as the asymmetry measure:

    asymmetry = ||W||_F / ||S||_F

If this ratio exceeds @p tolerance the method returns false and issues a
WARNING, indicating that the non-symmetric part is too large for the
symmetric approximation to be reliable.  The eigenvalues are still written
to the output arguments so the caller can inspect them, but they should be
treated with caution.

Eigenvalues are returned in descending order (largest to smallest),
consistent with EigenValuesPositiveDefiniteSymmetricMatrix and
EigenSymmetric.

@param eigenVals   Output: eigenvalues in descending order.
@param eigenVecs   Output: corresponding eigenvectors as columns.
@param tolerance   Maximum permitted ratio ||W||_F / ||S||_F.
                   Default is 1e-6.  Increase for more permissive checks.

@return true  if the asymmetry ratio is within tolerance.
        false if the asymmetry ratio exceeds tolerance (result unreliable).
*/
bool TensorVariable<3U>::EigenWeaklyNonSymmetric( VectorVariable<3U>& eigenVals, TensorVariable<3U>& eigenVecs, double tolerance ) const
{
    constexpr int n = 3;

    // ------------------------------------------------------------------
    //  Decompose into symmetric (S) and skew-symmetric (W) parts.
    // ------------------------------------------------------------------
    double S[3][3], W[3][3];

    for ( uint32_t i = 0; i < n; ++i )
        for ( uint32_t j = 0; j < n; ++j )
        {
            S[i][j] = ( (*this)(i,j) + (*this)(j,i) ) / 2.0;
            W[i][j] = ( (*this)(i,j) - (*this)(j,i) ) / 2.0;
        }

    // ------------------------------------------------------------------
    //  Compute Frobenius norms of S and W.
    //  ||A||_F = sqrt( sum_ij a_ij^2 )
    // ------------------------------------------------------------------
    double normS = 0.0;
    double normW = 0.0;

    for ( int i = 0; i < n; ++i )
        for ( int j = 0; j < n; ++j )
        {
            normS += S[i][j] * S[i][j];
            normW += W[i][j] * W[i][j];
        }
    normS = std::sqrt( normS );
    normW = std::sqrt( normW );

    // ------------------------------------------------------------------
    //  Check asymmetry ratio.
    //  If normS is effectively zero the tensor is purely skew-symmetric
    //  (all eigenvalues are zero) — handle as a special case.
    // ------------------------------------------------------------------
    const double asymmetry = ( normS > std::numeric_limits<double>::epsilon() )
                             ? normW / normS
                             : normW;

    bool within_tolerance = ( asymmetry <= tolerance );

    if ( !within_tolerance )
    {
        ErrorHandler::Instance().Note( WARNING,
            "TensorVariable<3U>::EigenWeaklyNonSymmetric",
            ( string("Asymmetry ratio ||W||_F / ||S||_F = ")
              + std::to_string( asymmetry )
              + " exceeds tolerance "
              + std::to_string( tolerance )
              + ". Eigenvalues computed from symmetric part only "
              + "and may not be representative of the full tensor."
            ).c_str() );
    }

    // ------------------------------------------------------------------
    //  Build a TensorVariable from S and delegate to EigenSymmetric.
    // ------------------------------------------------------------------
    TensorVariable<3U> symTensor;
    for ( uint32_t i = 0; i < n; ++i )
        for ( uint32_t j = 0; j < n; ++j )
            symTensor(i,j) = S[i][j];

    symTensor.EigenSymmetric( eigenVals, eigenVecs );

    return within_tolerance;
}





/**

Computes the Eigen values for a positive definite symmetric matrix.
A source for solving the cubic equation is: http://mathworld.wolfram.com/CubicFormula.html

*/
bool TensorVariable<3U>::EigenValuesPositiveDefiniteSymmetricMatrix(
    double& eigenValue0,
    double& eigenValue1,
    double& eigenValue2 ) const
{
    constexpr double PI = 3.14159265358979323846;

    const double f_0_1 = ( data[1][0] + data[0][1] ) / 2.0;
    const double f_0_2 = ( data[2][0] + data[0][2] ) / 2.0;
    const double f_1_2 = ( data[1][2] + data[2][1] ) / 2.0;

    const double a2 = -( data[0][0] + data[1][1] + data[2][2] );

    const double a1 = -( f_0_1 * f_0_1
                       + f_0_2 * f_0_2
                       + f_1_2 * f_1_2
                       - data[0][0] * data[1][1]
                       - data[0][0] * data[2][2]
                       - data[1][1] * data[2][2] );

    const double a0 = -( data[0][0] * data[1][1] * data[2][2]
                       - data[0][0] * f_1_2 * f_1_2
                       - f_0_1 * f_0_1 * data[2][2]
                       + f_0_1 * f_0_2 * f_1_2
                       + f_0_2 * f_0_1 * f_1_2
                       - f_0_2 * f_0_2 * data[1][1] );

    const double Q = ( 3.0 * a1 - a2 * a2 ) / 9.0;
    const double R = ( 9.0 * a2 * a1 - 27.0 * a0 - 2.0 * a2 * a2 * a2 ) / 54.0;
    const double D = Q * Q * Q + R * R;

    // ------------------------------------------------------------------
    //  Tolerance for branch selection.
    //
    //  D = Q^3 + R^2 suffers catastrophic cancellation when the two
    //  terms are nearly equal in magnitude (repeated eigenvalue case).
    //  We classify D as zero when |D| is small relative to the scale
    //  of the terms that compose it, using sqrt(epsilon) as the
    //  relative tolerance to give a wide enough band.
    // ------------------------------------------------------------------
    const double D_scale = std::fabs( Q * Q * Q ) + std::fabs( R * R );
    const double D_tol   = D_scale
                         * std::sqrt( std::numeric_limits<double>::epsilon() );

    auto sortDescending = []( double& v0, double& v1, double& v2 )
    {
        if ( v1 > v0 ) std::swap( v0, v1 );
        if ( v2 > v1 ) std::swap( v1, v2 );
        if ( v1 > v0 ) std::swap( v0, v1 );
    };

    if ( D > D_tol )
    {
        // ------------------------------------------------------------------
        //  One real root and a pair of complex conjugate roots.
        //  For a symmetric matrix this should not occur; return false.
        // ------------------------------------------------------------------
        double m = R + std::sqrt( D );
        double n = R - std::sqrt( D );
        m = ( m >= 0.0 ) ?  cbrt(  m ) : -cbrt( -m );
        n = ( n >= 0.0 ) ?  cbrt(  n ) : -cbrt( -n );
        eigenValue0 = m + n - a2 / 3.0;
        eigenValue1 = eigenValue2 = 0.0;
        return false;
    }
    else if ( std::fabs( D ) <= D_tol )
    {
        // ------------------------------------------------------------------
        //  Three real roots, at least two equal.
        // ------------------------------------------------------------------
        eigenValue0 =  2.0 * cbrt( R ) - a2 / 3.0;
        eigenValue1 = -cbrt( R ) - a2 / 3.0;
        eigenValue2 = eigenValue1;
        sortDescending( eigenValue0, eigenValue1, eigenValue2 );
        return true;
    }
    else
    {
        // ------------------------------------------------------------------
        //  Three distinct real roots  (D < -D_tol).
        // ------------------------------------------------------------------
        const double sqrtNegQ = std::sqrt( -Q );
        const double denom    = sqrtNegQ * sqrtNegQ * sqrtNegQ;

        // Guard against domain error in acos due to rounding.
        const double cosArg = std::max( -1.0,
                              std::min(  1.0, R / denom ) );
        const double theta  = std::acos( cosArg );
        const double two    = 2.0 * sqrtNegQ;

        eigenValue0 = two * std::cos( theta / 3.0 )
                    - a2 / 3.0;
        eigenValue1 = two * std::cos( theta / 3.0 + 2.0 * PI / 3.0 )
                    - a2 / 3.0;
        eigenValue2 = two * std::cos( theta / 3.0 + 4.0 * PI / 3.0 )
                    - a2 / 3.0;

        sortDescending( eigenValue0, eigenValue1, eigenValue2 );
        return true;
    }
}



template ostream&  operator<<( ostream&, const TensorVariable<1U>& );
template ostream&  operator<<( ostream&, const TensorVariable<2U>& );
template ostream&  operator<<( ostream&, const TensorVariable<3U>& );

} // end namespace csmp

