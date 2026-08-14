// TensorVariable2.cpp

#include "TensorVariable2.h"
#include "ErrorHandler.h"

#include <iostream>
#include <string>

using namespace std;

namespace csmp {

// ============================================================================
//  Interactive I/O
// ============================================================================

void TensorVariable<2U>::In()
{
    cout << "\nEnter [2] tensor variable status: ";
    cout.flush();
    string status;
    cin >> status;
    for ( uint32_t i = 0; i < 2U; ++i )
        flag[i] = parseStatus( status.c_str() );

    cout << "\nEnter first row of elements : ";
    cout.flush();
    for ( uint32_t i = 0; i < 2U; ++i ) cin >> data[0][i];
    cout << "Enter second row of elements: ";
    cout.flush();
    for ( uint32_t i = 0; i < 2U; ++i ) cin >> data[1][i];
}

void TensorVariable<2U>::Out() const noexcept
{
    cout << "\nStatus:\n";
    for ( uint32_t j = 0; j < 2U; ++j )
        cout << parseStatus( flag[j] ) << "  ";
    cout << "\nValues:\n";
    for ( uint32_t i = 0; i < 2U; ++i )
    {
        for ( uint32_t j = 0; j < 2U; ++j )
            cout << data[i][j] << "\t\t";
        cout << "\n";
    }
}

// ============================================================================
//  Binary I/O
// ============================================================================

bool TensorVariable<2U>::Out( fstream& fp ) const
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    const int32_t f0( flag[0] ), f1( flag[1] );
    fp.write( reinterpret_cast<const char*>( &f0 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f1 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &data[0][0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1][0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[0][1] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1][1] ), data_size );
    return true;
}

bool TensorVariable<2U>::In( fstream& fp )
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    fp.read( reinterpret_cast<char*>( &flag[0] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[1] ), flag_size );
    fp.read( reinterpret_cast<char*>( &data[0][0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1][0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[0][1] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1][1] ), data_size );
    return true;
}

// ============================================================================
//  Eigenvalue / eigenvector methods
// ============================================================================

bool TensorVariable<2U>::EigenValues( VectorVariable<2U>& ev ) const
{
    const double a1 = -( data[0][0] + data[1][1] );
    const double a0 =    data[0][0] * data[1][1]
                       - data[1][0] * data[0][1];
    const double D  = a1 * a1 - 4.0 * a0;

    if ( D >= 0.0 )
    {
        const double sqrtD = std::sqrt( D );
        ev(0) = ( -a1 + sqrtD ) / 2.0;
        ev(1) = ( -a1 - sqrtD ) / 2.0;
        if ( ev(1) > ev(0) ) std::swap( ev(0), ev(1) );
        return true;
    }
    else
    {
        ev(0) = std::numeric_limits<double>::quiet_NaN();
        ev(1) = std::numeric_limits<double>::quiet_NaN();
        ErrorHandler::Instance().Note( WARNING,
            "TensorVariable<2U>::EigenValues",
            "Discriminant is negative: eigenvalues are complex. "
            "NaN returned. Check that the tensor is symmetric "
            "positive definite." );
        return false;
    }
}

bool TensorVariable<2U>::EigenValues( std::vector<double>& ev ) const
{
    if ( ev.size() < 2 )
        ev.resize( 2, 0.0 );
    VectorVariable<2U> vv;
    const bool ok = EigenValues( vv );
    ev[0] = vv(0);
    ev[1] = vv(1);
    return ok;
}

bool TensorVariable<2U>::EigenSymmetric(
    VectorVariable<2U>& eigenVals,
    TensorVariable<2U>& eigenVecs,
    bool                bNormalize ) const
{
    const double f01 = ( data[0][1] + data[1][0] ) / 2.0;

    TensorVariable<2U> sym;
    sym(0,0) = data[0][0];
    sym(1,1) = data[1][1];
    sym(0,1) = sym(1,0) = f01;

    if ( !sym.EigenValues( eigenVals ) )
        return false;

    if ( std::fabs( f01 ) <= std::numeric_limits<double>::epsilon()
                             * std::max( std::fabs( data[0][0] ),
                                         std::fabs( data[1][1] ) ) )
    {
        eigenVecs(0,0) = 1.0;  eigenVecs(0,1) = 0.0;
        eigenVecs(1,0) = 0.0;  eigenVecs(1,1) = 1.0;
        return true;
    }

    {
        VectorVariable<2U> v;
        v(0) = -( data[1][1] - eigenVals(0) ) / f01;
        v(1) = 1.0;
        if ( bNormalize ) v.EuclideanNormalize();
        eigenVecs(0,0) = v(0);
        eigenVecs(1,0) = v(1);
    }
    {
        VectorVariable<2U> v;
        v(0) = -( data[1][1] - eigenVals(1) ) / f01;
        v(1) = 1.0;
        if ( bNormalize ) v.EuclideanNormalize();
        eigenVecs(0,1) = v(0);
        eigenVecs(1,1) = v(1);
    }

    return true;
}

bool TensorVariable<2U>::EigenWeaklyNonSymmetric(
    VectorVariable<2U>& eigenVals,
    TensorVariable<2U>& eigenVecs,
    double              tolerance ) const
{
    const double S00 =   data[0][0];
    const double S11 =   data[1][1];
    const double S01 = ( data[0][1] + data[1][0] ) / 2.0;
    const double W01 = ( data[0][1] - data[1][0] ) / 2.0;

    const double normS = std::sqrt( S00*S00 + S11*S11 + 2.0*S01*S01 );
    const double normW = std::sqrt( 2.0 * W01 * W01 );

    const double asymmetry = ( normS > std::numeric_limits<double>::epsilon() )
                             ? normW / normS
                             : normW;

    const bool within_tolerance = ( asymmetry <= tolerance );

    if ( !within_tolerance )
        ErrorHandler::Instance().Note( WARNING,
            "TensorVariable<2U>::EigenWeaklyNonSymmetric",
            ( std::string("Asymmetry ratio ||W||_F / ||S||_F = ")
              + std::to_string( asymmetry )
              + " exceeds tolerance "
              + std::to_string( tolerance )
              + ". Eigenvalues computed from symmetric part only "
              + "and may not be representative of the full tensor."
            ).c_str() );

    TensorVariable<2U> symTensor;
    symTensor(0,0) = S00;
    symTensor(1,1) = S11;
    symTensor(0,1) = symTensor(1,0) = S01;

    symTensor.EigenSymmetric( eigenVals, eigenVecs, true );

    return within_tolerance;
}

} // namespace csmp

