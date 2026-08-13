// VectorVariable2.cpp

#include "VectorVariable2.h"

using namespace std;

namespace csmp {

// ============================================================================
//  Flip — kept in .cpp: constructs a local temporary with swapped indices.
// ============================================================================

VectorVariable<2U> VectorVariable<2U>::Flip() noexcept
{
    VectorVariable<2U> temp;
    temp.flag[0] = flag[1];
    temp.flag[1] = flag[0];
    temp.data[0] = data[1];
    temp.data[1] = data[0];
    return temp;
}

// ============================================================================
//  AngleTo — kept in .cpp: multiple branches and trig call.
// ============================================================================

double VectorVariable<2U>::AngleTo(
    const VectorVariable<2U>& v ) const noexcept
{
    const double ab    = data[0]*v.data[0] + data[1]*v.data[1];
    const double denom = std::sqrt(
        ( data[0]*data[0]   + data[1]*data[1]   ) *
        ( v.data[0]*v.data[0] + v.data[1]*v.data[1] ) );

    if ( denom == 0.0 ) return 90.0;

    const double cos_angle = ab / denom;
    if ( cos_angle >  1.0 ) return   0.0;
    if ( cos_angle < -1.0 ) return 180.0;

    return ( 180.0 / 3.14159265358979323846 ) * std::acos( cos_angle );
}

// ============================================================================
//  Binary I/O
// ============================================================================

bool VectorVariable<2U>::Out( std::fstream& fp ) const
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    const int32_t f0( flag[0] ), f1( flag[1] );
    fp.write( reinterpret_cast<const char*>( &f0 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f1 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &data[0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1] ), data_size );
    return true;
}

bool VectorVariable<2U>::In( std::fstream& fp )
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    fp.read( reinterpret_cast<char*>( &flag[0] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[1] ), flag_size );
    fp.read( reinterpret_cast<char*>( &data[0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1] ), data_size );
    return true;
}

// ============================================================================
//  Interactive I/O
// ============================================================================

void VectorVariable<2U>::In()
{
    string status;
    for ( uint32_t i = 0; i < 2U; ++i )
    {
        if ( i == 0 ) cout << "\nEnter status for x-component: ";
        else          cout << "\nEnter status for y-component: ";
        cout.flush();
        cin >> status;
        flag[i] = parseStatus( status.c_str() );
    }
    cout << "\nEnter x and y vector elements: ";
    cout.flush();
    for ( uint32_t i = 0; i < 2U; ++i ) cin >> data[i];
}

void VectorVariable<2U>::Out() const noexcept
{
    cout << "\nStatus:\n";
    for ( uint32_t i = 0; i < 2U; ++i )
        cout << parseStatus( flag[i] ) << "\t\t";
    cout << "\nValues:\n";
    for ( uint32_t i = 0; i < 2U; ++i )
        cout << data[i] << "\t\t";
    cout << "\n";
}

} // namespace csmp

