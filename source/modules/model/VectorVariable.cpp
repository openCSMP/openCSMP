// VectorVariable.cpp

#include "VectorVariable.h"

#include <iostream>
#include <string>

using namespace std;

namespace csmp {

// ============================================================================
//  makeVector with std::vector arguments — cannot be inline because
//  it uses assert which may not be available in all header contexts.
// ============================================================================

VectorVariable<3U> makeVector(
    const std::vector<VARIABLE_FLAG>& flags,
    const std::vector<double>&        vals ) noexcept
{
    assert( flags.size() == 3U );
    assert( vals.size()  == 3U );
    return VectorVariable<3U>( flags[0], flags[1], flags[2],
                               vals[0],  vals[1],  vals[2] );
}

// ============================================================================
//  Flip — kept in .cpp because it constructs a local temporary with
//  swapped indices which is clearer to read outside the header.
// ============================================================================

VectorVariable<3U> VectorVariable<3U>::Flip() noexcept
{
    VectorVariable<3U> temp;
    temp.flag[0] = flag[2];
    temp.flag[1] = flag[1];
    temp.flag[2] = flag[0];
    temp.data[0] = data[2];
    temp.data[1] = data[1];
    temp.data[2] = data[0];
    return temp;
}

// ============================================================================
//  AngleTo — kept in .cpp due to multiple branches and trig call.
// ============================================================================

double VectorVariable<3U>::AngleTo(
    const VectorVariable<3U>& v ) const noexcept
{
    const double ab = data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2];
    const double denom = std::sqrt(
        ( data[0]*data[0]   + data[1]*data[1]   + data[2]*data[2] ) *
        ( v.data[0]*v.data[0] + v.data[1]*v.data[1] + v.data[2]*v.data[2] ) );

    if ( denom == 0.0 ) return 90.0;

    const double cos_angle = ab / denom;
    if ( cos_angle >  1.0 ) return   0.0;
    if ( cos_angle < -1.0 ) return 180.0;

    return ( 180.0 / 3.14159265358979323846 ) * std::acos( cos_angle );
}

// ============================================================================
//  Binary operators — 1D and 2D (kept in .cpp since 1D/2D types are
//  defined in VectorVariable1.h / VectorVariable2.h which are included
//  by the header but whose inline definitions live in their own files)
// ============================================================================

VectorVariable<1U> operator+( const VectorVariable<1U>& a, double b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) + b ); }

VectorVariable<1U> operator-( const VectorVariable<1U>& a, double b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) - b ); }

VectorVariable<1U> operator*( const VectorVariable<1U>& a, double b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) * b ); }

VectorVariable<1U> operator/( const VectorVariable<1U>& a, double b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) / b ); }

VectorVariable<2U> operator+( const VectorVariable<2U>& a, double b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)+b, a(1)+b ); }

VectorVariable<2U> operator-( const VectorVariable<2U>& a, double b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)-b, a(1)-b ); }

VectorVariable<2U> operator*( const VectorVariable<2U>& a, double b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)*b, a(1)*b ); }

VectorVariable<2U> operator/( const VectorVariable<2U>& a, double b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)/b, a(1)/b ); }

// ============================================================================
//  Binary operators — ScalarVariable (1D and 2D)
// ============================================================================

VectorVariable<1U> operator+( const VectorVariable<1U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) + b() ); }

VectorVariable<1U> operator-( const VectorVariable<1U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) - b() ); }

VectorVariable<1U> operator*( const VectorVariable<1U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) * b() ); }

VectorVariable<1U> operator/( const VectorVariable<1U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<1U>( a.Flag(0), a(0) / b() ); }

VectorVariable<2U> operator+( const VectorVariable<2U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)+b(), a(1)+b() ); }

VectorVariable<2U> operator-( const VectorVariable<2U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)-b(), a(1)-b() ); }

VectorVariable<2U> operator*( const VectorVariable<2U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)*b(), a(1)*b() ); }

VectorVariable<2U> operator/( const VectorVariable<2U>& a, const ScalarVariable& b ) noexcept
{ return VectorVariable<2U>( a.Flag(0), a.Flag(1), a(0)/b(), a(1)/b() ); }

// ============================================================================
//  Point op VectorVariable (1D and 2D)
// ============================================================================

Point<1U> operator+( const Point<1U>& p, const VectorVariable<1U>& vc ) noexcept
{ return Point<1U>( p[0] + vc[0] ); }

Point<1U> operator-( const Point<1U>& p, const VectorVariable<1U>& vc ) noexcept
{ return Point<1U>( p[0] - vc[0] ); }

Point<1U> operator*( const Point<1U>& p, const VectorVariable<1U>& vc ) noexcept
{ return Point<1U>( p[0] * vc[0] ); }

Point<1U> operator/( const Point<1U>& p, const VectorVariable<1U>& vc ) noexcept
{ return Point<1U>( p[0] / vc[0] ); }

Point<2U> operator+( const Point<2U>& p, const VectorVariable<2U>& vc ) noexcept
{ return Point<2U>( p[0]+vc[0], p[1]+vc[1] ); }

Point<2U> operator-( const Point<2U>& p, const VectorVariable<2U>& vc ) noexcept
{ return Point<2U>( p[0]-vc[0], p[1]-vc[1] ); }

Point<2U> operator*( const Point<2U>& p, const VectorVariable<2U>& vc ) noexcept
{ return Point<2U>( p[0]*vc[0], p[1]*vc[1] ); }

Point<2U> operator/( const Point<2U>& p, const VectorVariable<2U>& vc ) noexcept
{ return Point<2U>( p[0]/vc[0], p[1]/vc[1] ); }

// ============================================================================
//  Binary I/O
// ============================================================================

bool VectorVariable<3U>::Out( std::fstream& fp ) const
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    const int32_t f0(flag[0]), f1(flag[1]), f2(flag[2]);
    fp.write( reinterpret_cast<const char*>( &f0 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f1 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &f2 ), flag_size );
    fp.write( reinterpret_cast<const char*>( &data[0] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[1] ), data_size );
    fp.write( reinterpret_cast<const char*>( &data[2] ), data_size );
    return true;
}

bool VectorVariable<3U>::In( std::fstream& fp )
{
    const size_t flag_size = sizeof(int32_t);
    const size_t data_size = sizeof(double);
    fp.read( reinterpret_cast<char*>( &flag[0] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[1] ), flag_size );
    fp.read( reinterpret_cast<char*>( &flag[2] ), flag_size );
    fp.read( reinterpret_cast<char*>( &data[0] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[1] ), data_size );
    fp.read( reinterpret_cast<char*>( &data[2] ), data_size );
    return true;
}

// ============================================================================
//  Interactive I/O
// ============================================================================

void VectorVariable<3U>::In()
{
    string status;
    for ( uint32_t i = 0; i < 3U; ++i )
    {
        if      ( i == 0 ) cout << "\nEnter status for x-component: ";
        else if ( i == 1 ) cout << "\nEnter status for y-component: ";
        else               cout << "\nEnter status for z-component: ";
        cout.flush();
        cin >> status;
        flag[i] = parseStatus( status.c_str() );
    }
    cout << "\nEnter 3 vector elements: ";
    cout.flush();
    for ( uint32_t i = 0; i < 3U; ++i ) cin >> data[i];
}

void VectorVariable<3U>::Out() const noexcept
{
    cout << "\nStatus:\n";
    for ( uint32_t i = 0; i < 3U; ++i )
        cout << parseStatus( flag[i] ) << "\t\t";
    cout << "\nValues:\n";
    for ( uint32_t i = 0; i < 3U; ++i )
        cout << data[i] << "\t\t";
    cout << "\n";
}

// ============================================================================
//  Stream output operator (template — explicit instantiations)
// ============================================================================

template<uint32_t dim>
ostream& operator<<( ostream& stream, const VectorVariable<dim>& o )
{
    for ( uint32_t i = 0; i < dim; ++i )
        stream << o[i] << " (" << parseStatus( o.Flag(i) ) << ") ";
    return stream;
}

template ostream& operator<< <1U>( ostream&, const VectorVariable<1U>& );
template ostream& operator<< <2U>( ostream&, const VectorVariable<2U>& );
template ostream& operator<< <3U>( ostream&, const VectorVariable<3U>& );

} // namespace csmp

