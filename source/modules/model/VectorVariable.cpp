#include "VectorVariable.h"

using namespace std;

namespace csmp {


VectorVariable<3U>::VectorVariable()
  : flag{ { ANY,ANY,ANY } },
  data{ { std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN() } }
{
}





VectorVariable<3U>::VectorVariable( VARIABLE_FLAG f, double val )
  : flag{ { f,f,f } },
  data{ { val,val,val } }
{
}




VectorVariable<3U>::VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
                                    double  val1, double  val2, double  val3 )
  : flag{ { f1,f2,f3 } },
  data{ { val1,val2,val3 } }
{
}




// 1D specializations
VectorVariable<1U> makeVector( VARIABLE_FLAG fx, double vx )
{
  return VectorVariable<1U>( fx, vx );
}

// 2D specializations
VectorVariable<2U> makeVector( VARIABLE_FLAG fx, VARIABLE_FLAG fy, double vx, double vy )
{
  return VectorVariable<2U>( fx, fy, vx, vy );
}


// 3D specializations

VectorVariable<3U> makeVector( VARIABLE_FLAG fx, VARIABLE_FLAG fy, VARIABLE_FLAG fz, double vx, double vy, double vz )
{
  return VectorVariable<3U>( fx, fy, fz, vx, vy, vz );
}

VectorVariable<3U> makeVector( const std::array<VARIABLE_FLAG, 3U>& flags, const std::array<double, 3U>& vals )
{
  return VectorVariable<3U>( flags[0], flags[1], flags[2], vals[0], vals[1], vals[2] );
}

VectorVariable<3U> makeVector( const std::vector<VARIABLE_FLAG>& flags, const std::vector<double>& vals )
{
  assert( flags.size() == 3U );
  assert( vals.size() == 3U );
  return VectorVariable<3U>( flags[0], flags[1], flags[2], vals[0], vals[1], vals[2] );
}



double& VectorVariable<3U>::operator()( uint32_t i )
{
#ifdef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nVectorVariable<3U>::operator(): vector access violation, i=" << i << std::endl;
    return data[0];
  }
#endif
  return data[i];
}

const double& VectorVariable<3U>::operator()( uint32_t i ) const
{
#ifdef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nVectorVariable<3U>::operator() const: vector access violation, i=" << i << std::endl;
    return data[0];
  }
#endif
  return data[i];
}


double  VectorVariable<3U>::operator[]( uint32_t i ) const
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nVectorVariable<3U>::operator[]: vector access violation, i=" << i << std::endl;
    return data[0];
  }
#endif
  return data[i];
}



void  VectorVariable<3U>::Component( uint32_t i, double val )
{
  assert( i < 3U );
  data[i] = val;
}




double  VectorVariable<3U>::Component( uint32_t i ) const
{
  assert( i < 3U );
  return data[i];
}




VARIABLE_FLAG&  VectorVariable<3U>::Flag( uint32_t i )
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nVectorVariable<3U>::Flag(): access violation, i=" << i << std::endl;
    return flag[0];
  }
#endif
  return flag[i];
}


VARIABLE_FLAG   VectorVariable<3U>::Flag( uint32_t i ) const
{
#ifndef NDEBUG 
  if ( i >= 3U ) {
    std::cerr << "\nVectorVariable<3U>::Flag(): access violation, i=" << i << std::endl;
    return flag[0];
  }
#endif
  return flag[i];
}

uint32_t VectorVariable<3U>::Size() const
{
  return 3U;
}



VectorVariable<3U>::VectorVariable( const std::vector<double>& v )
  : flag{ { ANY,ANY,ANY } },
  data{ { v[0],v[1],v[2] } }
{
}


VectorVariable<3U>::VectorVariable( const csmp::Point<3U>& p )
  : flag{ { ANY,ANY,ANY } },
  data{ { p[0],p[1],p[2] } }
{
}




VectorVariable<3U>&  VectorVariable<3U>::operator+=( double val )
{
  data[0] += val;
  data[1] += val;
  data[2] += val;

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator-=( double val )
{
  data[0] -= val;
  data[1] -= val;
  data[2] -= val;

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator*=( double val )
{
  data[0] *= val;
  data[1] *= val;
  data[2] *= val;

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator/=( double val )
{
  data[0] /= val;
  data[1] /= val;
  data[2] /= val;

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator+=( const ScalarVariable& sc )
{
  data[0] += sc();
  data[1] += sc();
  data[2] += sc();

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator-=( const ScalarVariable& sc )
{
  data[0] -= sc();
  data[1] -= sc();
  data[2] -= sc();

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator*=( const ScalarVariable& sc )
{
  data[0] *= sc();
  data[1] *= sc();
  data[2] *= sc();

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator/=( const ScalarVariable& sc )
{
  data[0] /= sc();
  data[1] /= sc();
  data[2] /= sc();

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator+=( const VectorVariable<3U>& v )
{
  data[0] += v.data[0];
  data[1] += v.data[1];
  data[2] += v.data[2];

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator-=( const VectorVariable<3U>& v )
{
  data[0] -= v.data[0];
  data[1] -= v.data[1];
  data[2] -= v.data[2];

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator*=( const VectorVariable<3U>& v )
{
  data[0] *= v.data[0];
  data[1] *= v.data[1];
  data[2] *= v.data[2];

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator/=( const VectorVariable<3U>& v )
{
  data[0] /= v.data[0];
  data[1] /= v.data[1];
  data[2] /= v.data[2];

  return *this;
}




// --------------------
// ASSIGNMENT OPERATORS
// --------------------

VectorVariable<3U>&  VectorVariable<3U>::operator=( double val )
{
  data[0] = val;
  data[1] = val;
  data[2] = val;

  return *this;
}




VectorVariable<3U>&  VectorVariable<3U>::operator=( const ScalarVariable& sc )
{
  flag[0] = flag[1] = flag[2] = sc.Flag();
  data[0] = data[1] = data[2] = sc();

  return *this;
}



VectorVariable<3U>&  VectorVariable<3U>::operator=( const csmp::Point<3U>& pt )
{
  data[0] = pt[0];
  data[1] = pt[1];
  data[2] = pt[2];

  return *this;
}

// using the comparitor of the standard array
bool VectorVariable<3U>::operator==( const VectorVariable<3U>& v ) const
{
  return(flag == v.flag && data == v.data);
}


bool VectorVariable<3U>::operator!=( const VectorVariable<3U>& v ) const
{
  return(flag != v.flag || data != v.data);
}


/// compare the length of two vectors
bool VectorVariable<3U>::operator<( const VectorVariable<3U>& v ) const
{
  return (this < &v);
}



// -------
// METHODS
// -------

/// normalize L2

void VectorVariable<3U>::EuclideanNormalize()
{
  const double fNorm( std::sqrt( data[0] * data[0] + data[1] * data[1] + data[2] * data[2] ) );

  if ( fNorm == 0. ) return; //added AP

  data[0] /= fNorm;
  data[1] /= fNorm;
  data[2] /= fNorm;
}


double VectorVariable<3U>::DotProduct( const csmp::Point<3U>& p ) const
{
  return data[0] * p[0] + data[1] * p[1] + data[2] * p[2];
}

double VectorVariable<3U>::DotProduct( const VectorVariable& v ) const
{
  return data[0] * v[0] + data[1] * v[1] + data[2] * v[2];
}

VectorVariable<3U> VectorVariable<3U>::CrossProduct( const csmp::Point<3U>& p ) const
{
  return VectorVariable<3U>( flag[0], flag[1], flag[2],
                             data[1] * p[2] - data[2] * p[1],
                             data[2] * p[0] - data[0] * p[2],
                             data[0] * p[1] - data[1] * p[0] );
}

VectorVariable<3U> VectorVariable<3U>::CrossProduct( const VectorVariable& v ) const
{
  return VectorVariable<3U>( flag[0], flag[1], flag[2],
                             data[1] * v[2] - data[2] * v[1],
                             data[2] * v[0] - data[0] * v[2],
                             data[0] * v[1] - data[1] * v[0] );
}

double  VectorVariable<3U>::Length() const
{
  return std::hypot( data[0], data[1], data[2] );
}



Point<3U>  VectorVariable<3U>::P() const
{
  return csmp::Point<3U>( data[0], data[1], data[2] );
}



bool  VectorVariable<3U>::IsWithinRange( double vmin, double vmax ) const
{
  if ( data[0] < vmin || data[0] > vmax ) return false;
  if ( data[1] < vmin || data[1] > vmax ) return false;
  if ( data[2] < vmin || data[2] > vmax ) return false;

  return true;
}



  /// tests whether variable contains NaN value(s)
bool  VectorVariable<3U>::Has_NaN_Values() const {
   if ( isnan(data[0]) ) return true;
   if ( isnan(data[1]) ) return true;
   if ( isnan(data[2]) ) return true;
   return false;
}



bool VectorVariable<3U>::Out( std::fstream& fp ) const
{
  const int32_t flag_0( this->flag[0] );
  const int32_t flag_1( this->flag[1] );
  const int32_t flag_2( this->flag[2] );
  const size_t flag_size = sizeof( int32_t );
  fp.write( (char*)&flag_0, flag_size );
  fp.write( (char*)&flag_1, flag_size );
  fp.write( (char*)&flag_2, flag_size );
  const size_t data_size = sizeof( double );
  fp.write( (char*)&this->data[0], data_size );
  fp.write( (char*)&this->data[1], data_size );
  fp.write( (char*)&this->data[2], data_size );
  return true;
}

bool VectorVariable<3U>::In( std::fstream& fp )
{
  const size_t flag_size = sizeof( int32_t );
  fp.read( (char*)&this->flag[0], flag_size );
  fp.read( (char*)&this->flag[1], flag_size );
  fp.read( (char*)&this->flag[2], flag_size );
  const size_t data_size = sizeof( double );
  fp.read( (char*)&this->data[0], data_size );
  fp.read( (char*)&this->data[1], data_size );
  fp.read( (char*)&this->data[2], data_size );
  return true;
}

/// operators with Points
// 1D

Point<1U> operator+( const Point<1U>& p, const VectorVariable<1U>& vc )
{
  return Point<1U>( p[0] + vc[0] );
}


Point<1U> operator-( const Point<1U>& p, const VectorVariable<1U>& vc )
{
  return Point<1U>( p[0] - vc[0] );
}


Point<1U> operator*( const Point<1U>& p, const VectorVariable<1U>& vc )
{
  return Point<1U>( p[0] * vc[0] );
}


Point<1U> operator/( const Point<1U>& p, const VectorVariable<1U>& vc )
{
  return Point<1U>( p[0] / vc[0] );
}

// 2D

Point<2U> operator+( const Point<2U>& p, const VectorVariable<2U>& vc )
{
  return Point<2U>( p[0] + vc[0], p[1] + vc[1] );
}


Point<2U> operator-( const Point<2U>& p, const VectorVariable<2U>& vc )
{
  return Point<2U>( p[0] - vc[0], p[1] - vc[1] );
}


Point<2U> operator*( const Point<2U>& p, const VectorVariable<2U>& vc )
{
  return Point<2U>( p[0] * vc[0], p[1] * vc[1] );
}


Point<2U> operator/( const Point<2U>& p, const VectorVariable<2U>& vc )
{
  return Point<2U>( p[0] / vc[0], p[1] / vc[1] );
}

// 3D

Point<3U> operator+( const Point<3U>& p, const VectorVariable<3U>& vc )
{
  return Point<3U>( p[0] + vc[0], p[1] + vc[1], p[2] + vc[2] );
}


Point<3U> operator-( const Point<3U>& p, const VectorVariable<3U>& vc )
{
  return Point<3U>( p[0] - vc[0], p[1] - vc[1], p[2] - vc[2] );
}


Point<3U> operator*( const Point<3U>& p, const VectorVariable<3U>& vc )
{
  return Point<3U>( p[0] * vc[0], p[1] * vc[1], p[2] * vc[2] );
}


Point<3U> operator/( const Point<3U>& p, const VectorVariable<3U>& vc )
{
  return Point<3U>( p[0] / vc[0], p[1] / vc[1], p[2] / vc[2] );
}



VectorVariable<3U>  VectorVariable<3U>::operator+( const VectorVariable<3U>& v ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] + v.data[0], data[1] + v.data[1], data[2] + v.data[2] );
}




VectorVariable<3U>  VectorVariable<3U>::operator-( const VectorVariable<3U>& v ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] - v.data[0], data[1] - v.data[1], data[2] - v.data[2] );
}



VectorVariable<3U>  VectorVariable<3U>::operator*( const VectorVariable<3U>& v ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] * v.data[0], data[1] * v.data[1], data[2] * v.data[2] );
}



VectorVariable<3U>  VectorVariable<3U>::operator/( const VectorVariable<3U>& v ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] / v.data[0], data[1] / v.data[1], data[2] / v.data[2] );
}



VectorVariable<3U>  VectorVariable<3U>::operator+( double val ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] + val, data[1] + val, data[2] + val );
}



VectorVariable<3U>  VectorVariable<3U>::operator-( double val ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] - val, data[1] - val, data[2] - val );
}



VectorVariable<3U>  VectorVariable<3U>::operator*( double val ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] * val, data[1] * val, data[2] * val );
}



VectorVariable<3U>  VectorVariable<3U>::operator/( double val ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[0] / val, data[1] / val, data[2] / val );
}




VectorVariable<3U>  VectorVariable<3U>::operator^( double val ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    std::pow( data[0], val ), std::pow( data[1], val ), std::pow( data[2], val ) );
}




/// cross product ' % ' of two vectors
VectorVariable<3U>  VectorVariable<3U>::operator%( const VectorVariable<3U>& v ) const
{
  return VectorVariable( flag[0], flag[1], flag[2],
                    data[1] * v.data[2] - v.data[1] * data[2],
                    -(data[0] * v.data[2] - v.data[0] * data[2]),
                    data[0] * v.data[1] - v.data[0] * data[1] );
}




/// return angle in degrees
double  VectorVariable<3U>::AngleTo( const VectorVariable<3U>& v ) const
{
  double ab = data[0] * v.data[0] + data[1] * v.data[1] + data[2] * v.data[2];
  double a_dot_b = std::sqrt( (data[0] * data[0] + data[1] * data[1] + data[2] * data[2]) *
                                (v.data[0] * v.data[0] + v.data[1] * v.data[1] + v.data[2] * v.data[2]) );
  // a b
  // ---
  double cos_angle = ab / a_dot_b;

  // if zero intercept
  if ( cos_angle == 0.0 ) return  90.0;
  // if outside of range of 'acos' function
  if ( cos_angle >  1.0 ) return   0.0;
  if ( cos_angle < -1.0 ) return 180.0;

  return (static_cast<double>(180.) / static_cast<double>(3.14159265358979324)) * std::acos( cos_angle );
}



VectorVariable<3U>  VectorVariable<3U>::Flip()
{
  VectorVariable<3U>  temp;

  temp.flag[0] = flag[2];
  temp.flag[2] = flag[0];
  temp.flag[1] = flag[1];
  temp.data[0] = data[2];
  temp.data[2] = data[0];
  temp.data[1] = data[1];

  return temp;
}


/// Multiplies by negative unity vector
void  VectorVariable<3U>::Invert()
{
  data[0] *= -1.;
  data[1] *= -1.;
  data[2] *= -1.;
}




/// projects this vector variable onto the supplied vector v and returns the projection vector

VectorVariable<3U>  VectorVariable<3U>::ProjectOnto( const std::vector<double>& v ) const
{
  double ratio( (data[0] * v[0] + data[1] * v[1] + data[2] * v[2]) / (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) );

  return VectorVariable<3U>( flag[0], flag[1], flag[2], v[0] * ratio, v[1] * ratio, v[2] * ratio );
}



VectorVariable<3U>  VectorVariable<3U>::ProjectOnto( const VectorVariable<3U>& v ) const
{
  double ratio( (data[0] * v.data[0] + data[1] * v.data[1] + data[2] * v.data[2]) /
                  (v.data[0] * v.data[0] + v.data[1] * v.data[1] + v.data[2] * v.data[2]) );

  return VectorVariable<3U>( flag[0], flag[1], flag[2],
                    v.data[0] * ratio, v.data[1] * ratio, v.data[2] * ratio );
}



template<uint32_t dim>
ostream&  operator<<( ostream& stream, const VectorVariable<dim>& o )
{
  for ( auto i{0U}; i<dim; i++ )
    stream << o[i] << " (" << parseStatus( o.Flag( i ) ) << ") ";

  return stream;
}



void  VectorVariable<3U>::In()
{
  string  status;

  for ( auto i = 0; i<3U; i++ ) {
    if ( i == 0 )      cout << "\nEnter status for x-component of variable: ";
    else if ( i == 1 ) cout << "\nEnter status for y-component of variable: ";
    else               cout << "\nEnter status for z-component of variable: ";
    cout.flush();
    cin >> status;
    flag[i] = parseStatus( status.c_str() );
  }

  cout << "\nEnter 3 vector elements: ";
  cout.flush();
  for ( auto i = 0; i<3U; i++ ) cin >> data[i];

} // end In



void  VectorVariable<3U>::Out() const
{
  cout << "\nStatus: " << endl;
  for ( auto i = 0; i<3U; i++ )
    cout << parseStatus( flag[i] ) << "\t\t";
  cout << endl;
  for ( auto i = 0; i<3U; i++ ) cout << data[i] << "\t\t";
  cout << endl;

} // end Out


template ostream&  operator<< <1U>(ostream& stream, const VectorVariable<1U>& o);
template ostream&  operator<< <2U>(ostream& stream, const VectorVariable<2U>& o);
template ostream&  operator<< <3U>(ostream& stream, const VectorVariable<3U>& o);

} // end namespace csmp

