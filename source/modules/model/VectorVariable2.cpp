#include "VectorVariable2.h"

using namespace std;

namespace csmp {

VectorVariable<2U>::VectorVariable() noexcept
  : flag{ { ANY,ANY } },
  data{ { std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN() } }
{
}


VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f, double val ) noexcept
  : flag{ { f,f } },
  data{ { val,val } }
{
}


VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double val1, double val2 ) noexcept
  : flag{ { f1,f2 } },
  data{ { val1,val2 } }
{
}


VectorVariable<2U>::VectorVariable( const std::vector<double>& v ) noexcept
  : flag{ { ANY,ANY } },
  data{ { v[0],v[1] } }
{
}


VectorVariable<2U>::VectorVariable( const Point<2U>& p ) noexcept
  : flag{ { ANY,ANY } },
  data{ { p[0],p[1] } }
{
}



double& VectorVariable<2U>::operator()( uint32_t i ) noexcept
{
  if ( i == 0U ) return data[0];
  return              data[1];
}


const double& VectorVariable<2U>::operator()( uint32_t i ) const noexcept
{
  if ( i == 0U ) return data[0];
  return              data[1];
}



double  VectorVariable<2U>::operator[]( uint32_t i ) const noexcept
{
  if ( i == 0U ) return data[0];
  return              data[1];
}



void  VectorVariable<2U>::Component( uint32_t i, double val ) noexcept
{
  if ( i == 0U ) data[0] = val;
  else         data[1] = val;
}



double  VectorVariable<2U>::Component( uint32_t i ) const noexcept
{
  if ( i == 0U ) return data[0];
  return              data[1];
}



VectorVariable<2U>  VectorVariable<2U>::operator+( const VectorVariable<2U>& v ) const noexcept
{
  return VectorVariable( flag[0], flag[1], data[0] + v.data[0], data[1] + v.data[1] );
}




VectorVariable<2U>  VectorVariable<2U>::operator-( const VectorVariable<2U>& v ) const noexcept
{
  return VectorVariable( flag[0], flag[1], data[0] - v.data[0], data[1] - v.data[1] );
}



VectorVariable<2U>  VectorVariable<2U>::operator*( const VectorVariable<2U>& v ) const noexcept
{
  return VectorVariable( flag[0], flag[1], data[0] * v.data[0], data[1] * v.data[1] );
}



VectorVariable<2U>  VectorVariable<2U>::operator/( const VectorVariable<2U>& v ) const noexcept
{
  return VectorVariable( flag[0], flag[1], data[0] / v.data[0], data[1] / v.data[1] );
}



VectorVariable<2U>  VectorVariable<2U>::operator^( double val ) const noexcept
{
  return VectorVariable( flag[0], flag[1], std::pow( data[0], val ), std::pow( data[1], val ) );
}




VectorVariable<2U>&  VectorVariable<2U>::operator+=( double val ) noexcept
{
  data[0] += val;
  data[1] += val;

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator-=( double val ) noexcept
{
  data[0] -= val;
  data[1] -= val;

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator*=( double val ) noexcept
{
  data[0] *= val;
  data[1] *= val;

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator/=( double val ) noexcept
{
  data[0] /= val;
  data[1] /= val;

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator+=( const ScalarVariable& sc ) noexcept
{
  data[0] += sc();
  data[1] += sc();

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator-=( const ScalarVariable& sc ) noexcept
{
  data[0] -= sc();
  data[1] -= sc();

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator*=( const ScalarVariable& sc ) noexcept
{
  data[0] *= sc();
  data[1] *= sc();

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator/=( const ScalarVariable& sc ) noexcept
{
  data[0] /= sc();
  data[1] /= sc();

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator+=( const VectorVariable<2U>& v ) noexcept
{
  data[0] += v.data[0];
  data[1] += v.data[1];

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator-=( const VectorVariable<2U>& v ) noexcept
{
  data[0] -= v.data[0];
  data[1] -= v.data[1];

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator*=( const VectorVariable<2U>& v ) noexcept
{
  data[0] *= v.data[0];
  data[1] *= v.data[1];

  return *this;
}




VectorVariable<2U>&  VectorVariable<2U>::operator/=( const VectorVariable<2U>& v ) noexcept
{
  data[0] /= v.data[0];
  data[1] /= v.data[1];

  return *this;
}




// --------------------
// ASSIGNMENT OPERATORS
// --------------------

VectorVariable<2U>&  VectorVariable<2U>::operator=( double val ) noexcept
{
  data[0] = val;
  data[1] = val;

  return *this;
}



VectorVariable<2U>&  VectorVariable<2U>::operator=( const csmp::Point<2U>& p ) noexcept
{
  data[0] = p[0];
  data[1] = p[1];

  return *this;
}



VectorVariable<2U>&  VectorVariable<2U>::operator=( const ScalarVariable& sc ) noexcept
{
  flag[0] = flag[1] = sc.Flag();
  data[0] = data[1] = sc();

  return *this;
}


bool VectorVariable<2U>::operator==( const VectorVariable<2U>& v ) const noexcept
{
  return(flag == v.flag && data == v.data);
}


bool VectorVariable<2U>::operator!=( const VectorVariable<2U>& v ) const noexcept
{
  return(flag != v.flag || data != v.data);
}


// compare the length of two vectors

bool VectorVariable<2U>::operator<( const VectorVariable<2U>& v ) const noexcept
{
  return (this < &v);
}





// -------
// METHODS
// -------

/// L2 norm
void VectorVariable<2U>::EuclideanNormalize() noexcept
{
  const double fNorm( std::hypot( data[0], data[1] ) );

  if ( fNorm == 0. ) return; //added AP

  data[0] /= fNorm;
  data[1] /= fNorm;
}



double VectorVariable<2U>::DotProduct( const csmp::Point<2U>& p ) const noexcept
{
  return data[0] * p[0] + data[1] * p[1];
}

double VectorVariable<2U>::DotProduct( const VectorVariable& v ) const noexcept
{
  return data[0] * v[0] + data[1] * v[1];
}


VectorVariable<2U> VectorVariable<2U>::CrossProduct( const csmp::Point<2U>& p ) const noexcept
{
  return VectorVariable<2U>( flag[0], flag[1], 0., data[0] * p[1] - data[1] * p[0] );
}

VectorVariable<2U> VectorVariable<2U>::CrossProduct( const VectorVariable& v ) const noexcept
{
  return VectorVariable<2U>( flag[0], flag[1], 0., data[0] * v[1] - data[1] * v[0] );
}


double  VectorVariable<2U>::Length() const noexcept
{
  return std::hypot( data[0], data[1] );
}



VARIABLE_FLAG&  VectorVariable<2U>::Flag( uint32_t i ) noexcept
{
  if ( i == 0U ) return flag[0];
  return flag[1];
}


VARIABLE_FLAG  VectorVariable<2U>::Flag( uint32_t i ) const noexcept
{
  if ( i == 0U ) return flag[0];
  return flag[1];
}


Point<2U>  VectorVariable<2U>::P() const noexcept
{
  return csmp::Point<2U>( data[0], data[1] );
}



bool  VectorVariable<2U>::IsWithinRange( double vmin, double vmax ) const noexcept
{
  if ( data[0] < vmin || data[0] > vmax ) return false;
  if ( data[1] < vmin || data[1] > vmax ) return false;

  return true;
}


  /// tests whether variable contains NaN value(s)
bool  VectorVariable<2U>::Has_NaN_Values() const noexcept {
   if ( isnan(data[0]) ) return true;
   if ( isnan(data[1]) ) return true;
   return false;
}



bool VectorVariable<2U>::Out( std::fstream& fp ) const
{
  const int32_t flag_0( this->flag[0] );
  const int32_t flag_1( this->flag[1] );
  const size_t flag_size = sizeof( int32_t );
  fp.write( (char*)&flag_0, flag_size );
  fp.write( (char*)&flag_1, flag_size );
  const size_t data_size = sizeof( double );
  fp.write( (char*)&this->data[0], data_size );
  fp.write( (char*)&this->data[1], data_size );
  return true;
}

bool VectorVariable<2U>::In( std::fstream& fp )
{
  const size_t flag_size = sizeof( int32_t );
  fp.read( (char*)&this->flag[0], flag_size );
  fp.read( (char*)&this->flag[1], flag_size );
  const size_t data_size = sizeof( double );
  fp.read( (char*)&this->data[0], data_size );
  fp.read( (char*)&this->data[1], data_size );
  return true;
}

/** return angle in degrees

*/
double  VectorVariable<2U>::AngleTo( const VectorVariable<2U>& v ) const noexcept
{
  double ab, a_dot_b;

  // a b
  // ---
  ab = data[0] * v.data[0] + data[1] * v.data[1];
  // |a| . |b|
  // ---------
  a_dot_b = std::sqrt( (data[0] * data[0] + data[1] * data[1]) *
                       (v.data[0] * v.data[0] + v.data[1] * v.data[1]) );

  // a b
  // ---
  double cos_angle = ab / a_dot_b;

  // if zero intercept
  if ( cos_angle == 0.0 ) return  90.0;
  // if outside of range of 'acos' function
  if ( cos_angle >  1.0 ) return   0.0;
  if ( cos_angle < -1.0 ) return 180.0;

  return (180.0 / 3.14159265358979324) * std::acos( cos_angle );
}



VectorVariable<2U>  VectorVariable<2U>::Flip() noexcept
{
  VectorVariable<2U>  temp;

  temp.flag[0] = flag[1];
  temp.flag[1] = flag[0];
  temp.data[0] = data[1];
  temp.data[1] = data[0];

  return temp;
}

/// Multiplies by negative unity vector
void  VectorVariable<2U>::Invert() noexcept
{
  data[0] *= -1.;
  data[1] *= -1.;
}

VectorVariable<2U>  VectorVariable<2U>::ProjectOnto( const std::vector<double>& v ) const noexcept
{
  const double ratio( (data[0] * v[0] + data[1] * v[1]) / (v[0] * v[0] + v[1] * v[1]) );

  return VectorVariable<2U>( flag[0], flag[1], v[0] * ratio, v[1] * ratio );
}



VectorVariable<2U>  VectorVariable<2U>::ProjectOnto( const VectorVariable<2U>& v ) const noexcept
{
  const double ratio( (data[0] * v.data[0] + data[1] * v.data[1]) / (v.data[0] * v.data[0] + v.data[1] * v.data[1]) );

  return VectorVariable<2U>( flag[0], flag[1], v.data[0] * ratio, v.data[1] * ratio );
}



void  VectorVariable<2U>::In()
{
  string  status;

  cout.flush();
  for ( uint32_t i = 0; i<2U; i++ )
  {
    if ( i == 0 ) cout << "\nEnter status for x-component of variable: ";
    else          cout << "\nEnter status for y-component of variable: ";
    cout.flush();
    cin >> status;
    flag[i] = parseStatus( status.c_str() );
  }

  cout << "\nEnter x=0 and y=1 vector variable elements: ";
  cout.flush();
  for ( uint32_t i = 0; i<2U; i++ ) cin >> data[i];

} // end In




void  VectorVariable<2U>::Out() const noexcept
{
  string  status;

  cout << "\nStatus: " << endl;
  for ( uint32_t i = 0; i<2U; i++ )
    cout << (status = parseStatus( flag[i] )) << "\t\t";
  cout << endl;
  for ( uint32_t i = 0; i<2U; i++ ) cout << data[i] << "\t\t";
  cout << endl;

} // end Out


} // end namespace csmp

