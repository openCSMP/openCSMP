#include "VectorVariable1.h"

using namespace std;

namespace csmp {

VectorVariable<1U>::VectorVariable()
  : flag( ANY ), data( std::numeric_limits<double>::quiet_NaN() )
{
}



VectorVariable<1U>::VectorVariable( const VectorVariable& vc )
  : flag( vc.flag ), data( vc.data )
{
}




VectorVariable<1U>::~VectorVariable()
{
}




VectorVariable<1U>&  VectorVariable<1U>::operator=( const VectorVariable<1U>& v )
{
  if ( &v != this ) {
    flag = v.flag;
    data = v.data;
  }
  return *this;
}



double& VectorVariable<1U>::operator()( uint32_t )
{
  return data;
}


const double& VectorVariable<1U>::operator()( uint32_t ) const
{
  return data;
}


double  VectorVariable<1U>::operator[]( uint32_t ) const
{
  return data;
}



void  VectorVariable<1U>::Component( uint32_t, double val )
{
  data = val;
}



double  VectorVariable<1U>::Component( uint32_t ) const
{
  return data;
}

VARIABLE_FLAG& VectorVariable<1U>::Flag( uint32_t )
{
  return flag;
}

VARIABLE_FLAG VectorVariable<1U>::Flag( uint32_t ) const
{
  return flag;
}

uint32_t VectorVariable<1U>::Size() const
{
  return 1U;
}

VectorVariable<1U>::VectorVariable( VARIABLE_FLAG f, double val )
  : flag( f ), data( val )
{
}




VectorVariable<1U>::VectorVariable( const std::vector<double>& v )
  : flag( ANY ), data( v[0] )
{
}

VectorVariable<1U>::VectorVariable( const csmp::Point<1U>& p )
  : flag( ANY ), data( p[0] )
{
}



VectorVariable<1U>  VectorVariable<1U>::operator+( const VectorVariable<1U>& v ) const
{
  return VectorVariable( flag, data + v.data );
}




VectorVariable<1U>  VectorVariable<1U>::operator-( const VectorVariable<1U>& v ) const
{
  return VectorVariable( flag, data - v.data );
}



VectorVariable<1U>  VectorVariable<1U>::operator*( const VectorVariable<1U>& v ) const
{
  return VectorVariable( flag, data * v.data );
}



VectorVariable<1U>  VectorVariable<1U>::operator/( const VectorVariable<1U>& v ) const
{
  return VectorVariable( flag, data / v.data );
}



VectorVariable<1U>  VectorVariable<1U>::operator+( double val ) const
{
  return VectorVariable( flag, data + val );
}



VectorVariable<1U>  VectorVariable<1U>::operator-( double val ) const
{
  return VectorVariable( flag, data - val );
}



VectorVariable<1U>  VectorVariable<1U>::operator*( double val ) const
{
  return VectorVariable( flag, data * val );
}



VectorVariable<1U>  VectorVariable<1U>::operator/( double val ) const
{
  return VectorVariable( flag, data / val );
}




VectorVariable<1U>  VectorVariable<1U>::operator^( double val ) const
{
  return VectorVariable( flag, std::pow( data, val ) );
}




VectorVariable<1U>&  VectorVariable<1U>::operator+=( double val )
{
  data += val;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator-=( double val )
{
  data -= val;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator*=( double val )
{
  data *= val;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator/=( double val )
{
  data /= val;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator+=( const ScalarVariable& sc )
{
  data += sc();

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator-=( const ScalarVariable& sc )
{
  data -= sc();

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator*=( const ScalarVariable& sc )
{
  data *= sc();

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator/=( const ScalarVariable& sc )
{
  data /= sc();

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator+=( const VectorVariable<1U>& v )
{
  data += v.data;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator-=( const VectorVariable<1U>& v )
{
  data -= v.data;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator*=( const VectorVariable<1U>& v )
{
  data *= v.data;

  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator/=( const VectorVariable<1U>& v )
{
  data /= v.data;

  return *this;
}




// --------------------
// ASSIGNMENT OPERATORS
// --------------------
VectorVariable<1U>&  VectorVariable<1U>::operator=( double val )
{
  data = val;
  return *this;
}




VectorVariable<1U>&  VectorVariable<1U>::operator=( const csmp::Point<1U>& p )
{
  data = p[0];
  return *this;
}



VectorVariable<1U>&  VectorVariable<1U>::operator=( const ScalarVariable& sc )
{
  flag = sc.Flag();
  data = sc();

  return *this;
}



bool VectorVariable<1U>::operator==( const VectorVariable<1U>& v ) const
{
  return((v.flag == flag && v.data == data));
}



bool VectorVariable<1U>::operator!=( const VectorVariable<1U>& v ) const
{
  return(v.flag != flag || v.data != data);
}



/// compare magnitude

bool VectorVariable<1U>::operator<( const VectorVariable<1U>& v ) const
{
  return(this < &v);
}





// -------
// METHODS
// -------


double VectorVariable<1U>::DotProduct( const csmp::Point<1U>& p ) const
{
  return data * p[0];
}

double VectorVariable<1U>::DotProduct( const VectorVariable& v ) const
{
  return data * v[0];
}


// not defined in 1D -> degenerate result = 0.

VectorVariable<1U> VectorVariable<1U>::CrossProduct( const csmp::Point<1U>& ) const
{
  return VectorVariable<1U>( flag, 0. );
}

VectorVariable<1U> VectorVariable<1U>::CrossProduct( const VectorVariable& ) const
{
  return VectorVariable<1U>( flag, 0. );
}


VectorVariable<1U>   VectorVariable<1U>::ProjectOnto( const std::vector<double>& v ) const
{
  return VectorVariable<1U>( flag, v[0] );
}



VectorVariable<1U>   VectorVariable<1U>::ProjectOnto( const VectorVariable& v ) const
{
  return VectorVariable<1U>( flag, v.data );
}



double  VectorVariable<1U>::Length() const
{
  return std::fabs( data );
}


VectorVariable<1U>  VectorVariable<1U>::Flip()
{
  return VectorVariable<1U>( flag, -data );
}



Point<1U>  VectorVariable<1U>::P() const
{
  return csmp::Point<1U>( data );
}




bool  VectorVariable<1U>::IsWithinRange( double vmin, double vmax ) const
{
  if ( data < vmin || data > vmax ) return false;

  return true;
}


bool VectorVariable<1U>::Out( std::fstream& fp ) const
{
  const int32_t flag_0( flag );
  fp.write( (char*)&flag_0, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.write( (char*)&data, sizeof( double ) );
  return true;
}

bool VectorVariable<1U>::In( std::fstream& fp )
{
  fp.read( (char*)&flag, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.read( (char*)&data, sizeof( double ) );
  return true;
}


void  VectorVariable<1U>::In()
{
  string  status;
  cout << "\nEnter status for x-component of variable: ";
  cin >> status;
  flag = parseStatus( status.c_str() );

  cout << "\nEnter vector element: ";
  cin >> data;

} // end In


  /// Multiplies by negative unity vector
void  VectorVariable<1U>::Invert()
{
  data *= -1.;
}



void  VectorVariable<1U>::Out() const
{
  cout << "\nStatus: " << parseStatus( flag ) << "\t\t";
  cout << endl;
  cout << data << "\t\t";
  cout << endl;

} // end Out
 
 
 
} // end namespace csmp

