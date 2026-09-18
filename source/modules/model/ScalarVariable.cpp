// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ScalarVariable.h"
#include "compareFloats.h"

namespace csmp {

double& ScalarVariable::operator()( void ) noexcept { return data_; }
double  ScalarVariable::operator()( void ) const noexcept { return data_; }

VARIABLE_FLAG&  ScalarVariable::Flag() noexcept { return flag_; }


VARIABLE_FLAG   ScalarVariable::Flag() const noexcept { return flag_; }


void ScalarVariable::Resize( uint32_t, double newValue ) noexcept { data_ = newValue; }



ScalarVariable::ScalarVariable() noexcept : flag_( ANY ), data_( std::numeric_limits<double>::quiet_NaN() ) {}



ScalarVariable::ScalarVariable( VARIABLE_FLAG f, double val ) noexcept
  : flag_( f ), data_( val )
{
}


ScalarVariable&  ScalarVariable::operator+=( double val ) noexcept { data_ += val; return *this; }
ScalarVariable&  ScalarVariable::operator-=( double val ) noexcept { data_ -= val; return *this; }
ScalarVariable&  ScalarVariable::operator*=( double val ) noexcept { data_ *= val; return *this; }
ScalarVariable&  ScalarVariable::operator/=( double val ) noexcept { data_ /= val; return *this; }



ScalarVariable&  ScalarVariable::operator+=( const ScalarVariable& s ) noexcept
{
  data_ += s.data_;
  return(*this);
}

ScalarVariable&  ScalarVariable::operator-=( const ScalarVariable& s ) noexcept
{
  data_ -= s.data_;
  return(*this);
}

ScalarVariable&  ScalarVariable::operator*=( const ScalarVariable& s ) noexcept
{
  data_ *= s.data_;
  return(*this);
}

ScalarVariable&  ScalarVariable::operator/=( const ScalarVariable& s ) noexcept
{
  data_ /= s.data_;
  return(*this);
}




bool  ScalarVariable::operator<( const ScalarVariable& s ) const noexcept
{
  return(s.data_ > data_);
}



bool  ScalarVariable::operator>( const ScalarVariable& s ) const noexcept
{
  return(s.data_ < data_);
}



bool  ScalarVariable::operator<=( const ScalarVariable& s ) const noexcept
{
  return(s.data_ >= data_);
}



bool  ScalarVariable::operator>=( const ScalarVariable& s ) const noexcept
{
  return(s.data_ <= data_);
}



// keep for associative containers
bool  ScalarVariable::operator==( const ScalarVariable& s ) const noexcept
{
  if ( s.flag_ != flag_ ) return false;
  return essentiallyEqual(data_,s.data_);
}




// keep for associative containers
bool  ScalarVariable::operator!=( const ScalarVariable& s ) const noexcept
{
  return !(*this == s);
}





bool  ScalarVariable::IsWithinRange( double vmin, double vmax ) const noexcept
{
  if ( data_ < vmin || data_ > vmax ) return false;
  return true;
}



ScalarVariable  operator+( const ScalarVariable& l, const double& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() + r );
}



ScalarVariable  operator-( const ScalarVariable& l, const double& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() - r );
}



ScalarVariable  operator*( const ScalarVariable& l, const double& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() * r );
}



ScalarVariable  operator/( const ScalarVariable& l, const double& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() / r );
}



ScalarVariable  operator+( double l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( r.Flag(), l + r() );
}



ScalarVariable  operator-( double l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( r.Flag(), l - r() );
}



ScalarVariable  operator*( double l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( r.Flag(), l * r() );
}



ScalarVariable  operator/( double l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( r.Flag(), l / r() );
}



ScalarVariable  operator+( const ScalarVariable& l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() + r() );
}



ScalarVariable  operator-( const ScalarVariable& l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() - r() );
}



ScalarVariable  operator*( const ScalarVariable& l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() * r() );
}



ScalarVariable  operator/( const ScalarVariable& l, const ScalarVariable& r ) noexcept
{
  return ScalarVariable( l.Flag(), l() / r() );
}


template<uint32_t dim>
VectorVariable<dim>  operator*( const ScalarVariable& l, const VectorVariable<dim>& r ) noexcept
{
  return r * l();
}

template<uint32_t dim>
TensorVariable<dim>  operator*( const ScalarVariable& l, const TensorVariable<dim>& r ) noexcept
{
  return r * l();
}

// extensively tested fastest version that does not generate any temporaries
ScalarVariable  makeScalar( VARIABLE_FLAG flag, double val ) noexcept
{
  return ScalarVariable( flag, val );
}


bool ScalarVariable::Out( std::fstream& fp ) const
{
  const int32_t flag( flag_ );
  fp.write( (char*)&flag, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.write( (char*)&data_, sizeof( double ) );
  return true;
}


bool ScalarVariable::In( std::fstream& fp )
{
  fp.read( (char*)&flag_, sizeof( int32_t ) ); // VARIABLE_FLAG
  fp.read( (char*)&data_, sizeof( double ) );
  return true;
}


void  ScalarVariable::Out() const noexcept
{
  std::cout << "\nStatus: " << parseStatus( flag_ );
  if ( std::isnan( data_ ) )
    std::cout << ", value: NAN\n";
  else
    std::cout << ", value: " << data_ << std::endl;
}



std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& o ) noexcept
{
  stream << o() << " (" << parseStatus( o.Flag() ) << ")";
  return stream;
}


} // end namespace csmp
