#include "ScalarVariable.h"

namespace csmp {

double64& ScalarVariable::operator()(void) { return data_; }
double64  ScalarVariable::operator()(void) const { return data_; }

VARIABLE_FLAG&  ScalarVariable::Flag() { return flag_; }


VARIABLE_FLAG   ScalarVariable::Flag() const { return flag_; }


size_t ScalarVariable::Size() const  { return 1U;  }


void ScalarVariable::Resize( size_t, double64 newValue ) { data_ = newValue; }



ScalarVariable::ScalarVariable() : flag_(ANY), data_(std::numeric_limits<double64>::quiet_NaN()) {}



ScalarVariable::ScalarVariable(  VARIABLE_FLAG f, double64 val )
 : flag_(f), data_(val)
 {
 }




ScalarVariable::ScalarVariable( const ScalarVariable& s ) 
 : flag_(s.flag_),
   data_(s.data_)
 { 
 }
 


ScalarVariable::~ScalarVariable() {}
 
 




ScalarVariable&  ScalarVariable::operator+=( const ScalarVariable& s )
 {
    data_ += s.data_;
    return( *this );
 }
 
 


ScalarVariable&  ScalarVariable::operator-=( const ScalarVariable& s )
 {
    data_ -= s.data_;
    return( *this );
 }
 
 


ScalarVariable&  ScalarVariable::operator*=( const ScalarVariable& s )
 {
    data_ *= s.data_;
    return( *this );
 }
 
 


ScalarVariable&  ScalarVariable::operator/=( const ScalarVariable& s )
 {
    data_ /= s.data_;
    return( *this );
 }
 
 




ScalarVariable&  ScalarVariable::operator=( const ScalarVariable& s )
 {
    if ( &s == this ) return *this;
    flag_ = s.flag_;
    data_ = s.data_;
    return( *this );
 }
 

bool  ScalarVariable::operator<( const ScalarVariable& s ) const
 {
    return( s.data_ > data_ );
 }



bool  ScalarVariable::operator>( const ScalarVariable& s ) const
 {
    return( s.data_ < data_ );
 }



bool  ScalarVariable::operator<=( const ScalarVariable& s ) const
 {
    return( s.data_ >= data_ );
 }



bool  ScalarVariable::operator>=( const ScalarVariable& s ) const
 {
    return( s.data_ <= data_ );
 }
 


// keep for associative containers
bool  ScalarVariable::operator==( const ScalarVariable& s ) const
 {
    if ( s.flag_ != flag_ ) return false;
    return !(data_ > s.data_ and data_ < s.data_);
 }
 
  


// keep for associative containers
bool  ScalarVariable::operator!=( const ScalarVariable& s ) const
 {
    return !(*this == s);
 } 





bool  ScalarVariable::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data_ < vmin || data_ > vmax ) return false;
    return true;
 }



ScalarVariable  operator+( const ScalarVariable& l, const double64& r )
 {
    return ScalarVariable( l.Flag(), l() + r );
 }
 


ScalarVariable  operator-( const ScalarVariable& l, const double64& r )
 {
    return ScalarVariable( l.Flag(), l() - r );
 }
 


ScalarVariable  operator*( const ScalarVariable& l, const double64& r )
 {
    return ScalarVariable( l.Flag(), l() * r );
 }
 


ScalarVariable  operator/( const ScalarVariable& l, const double64& r )
 {
    return ScalarVariable( l.Flag(), l() / r );
 }



ScalarVariable  operator+( const double64& l, const ScalarVariable& r )
 {
    return ScalarVariable( r.Flag(), l + r() );
 }
 


ScalarVariable  operator-( const double64& l, const ScalarVariable& r )
 {
    return ScalarVariable( r.Flag(), l - r() );
 }
 


ScalarVariable  operator*( const double64& l, const ScalarVariable& r )
 {
    return ScalarVariable( r.Flag(), l * r() );
 }
 


ScalarVariable  operator/( const double64& l, const ScalarVariable& r )
 {
    return ScalarVariable( r.Flag(), l / r() );
 }



ScalarVariable  operator+( const ScalarVariable& l, const ScalarVariable& r )
 {
    return ScalarVariable( l.Flag(), l() + r() );
 }
 


ScalarVariable  operator-( const ScalarVariable& l, const ScalarVariable& r )
 {
    return ScalarVariable( l.Flag(), l() - r() );
 }
 


ScalarVariable  operator*( const ScalarVariable& l, const ScalarVariable& r )
 {
    return ScalarVariable( l.Flag(), l() * r() );
 }
 


ScalarVariable  operator/( const ScalarVariable& l, const ScalarVariable& r )
 {
    return ScalarVariable( l.Flag(), l() / r() );
 }


template<size_t dim>
VectorVariable<dim>  operator*( const ScalarVariable& l, const VectorVariable<dim>& r )
{
    return r * l();
}

template<size_t dim>
TensorVariable<dim>  operator*( const ScalarVariable& l, const TensorVariable<dim>& r )
{
    return r * l();
}

// extensively tested fastest version that does not generate any temporaries
ScalarVariable  makeScalar( VARIABLE_FLAG flag, double64 val )
  {
     return ScalarVariable(flag,val);
  }


bool ScalarVariable::Out( FILE* fp ) const
  {
  fwrite( (void*)&flag_, sizeof(VARIABLE_FLAG), 1, fp);
  fwrite( (void*)&data_, sizeof(double64), 1, fp);
  return true;
  }


bool ScalarVariable::In( FILE* fp )
  {
  fread( (void*)&flag_, sizeof(VARIABLE_FLAG), 1, fp);
  fread( (void*)&data_, sizeof(double64), 1, fp);
  return true;
  }


void  ScalarVariable::Out() const
 {
    std::cout <<"\nStatus: "<< parseStatus(flag_);
    if ( isnan(data_) )
      std::cout <<", value: NAN\n";
    else
      std::cout <<", value: " << data_ << std::endl;
 }

 

std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& o )
 {
    stream << o() <<" ("<< parseStatus(o.Flag()) <<")";
    return stream;
 }


} // end namespace csmp
