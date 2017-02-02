#ifndef CSMP_SCALAR_VARIABLE_H
#define CSMP_SCALAR_VARIABLE_H

#include "CSMP_number_types.h"
#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class VectorVariable;
template<size_t> class TensorVariable;

/**
 
@brief ScalarVariable combines a floating point value with a flag
that specifies its use during computations.

@author S.K. Matthai
@author Stephen G. Roberts
@author S. Geiger
@date 2001

Associates a floating point variable with a flag for the purpose of
special treatment in finite-element computations. Operators are 
defined so that ScalarVariable behaves essentially like a double,
but flags can be set, their values compared and the variable 
can be used together with VectorVariable and TensorVariable for the
representation of physical variables in CSMP computations.  
The ScalarVariable type is templatized to allow its use with variables
of different precision. The default use is with 'double64' variables
defined as 8-byte doubles.  

 
@section motivation Motivation

For the purpose of multiphysics computations, variable flags must be
combined with variable values and not with nodes, faces, or elements,
because several variables may be specified in such locations and their
flags may differ. Also range checking against the property database 
needs to be possible for physical variables.
 
 
@section applicability Applicability

Scalar- are used together with Vector- and TensorVariables to implement
CSP property computations.
 
 
@section implementation Implementation

The types of Scalar- are used together with Vector- and TensorVariables 
are linked with the enum VARIABLE_TYPE (SCALAR,VECTOR,TENSOR) for 
efficient discrimination (not depending on RTTI). They are stored in the 
MemoryManager object for each Model. The storage requirement
for the flag is 1-byte. Thus, the total storage required for a double 
scalar is 9 bytes (8 for the floating-point value and 1 for the flag).
 
 
@section examples Application Examples

This constructs a double ScalarVariable flagged to represent a NEUMANN
boundary condition with a value of 20. This value is then divided
by a double = 15.: 

@code
ScalarVariable<double>  my_value(NEUMANN,20.);

my_value /= 15.;
@endcode
 
*/
class ScalarVariable {
  public:
    ScalarVariable();
    ScalarVariable( VARIABLE_FLAG f, double64 val );
    ScalarVariable( const ScalarVariable& s );
    ScalarVariable( ScalarVariable&& s );
    ~ScalarVariable();
  
    ScalarVariable&  operator+=( double64 val );
    ScalarVariable&  operator-=( double64 val );
    ScalarVariable&  operator*=( double64 val );
    ScalarVariable&  operator/=( double64 val );

    ScalarVariable&  operator+=( const ScalarVariable& s );
    ScalarVariable&  operator-=( const ScalarVariable& s );
    ScalarVariable&  operator*=( const ScalarVariable& s );
    ScalarVariable&  operator/=( const ScalarVariable& s );
    
    ScalarVariable&  operator=( double64 val );
    ScalarVariable&  operator=( const ScalarVariable& s );

    bool             operator<(  double64 val ) const; 
    bool             operator>(  double64 val ) const; 
    bool             operator<=( double64 val ) const; 
    bool             operator>=( double64 val ) const; 
  
    bool             operator<(  const ScalarVariable& s ) const; 
    bool             operator>(  const ScalarVariable& s ) const; 
    bool             operator<=( const ScalarVariable& s ) const; 
    bool             operator>=( const ScalarVariable& s ) const;
  
    /// for storing scalars in associative containers with respective predicates
    bool             operator==( const ScalarVariable& s ) const; 
    bool             operator!=( const ScalarVariable& s ) const; 
  
    /// assignment as an lvalue
    double64&        operator()(void);
    double64         operator()(void) const;
  
    void             Component( size_t, double64 );
    double64         Component( size_t ) const;
    
    size_t           Components() const;
    double64         Value()   const;
    double64         Average() const;
  
    /// tests whether the variable value lies within the given bounds
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
  
    /// returns size of variable (=1 for scalar)
    size_t           Size() const;
  
    /// value assignment to scalar: cannot resize, but assigns user-defined or default value
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
  
    VARIABLE_FLAG&   Flag();
    VARIABLE_FLAG    Flag() const;
    void             Zero();
    void             Sqrt();
    /// natural logarithm of scalar
    void             Ln();
    /// base 10 logarithm
    void             Log10();
    void             Out() const;
    /// reading and writing scalars binary files
    bool             Out( FILE* fp ) const;
    bool             In( FILE* fp );

    friend class VectorVariable<2U>;
    friend class TensorVariable<2U>;
    friend class VectorVariable<3U>;
    friend class TensorVariable<3U>;

  private:
    VARIABLE_FLAG flag_;
    double64      data_;
 };

 ScalarVariable  operator+( const ScalarVariable&, const double64& );
 ScalarVariable  operator-( const ScalarVariable&, const double64& );
 ScalarVariable  operator*( const ScalarVariable&, const double64& );
 ScalarVariable  operator/( const ScalarVariable&, const double64& );
 ScalarVariable  operator+( const double64&, const ScalarVariable& );
 ScalarVariable  operator-( const double64&, const ScalarVariable& );
 ScalarVariable  operator*( const double64&, const ScalarVariable& );
 ScalarVariable  operator/( const double64&, const ScalarVariable& );
 ScalarVariable  operator+( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator-( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator*( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator/( const ScalarVariable&, const ScalarVariable& );

 template<size_t dim>
 VectorVariable<dim>  operator*( const ScalarVariable&, const VectorVariable<dim>& );

 template<size_t dim>
 TensorVariable<dim>  operator*( const ScalarVariable&, const TensorVariable<dim>& );


// *******************************************************************
//
//             INLINE FUNCTIONS
//
// *******************************************************************


/// for printing scalars using the standard streams cout, cerr, clog
std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& );

/// helper functions
const ScalarVariable&  makeScalar( VARIABLE_FLAG, double64 );


inline  double64& ScalarVariable::operator()(void) { return data_; }
inline  double64  ScalarVariable::operator()(void) const { return data_; }


inline  void      ScalarVariable::Component(size_t, double64 val) { data_ = val; }

inline  double64      ScalarVariable::Component(size_t) const { return data_; }


inline  size_t  ScalarVariable::Components() const { return 1U; }


inline double64       ScalarVariable::Value()   const  { return data_; }


inline double64       ScalarVariable::Average() const  { return data_; }


inline VARIABLE_FLAG&  ScalarVariable::Flag() { return flag_; }


inline VARIABLE_FLAG   ScalarVariable::Flag() const { return flag_; }


inline size_t ScalarVariable::Size() const  { return 1U;  }


inline void ScalarVariable::Resize( size_t, double64 newValue ) { data_ = newValue; }


inline void  ScalarVariable::Zero() { data_ = static_cast<double64>(0.);  }



inline void  ScalarVariable::Sqrt() 
  { 
     data_ = std::sqrt(data_);
  }



inline void  ScalarVariable::Ln()
  { 
     data_ = std::log(data_);
  }


inline void  ScalarVariable::Log10() 
  { 
     data_ = std::log10(data_);
  }



inline ScalarVariable::ScalarVariable() : flag_(ANY), data_(std::numeric_limits<double64>::quiet_NaN()) {}



inline ScalarVariable::ScalarVariable(  VARIABLE_FLAG f, double64 val )
 : flag_(f), data_(val)
 {
 }




inline ScalarVariable::ScalarVariable( const ScalarVariable& s ) 
 : flag_(s.flag_),
   data_(s.data_)
 { 
 }
 


inline ScalarVariable::ScalarVariable( ScalarVariable&& s )
 : flag_{s.flag_},
   data_{s.data_}
 { 
 }


inline ScalarVariable::~ScalarVariable() {}
 
 




inline ScalarVariable&  ScalarVariable::operator+=( const ScalarVariable& s )
 {
    data_ += s.data_;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator-=( const ScalarVariable& s )
 {
    data_ -= s.data_;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator*=( const ScalarVariable& s )
 {
    data_ *= s.data_;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator/=( const ScalarVariable& s )
 {
    data_ /= s.data_;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator+=( double64 val )
 {
    data_ += val;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator-=( double64 val )
 {
    data_ -= val;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator*=( double64 val )
 {
    data_ *= val;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator/=( double64 val )
 {
    data_ /= val;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator=( const ScalarVariable& s )
 {
    if ( &s == this ) return *this;
    flag_ = s.flag_;
    data_ = s.data_;
    return( *this );
 }
 
 


inline ScalarVariable&  ScalarVariable::operator=( double64 val )
 {
    data_ = val;
    return( *this );
 }




inline bool  ScalarVariable::operator<( double64 val ) const
 {
    return( data_ < val );
 }



inline bool  ScalarVariable::operator>( double64 val ) const
 {
    return( data_ > val );
 }



inline bool  ScalarVariable::operator<=( double64 val ) const
 {
    return( data_ <= val );
 }



inline bool  ScalarVariable::operator>=( double64 val ) const
 {
    return( data_ >= val );
 }


inline bool  ScalarVariable::operator<( const ScalarVariable& s ) const
 {
    return( s.data_ > data_ );
 }



inline bool  ScalarVariable::operator>( const ScalarVariable& s ) const
 {
    return( s.data_ < data_ );
 }



inline bool  ScalarVariable::operator<=( const ScalarVariable& s ) const
 {
    return( s.data_ >= data_ );
 }



inline bool  ScalarVariable::operator>=( const ScalarVariable& s ) const
 {
    return( s.data_ <= data_ );
 }
 


// keep for associative containers
inline bool  ScalarVariable::operator==( const ScalarVariable& s ) const
 {
    if ( s.flag_ != flag_ ) return false;
    return !(data_ > s.data_ and data_ < s.data_);
 }
 
  


// keep for associative containers
inline bool  ScalarVariable::operator!=( const ScalarVariable& s ) const
 {
    return !(*this == s);
 } 





inline bool  ScalarVariable::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data_ < vmin || data_ > vmax ) return false;
    return true;
 }



inline ScalarVariable  operator+( const ScalarVariable& l, const double64& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() + r ));
 }
 


inline ScalarVariable  operator-( const ScalarVariable& l, const double64& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() - r ));
 }
 


inline ScalarVariable  operator*( const ScalarVariable& l, const double64& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() * r ));
 }
 


inline ScalarVariable  operator/( const ScalarVariable& l, const double64& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() / r ));
 }



inline ScalarVariable  operator+( const double64& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( r.Flag(), l + r.Value() ));
 }
 


inline ScalarVariable  operator-( const double64& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( r.Flag(), l - r.Value() ));
 }
 


inline ScalarVariable  operator*( const double64& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( r.Flag(), l * r.Value() ));
 }
 


inline ScalarVariable  operator/( const double64& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( r.Flag(), l / r.Value() ));
 }



inline ScalarVariable  operator+( const ScalarVariable& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() + r.Value() ));
 }
 


inline ScalarVariable  operator-( const ScalarVariable& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() - r.Value() ));
 }
 


inline ScalarVariable  operator*( const ScalarVariable& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() * r.Value() ));
 }
 


inline ScalarVariable  operator/( const ScalarVariable& l, const ScalarVariable& r )
 {
    return std::move(ScalarVariable( l.Flag(), l.Value() / r.Value() ));
 }


template<size_t dim>
inline VectorVariable<dim>  operator*( const ScalarVariable& l, const VectorVariable<dim>& r )
{
    return r*l.Value();
}

template<size_t dim>
inline TensorVariable<dim>  operator*( const ScalarVariable& l, const TensorVariable<dim>& r )
{
    return r*l.Value();
}

// extensively tested fastest version that does not generate any temporaries
inline const ScalarVariable&  makeScalar( VARIABLE_FLAG flag, double64 val )
  {
     return std::move(ScalarVariable(flag,val));
  }

inline bool ScalarVariable::Out( FILE* fp ) const
  {
  fwrite( (void*)&flag_, sizeof(VARIABLE_FLAG), 1, fp);
  fwrite( (void*)&data_, sizeof(double64), 1, fp);
  return true;
  }

inline bool ScalarVariable::In( FILE* fp )
  {
  fread( (void*)&flag_, sizeof(VARIABLE_FLAG), 1, fp);
  fread( (void*)&data_, sizeof(double64), 1, fp);
  return true;
  }

} // csmp

#endif









