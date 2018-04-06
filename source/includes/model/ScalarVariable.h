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
    static constexpr VARIABLE_TYPE VariableType = SCALAR;

    ScalarVariable();
    ScalarVariable( VARIABLE_FLAG f, double64 val );
    ScalarVariable( const ScalarVariable& );
    ScalarVariable( ScalarVariable&& ) = default;
    ~ScalarVariable();
  
    ScalarVariable&  operator+=( const ScalarVariable& );
    ScalarVariable&  operator-=( const ScalarVariable& );
    ScalarVariable&  operator*=( const ScalarVariable& );
    ScalarVariable&  operator/=( const ScalarVariable& );
  
    // for conformance with the interfaces of the other CSMP variables
    ScalarVariable&  operator=( double64 val )  { data_ = val; return *this; }
    ScalarVariable&  operator+=( double64 val ) { data_ += val; return *this; }
    ScalarVariable&  operator-=( double64 val ) { data_ -= val; return *this; }
    ScalarVariable&  operator*=( double64 val ) { data_ *= val; return *this; }
    ScalarVariable&  operator/=( double64 val ) { data_ /= val; return *this; }

    /// assignment operator
    ScalarVariable&  operator=( const ScalarVariable& );
  
    /// move assignment
    ScalarVariable&  operator=( ScalarVariable&& ) = default;

    /// comparitor that is used by less<> predicate in STL
    bool             operator<(  const ScalarVariable& ) const;
    /// comparitor that is used by greater than<> predicate in STL
    bool             operator>(  const ScalarVariable& ) const;
    bool             operator<=( const ScalarVariable& ) const;
    bool             operator>=( const ScalarVariable& ) const;
  
    /// for storing scalars in associative containers with respective predicates
    bool             operator==( const ScalarVariable& ) const;
    bool             operator!=( const ScalarVariable& ) const;
  
    /// assignment as an lvalue
    double64&        operator()(void);
    double64         operator()(void) const;
  
    /// tests whether the variable value lies within the given bounds
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;

    /// universal way of assigning values to all CSMP variable types
    void             Component( size_t, double64 val ) { data_ = val; }

    /// universal accessor of CSMP variable values which works for all variable types
    double64         Component( size_t ) const { return data_; }
  
    /// returns size = number of components of the variable (=1 for scalar)
    size_t           Size() const;
  
    /// value assignment to scalar: cannot resize, but assigns user-defined or default value
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
  
    /// assigment: status of variable which determines how it is used in computations
    VARIABLE_FLAG&   Flag();

    /// accessor: status of variable which determines how it is used in computations
    VARIABLE_FLAG    Flag() const;

    /// prints flag/value pair to screen
    void             Out() const;
  
    /// reading and writing of scalar variables to binary files
    bool             Out( FILE* fp ) const;
    bool             In( FILE* fp );

  private:
    VARIABLE_FLAG flag_;
    double64      data_;
 };

 /// creates scalar and returns; use for inserting scalars into functions
 ScalarVariable  makeScalar( VARIABLE_FLAG, double64 );

 ScalarVariable  operator+( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator-( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator*( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator/( const ScalarVariable&, const ScalarVariable& );
 ScalarVariable  operator+( const ScalarVariable&, const double64& );
 ScalarVariable  operator-( const ScalarVariable&, const double64& );
 ScalarVariable  operator*( const ScalarVariable&, const double64& );
 ScalarVariable  operator/( const ScalarVariable&, const double64& );
 ScalarVariable  operator+( const double64&, const ScalarVariable& );
 ScalarVariable  operator-( const double64&, const ScalarVariable& );
 ScalarVariable  operator*( const double64&, const ScalarVariable& );
 ScalarVariable  operator/( const double64&, const ScalarVariable& );

 /// multiplies each element of vector variable with scalar
 template<size_t dim>
 VectorVariable<dim>  operator*( const ScalarVariable&, const VectorVariable<dim>& );

 /// multiplies each element of tensor variable with scalar
 template<size_t dim>
 TensorVariable<dim>  operator*( const ScalarVariable&, const TensorVariable<dim>& );
 
  /// for printing scalars using the standard streams cout, cerr, clog
  std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& );


} // csmp

#endif









