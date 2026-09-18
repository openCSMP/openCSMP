// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_SCALAR_VARIABLE_H
#define CSMP_SCALAR_VARIABLE_H

#include "CSMP_definitions.h"

namespace csmp {

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
of different precision. The default use is with 'double' variables
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

  ScalarVariable() noexcept;
  ScalarVariable( VARIABLE_FLAG f, double val ) noexcept;

  ScalarVariable&  operator+=( double ) noexcept;
  ScalarVariable&  operator-=( double ) noexcept;
  ScalarVariable&  operator*=( double ) noexcept;
  ScalarVariable&  operator/=( double ) noexcept;

  ScalarVariable&  operator+=( const ScalarVariable& ) noexcept;
  ScalarVariable&  operator-=( const ScalarVariable& ) noexcept;
  ScalarVariable&  operator*=( const ScalarVariable& ) noexcept;
  ScalarVariable&  operator/=( const ScalarVariable& ) noexcept;

  // for conformance with the interfaces of the other CSMP variables
  ScalarVariable&  operator=( double val ) noexcept { data_ = val; return *this; }

  /// comparitor that is used by less<> predicate in STL
  bool             operator<( const ScalarVariable& ) const noexcept;
  /// comparitor that is used by greater than<> predicate in STL
  bool             operator>( const ScalarVariable& ) const noexcept;
  bool             operator<=( const ScalarVariable& ) const noexcept;
  bool             operator>=( const ScalarVariable& ) const noexcept;

  /// for storing scalars in associative containers with respective predicates
  bool             operator==( const ScalarVariable& ) const noexcept;
  bool             operator!=( const ScalarVariable& ) const noexcept;

  /// assignment as an lvalue
  double&          operator()( void ) noexcept;
  double           operator()( void ) const noexcept;

  /// tests whether the variable value lies within the given bounds
  bool             IsWithinRange( double vmin, double vmax ) const noexcept;
  
  /// tests whether variable contains NaN value(s)
  bool             Has_NaN_Values() const noexcept { return std::isnan(data_); }

  /// universal way of assigning values to all CSMP variable types
  void             Component( uint32_t, double val ) noexcept { data_ = val; }

  /// universal accessor of CSMP variable values which works for all variable types
  [[nodiscard]] double Component( uint32_t ) const  noexcept{ return data_; }

  /// returns size = number of components of the variable (=1 for scalar)
  static constexpr uint32_t Size() noexcept { return 1u; };

  /// value assignment to scalar: cannot resize, but assigns user-defined or default value
  void             Resize( uint32_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() ) noexcept;

  /// assigment: status of variable which determines how it is used in computations
  VARIABLE_FLAG&   Flag() noexcept;

  /// accessor: status of variable which determines how it is used in computations
  VARIABLE_FLAG    Flag() const noexcept;

  /// prints flag/value pair to screen
  void             Out() const noexcept;

  /// reading and writing of scalar variables to binary files
  bool             Out( std::fstream& fp ) const;
  bool             In( std::fstream& fp );

private:
  VARIABLE_FLAG flag_;
  double      data_;
};

/// creates scalar and returns; use for inserting scalars into functions
ScalarVariable  makeScalar( VARIABLE_FLAG, double ) noexcept;

ScalarVariable  operator+( const ScalarVariable&, const ScalarVariable& ) noexcept;
ScalarVariable  operator-( const ScalarVariable&, const ScalarVariable& ) noexcept;
ScalarVariable  operator*( const ScalarVariable&, const ScalarVariable& ) noexcept;
ScalarVariable  operator/( const ScalarVariable&, const ScalarVariable& ) noexcept;

ScalarVariable  operator+( const ScalarVariable&, const double& ) noexcept;
ScalarVariable  operator-( const ScalarVariable&, const double& ) noexcept;
ScalarVariable  operator*( const ScalarVariable&, const double& ) noexcept;
ScalarVariable  operator/( const ScalarVariable&, const double& ) noexcept;

ScalarVariable  operator+( double, const ScalarVariable& ) noexcept;
ScalarVariable  operator-( double, const ScalarVariable& ) noexcept;
ScalarVariable  operator*( double, const ScalarVariable& ) noexcept;
ScalarVariable  operator/( double, const ScalarVariable& ) noexcept;

/// multiplies each element of vector variable with scalar
template<uint32_t dim>
VectorVariable<dim>  operator*( const ScalarVariable&, const VectorVariable<dim>& ) noexcept;

/// multiplies each element of tensor variable with scalar
template<uint32_t dim>
TensorVariable<dim>  operator*( const ScalarVariable&, const TensorVariable<dim>& ) noexcept;

/// for printing scalars using the standard streams cout, cerr, clog
std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& ) noexcept;


} // csmp

#endif
