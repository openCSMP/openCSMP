#ifndef CSMP_TENSOR_VARIABLE_H
#define CSMP_TENSOR_VARIABLE_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable1.h"
#include "TensorVariable2.h"

namespace csmp {

/**
@brief full 3D specialisation of the TensorVariable class template

@author S.K. Matthai
@author Stephen G. Roberts
@date 2001

@section motivation Motivation

Template specialization implementing optimized tensors in CSP. Depending
on the spatial dimension of the model, tensor variables are either
2 x 2 or 3 x 3 matrices. Each value is associated with a flag which
can be set to indicate to finite-element computations that the
corresponding degree of freedom is enabled or restricted.


@section functionality Functionality

Apart from the standard constructors with a default where all elements
are initialized to NaN (not-a-number) and the flags are set to plain,
assignment is defined to allow the storage of tensors in STL containers.
Ordering using STL's less() is facilitated through an operator<() that
compares the determinants of matrices. All costly operations like
this are unrolled.

TensorVariable implements typical accessors, mutators for tensor elements
and operations on entire tensors. To get at element 1,2 of 3D tensor, ts,
use ts(1,2), for instance.

Value-by-value operations such as matrix addition, subtraction, division
are implemented by the operators +, -, /. operator*() assumes different
meanings dependent on the type of its arguments. Thus, tensors can be
multiplied with scalars, vectors, and tensors to produce the standard
results expected for matrix algebra. Other matrix algebraic
operations include finding the Transposed(), Determinant(), Inverse(),
Identity(), and Adjoint() matrices.

The *= operator will always produce a temporary matrix. Therefore,
* should be preferred if the result is to be assigned to a new variable.

Assignent of a vector to a tensor is done by placing the vector in the
tensors diagonal and zeroing its off-diagonal elements.

Matrices can be input and output via the console, using the In() and
Out() member functions.


@section motivation Motivation

Part of a consistent framework for physical variables for which there
exist overloaded operators so that vector-matrix multiplications and
other operations can be written in a straight forward way.


@section implementation Implementation

Using the template specialization mechanism of C++. Specializations
exist only for the 2 and 3D cases. These are declared in this (3D)
and the header file 'TensorVariable2.h' (2D), respectively.


@section examples Application Examples

To diagonalize a tensor variable, one could write:

@code
VectorVariable<double,2> vc(PLAIN,PLAIN,1,1);
TensorVariable<double,2> ts;

ts.Identity();
vc = ts * vc;

vc.Out();
@endcode

@todo (2-D) Const/Non const access operator convention not consistent with VectorVariable/ArrayVariable
*/
template<>
class TensorVariable<3U> {
public:
  static constexpr VARIABLE_TYPE VariableType = TENSOR;

  TensorVariable();
  TensorVariable( const TensorVariable& );
  TensorVariable( TensorVariable&& ) = default;

  /// creates isotropic diagonal tensor with diagonal elements equal to supplied value
  TensorVariable( VARIABLE_FLAG flag, double val );

  /// full initialisation where all diagonal elements have the same flag
  TensorVariable( VARIABLE_FLAG f,
                  double v11, double v12, double v13,
                  double v21, double v22, double v23,
                  double v31, double v32, double v33 );

  /// full initialisation
  TensorVariable( const VARIABLE_FLAG f11, const VARIABLE_FLAG f22, const VARIABLE_FLAG f33,
                  const double  v11, const double  v12, const double  v13,
                  const double  v21, const double  v22, const double  v23,
                  const double  v31, const double  v32, const double  v33 );

  ~TensorVariable();

  /// basic assignment
  TensorVariable&  operator=( const TensorVariable& );
  TensorVariable&  operator=( TensorVariable&& ) = default;

  /// read/write access to the elements of the tensor
  double&        operator()( uint32_t i, uint32_t j );

  /// read-only access to the elements of the tensor
  const double&  operator()( uint32_t i, uint32_t j ) const;

  /// alternative mutator of tensor elements 0..8 accessing them sequentially row by row
  void      Component( uint32_t, double );

  /// alternative accessor of tensor elements 0..8 accessing them sequentially row by row
  double    Component( uint32_t i ) const;

  /// number of entries in tensor (dim x dim = 9 in this 3D case)
  uint32_t  Size() const { return 9U; }

  /// assigns second argument to all elements of the tensor, first argument is not used; @todo remove
  void Resize( uint32_t, double newValue = std::numeric_limits<double>::quiet_NaN() );

  TensorVariable   operator+( double val ) const;
  TensorVariable   operator-( double val ) const;
  TensorVariable   operator*( double val ) const;
  TensorVariable   operator/( double val ) const;

  TensorVariable&  operator+=( double val );
  TensorVariable&  operator-=( double val );
  TensorVariable&  operator*=( double val );
  TensorVariable&  operator/=( double val );

  TensorVariable&  operator+=( const ScalarVariable& );
  TensorVariable&  operator-=( const ScalarVariable& );
  TensorVariable&  operator*=( const ScalarVariable& );
  TensorVariable&  operator/=( const ScalarVariable& );
  TensorVariable   operator+( const TensorVariable& ) const;
  TensorVariable   operator-( const TensorVariable& ) const;
  TensorVariable   operator*( const TensorVariable& ) const;

  /// (A x) matrix-vector multiplication -> vector (v treated as column vector)
  VectorVariable<3U>  operator*( const VectorVariable<3U>& ) const;

  /// (A x) matrix-vector multiplication -> vector (v treated as column vector)
  Point<3U>  operator*( const Point<3U>& ) const;

  /// value by value division of the elements of the tensor with another one
  TensorVariable   operator/( const TensorVariable& ) const;

  /// value by value addition of the elements of the tensors
  TensorVariable&  operator+=( const TensorVariable& );

  /// value by value subtraction from the lh-tensor by the rh tensor
  TensorVariable&  operator-=( const TensorVariable& );

  /// value by value division of the elements of the lh-tensor by those of the righthand tensor
  TensorVariable&  operator/=( const TensorVariable& );

  /// matrix multiplication -> matrix (A.rows, B.cols); creates temporary tensor, avoid this using 'operator*()'
  TensorVariable&  operator*=( const TensorVariable& );

  /// sets all ij values to val or sc
  TensorVariable&  operator=( double val );

  /// sets all ii values to sc (value and flags of diagonal elements)
  TensorVariable&  operator=( const ScalarVariable& );

  /// assigns vector to diagonal elements of tensor, off-diagonal elements are set to zero
  TensorVariable&  operator=( const VectorVariable<3U>& );

  /// element-by-element comparison of ij values
  bool             operator==( const TensorVariable& ) const;

  /// element-by-element comparison of ij values
  bool             operator!=( const TensorVariable& ) const;

  /// comparison of the determinants of the tensors
  bool             operator<( const TensorVariable& ) const;

  /// compares the individual elements of the tensor with the ranges specified in PropertyDatabase file
  bool             IsWithinRange( double vmin, double vmax ) const;
  
  bool             Has_NaN_Values() const;


  /// returns the flag of the diagonal tensor element of choice
  VARIABLE_FLAG    Flag( uint32_t i = 0 ) const;

  /// returns smallest element in tensor (this is zero if the tensor is diagonal)
  double           MinElement() const;

  /// returns largest element in tensor
  double           MaxElement() const;

  double           Determinant() const;

  /// returns the sum of the diagonal values of the tensor
  double           Trace() const;

  /// returns the conjugate transpose of the tensor (commonly denoted M^(H))
  TensorVariable   Adjoint() const;

  /// returns a tensor the product of which with the original tensor gives the identity matrix
  TensorVariable   Inverse() const;

  /// returns the transposed of the tensor into a new tensor variable
  TensorVariable   Transposed()  const;

  /// assuming that the tensor is symmetric and positive definite, method returns its sorted Eigen values (largest to smallest)
  bool             EigenValuesPositiveDefiniteSymmetricMatrix( double& eigenValue0,
                                                               double& eigenValue1,
                                                               double& eigenValue2 ) const;

  /// alternative Eigen decomposition that should also work for non-symmetric matrices
  bool              EigenNonSymmetric( VectorVariable<3U>& eigenVals, TensorVariable<3U>& eigenVecs ) const;

  /// assuming that the tensor is symmetric and positive definite, its Eigenvalues are returned into the supplied VectorVariable
  bool              EigenValues( VectorVariable<3U>& vecEigenvalues ) const;

  /// assuming that the tensor is symmetric and positive definite, its Eigenvalues are returned into the supplied STL vector
  bool              EigenValues( std::vector<double>& vecEigenvalues ) const;

  /// assuming that the tensor is symmetric and positive definite, method returns Eigenvalues and vectors that can be normalised to 1
  bool              Eigen( VectorVariable<3U>& vvEigenvalues, TensorVariable<3U>& tvEigenvectors, bool bNormalize ) const;

  /// assigns diagonal values to tensor (off-diagonal elements are not touched)
  void              DiagonalValues( double f_00, double f_11, double f_22 );

  /// assigns diagonal values to tensor from STL vector
  void              DiagonalValues( const std::vector<double>& );

  /// assigns diagonal values to tensor
  void              DiagonalValues( const VectorVariable<3U>& );

  /// assigns vector variable to row i of the tensor
  void              AssignToRow( uint32_t i, VectorVariable<3U>& );

  /// assigns vector variable to column j of the tensor
  void              AssignToColumn( uint32_t j, VectorVariable<3U>& );

  /// returns row iRow into the argument VectorVariable
  VectorVariable<3U> Row( uint32_t iRow ) const;

  /// returns column iCol into the argument VectorVariable
  VectorVariable<3U> Column( uint32_t iCol ) const;

  /// accessor/mutator of the flags of the diagonal elements of the tensor
  VARIABLE_FLAG&    Flag( uint32_t i = 0 );

  /// converts tensor into identity matrix
  void              Identity();

  /// prompts user to initialise the tensor from the command line
  void              In();

  /// prints the tensor to the command line
  void              Out() const;

  /// reads the tensor from the supplied input file
  bool              In( std::fstream& fp );

  /// writes the tensor to the supplied output file
  bool              Out( std::fstream& fp ) const;

private:
  std::array<VARIABLE_FLAG, 3U>           flag;
  std::array<std::array<double, 3U>, 3U>  data;
};


/// (x^T A)^T = A^T x  (vector - matrix multiplication -> vector(A.cols))
VectorVariable<3U>  operator*( const VectorVariable<3U>& vc, const TensorVariable<3U>& ts );

/// (x^T A)^T = A^T x  (vector - matrix multiplication -> vector(A.cols))
Point<3U>  operator*( const Point<3U>& vc, const TensorVariable<3U>& ts );


/// ostream operator for exporting the tensor to cout or file streams
template<uint32_t dim>
std::ostream&  operator<<( std::ostream& stream, const TensorVariable<dim>& o );


/// fastest way to insert a tensor into an STL container
TensorVariable<2U> makeTensor( VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double,
                               double, double );

TensorVariable<3U> makeTensor( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG,
                               double, double, double,
                               double, double, double,
                               double, double, double );


} // csmp 

#endif

