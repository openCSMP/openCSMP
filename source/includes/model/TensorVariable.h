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

Matrices can be input and output from streams, using the In() and
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

vc.Out(std::ostream& os);
@endcode

@todo (2-D) Const/Non const access operator convention not consistent with VectorVariable/ArrayVariable

*/
template<>
class TensorVariable<3U> {
  public:
    TensorVariable();
    TensorVariable( const TensorVariable& );
    TensorVariable( TensorVariable&& );
  
    /// creates isotropic diagonal tensor with diagonal elements equal to supplied value
    TensorVariable( VARIABLE_FLAG flag, double64 val );
  
    /// full initialisation where all diagonal elements have the same flag
    TensorVariable( VARIABLE_FLAG f,
                    double64 v11, double64 v12, double64 v13,
                    double64 v21, double64 v22, double64 v23,
                    double64 v31, double64 v32, double64 v33 );
  
    /// full initialisation
    TensorVariable( const VARIABLE_FLAG f11, const VARIABLE_FLAG f22, const VARIABLE_FLAG f33,
                    const double64  v11, const double64  v12, const double64  v13,
                    const double64  v21, const double64  v22, const double64  v23,
                    const double64  v31, const double64  v32, const double64  v33 );
                    
    ~TensorVariable();
  
    /// basic assignment
    TensorVariable&  operator=( const TensorVariable& );
    TensorVariable&  operator=( TensorVariable&& );
 
    /// read/write access to the elements of the tensor
    double64&        operator()( size_t i, size_t j );

    /// read-only access to the elements of the tensor
    const double64&  operator()( size_t i, size_t j ) const;
  
    /// alternative mutator of tensor elements 0..8 accessing them sequentially row by row
    void      Component( size_t, double64 );

    /// alternative accessor of tensor elements 0..8 accessing them sequentially row by row
    double64  Component( size_t i ) const;

    /// number of entries in tensor (9 in this 3D case
    size_t  Components() const;
    size_t  Size() const { return Components(); }

    /// assigns second argument to all elements of the tensor, first argument is not used; @todo remove
    void Resize( size_t, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );

    /// addition creating a temporary that is assigned to the lefthand Tensor
    TensorVariable   operator+( double64 val ) const;
    TensorVariable   operator-( double64 val ) const;
    TensorVariable   operator*( double64 val ) const;
    TensorVariable   operator/( double64 val ) const;
    
    TensorVariable&  operator+=( double64 val );
    TensorVariable&  operator-=( double64 val );
    TensorVariable&  operator*=( double64 val );
    TensorVariable&  operator/=( double64 val );

    TensorVariable&  operator+=( const ScalarVariable& );
    TensorVariable&  operator-=( const ScalarVariable& );
    TensorVariable&  operator*=( const ScalarVariable& );
    TensorVariable&  operator/=( const ScalarVariable& );
    TensorVariable   operator+(  const TensorVariable& ) const;
    TensorVariable   operator-(  const TensorVariable& ) const;
    TensorVariable   operator*(  const TensorVariable& ) const;
    
    /// (A x) matrix-vector multiplication -> vector (v treated as column vector)
    VectorVariable<3U>  operator*( const VectorVariable<3U>& ) const;
    
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
    TensorVariable&  operator=( double64 val );

    /// sets all ij values to sc (value and flags of diagonal elements)
    TensorVariable&  operator=( const ScalarVariable& );
  
    /// assigns vector to diagonal elements of tensor, off-diagonal elements are set to zero
    TensorVariable&  operator=( const VectorVariable<3U>& );
    
    /// element-by-element comparison of ij values
    bool             operator==( const TensorVariable& ) const;
    
    /// element-by-element comparison of ij values
    bool             operator!=( const TensorVariable& ) const;
    
    /// comparison of the determinants of the tensors
    bool             operator<( const TensorVariable& ) const; 
    
    double64          Value( size_t i, size_t j ) const;
  
    /// averages all elements of tensor (may make sense only for diagonalised tensor)
    double64          Average() const;
  
    /// compares the individual elements of the tensor with the ranges specified in PropertyDatabase file
    bool              IsWithinRange( double64 vmin, double64 vmax ) const;
  
    /// returns the flag of the diagonal tensor element of choice
    VARIABLE_FLAG     Flag( size_t i=0 ) const;
  
    /// returns smallest element in tensor (this is zero if the tensor is diagonal)
    double64          MinElement() const;

    /// returns largest element in tensor
    double64          MaxElement() const;
    double64          Determinant() const;
  
    /// returns the sum of the diagonal values of the tensor
    double64          Trace() const;
  
    /// returns the conjugate transpose of the tensor (commonly denoted M^(H))
    TensorVariable    Adjoint()     const;
  
    /// returns a tensor the product of which with the original tensor gives the identity matrix
    TensorVariable    Inverse()     const;
  
    /// returns the transposed of the tensor into a new tensor variable
    TensorVariable    Transposed()  const;
  
    /// assuming that the tensor is symmetric and positive definite, method returns its sorted Eigen values (largest to smallest)
		bool              EigenValuesPositiveDefiniteSymmetricMatrix( double64& eigenValue0, double64& eigenValue1, double64& eigenValue2 ) const;
  
    /// alternative Eigen decomposition that should also work for non-symmetric matrices
    bool              EigenNonSymmetric( VectorVariable<3U>& eigenVals, TensorVariable<3U>& eigenVecs ) const;

    /// assuming that the tensor is symmetric and positive definite, its Eigenvalues are returned into the supplied VectorVariable
    bool              EigenValues( VectorVariable<3U>& vecEigenvalues ) const;

    /// assuming that the tensor is symmetric and positive definite, its Eigenvalues are returned into the supplied STL vector
    bool              EigenValues( std::vector<double64>& vecEigenvalues ) const;
  
    /// assuming that the tensor is symmetric and positive definite, method returns Eigenvalues and vectors that can be normalised to 1
    bool              Eigen( VectorVariable<3U>& vvEigenvalues, TensorVariable<3U>& tvEigenvectors, bool bNormalize ) const;

    /// assigns diagonal values to tensor (off-diagonal elements are not touched)
    void              DiagonalValues( double64 f_00, double64 f_11, double64 f_22 );

    /// assigns diagonal values to tensor from STL vector
    void              DiagonalValues( const std::vector<double64>& );

    /// assigns diagonal values to tensor
    void              DiagonalValues( const VectorVariable<3U>& );
  
    /// assigns vector variable to row i of the tensor
    void              AssignToRow( size_t i, VectorVariable<3U>& );

    /// assigns vector variable to column j of the tensor
    void              AssignToColumn( size_t j, VectorVariable<3U>& );
  
    /// returns row iRow into the argument VectorVariable
    VectorVariable<3U> Row( size_t iRow ) const;

    /// returns column iCol into the argument VectorVariable
    VectorVariable<3U> Column( size_t iCol ) const;

    /// accessor/mutator of the flags of the diagonal elements of the tensor
    VARIABLE_FLAG&    Flag( size_t i=0 );
  
    /// replaces all elements by 0.
    void              Zero();
  
    /// converts tensor into identity matrix
    void              Identity();
  
    /// shorthands for efficient tensor modification for display
    void              Fabs();
    void              Sqrt( bool from_absolute_value=false );
    void              Ln( bool from_absolute_value=false );
    void              Log10( bool from_absolute_value=false );

    /// prompts user to initialise the tensor from the console
    void              In();
  
    /// prints the tensor to a stream
    void              Out() const { Out(std::cout); }
    void              Out(std::ostream& os) const;
  
    /// reads the tensor from the supplied input file
    bool              In( FILE* fp );
  
    /// writes the tensor to the supplied output file
    bool              Out( FILE* fp ) const;

  private:
    std::array<VARIABLE_FLAG,3U>            flag;
    std::array<std::array<double64,3U>,3U>  data;
};


/// (x^T A)^T = A^T x  (vector - matrix multiplication -> vector(A.cols))
VectorVariable<3U>  operator*( const VectorVariable<3U>& vc, const TensorVariable<3U>& ts );

/// ostream operator for exporting the tensor to cout or file streams
template<size_t dim>
std::ostream&  operator<<( std::ostream& stream, const TensorVariable<dim>& o );


/// fastest way to insert a tensor into an STL container
const TensorVariable<2U>& makeTensor( VARIABLE_FLAG, VARIABLE_FLAG,
                                      double64, double64,
                                      double64, double64 );

const TensorVariable<3U>& makeTensor( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG,
                                      double64, double64, double64,
                                      double64, double64, double64,
                                      double64, double64, double64 );



// ******************************************************************************************
//
//            INLINE METHODS START HERE
//
// ******************************************************************************************

inline TensorVariable<3U>::~TensorVariable() {}


inline double64& TensorVariable<3U>::operator()( size_t i, size_t j ) 
 {
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }




inline const double64& TensorVariable<3U>::operator()( size_t i, size_t j ) const
 {
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }



inline size_t  TensorVariable<3U>::Components() const { return 9U; } 




inline void TensorVariable<3U>::Component( size_t i, double64 val )
 { 
    assert( i < Components() );
    // row by row
    if ( i == 0U )      data[0U][0U] = val;
    else if ( i == 1U ) data[0U][1U] = val;
    else if ( i == 2U ) data[0U][2U] = val;

    else if ( i == 3U ) data[1U][0U] = val;
    else if ( i == 4U ) data[1U][1U] = val;
    else if ( i == 5U ) data[1U][2U] = val;

    else if ( i == 6U ) data[2U][0U] = val;
    else if ( i == 7U ) data[2U][1U] = val;

    else data[2U][2U] = val; // remaining case
 }




inline double64 TensorVariable<3U>::Component( size_t i ) const
 { 
    assert( i < Components() );
    // row by row
    if ( i == 0U ) return data[0U][0U];
    if ( i == 1U ) return data[0U][1U];
    if ( i == 2U ) return data[0U][2U];

    if ( i == 3U ) return data[1U][0U];
    if ( i == 4U ) return data[1U][1U];
    if ( i == 5U ) return data[1U][2U];

    if ( i == 6U ) return data[2U][0U];
    if ( i == 7U ) return data[2U][1U];

    return data[2U][2U]; // remaining case
 }




inline double64 TensorVariable<3U>::Value( size_t i, size_t j ) const
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::Value(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::Value(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j];
 }
 
 

inline VARIABLE_FLAG& TensorVariable<3U>::Flag( size_t i )      
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }



inline VARIABLE_FLAG  TensorVariable<3U>::Flag( size_t i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nTensorVariable<3U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }




inline bool TensorVariable<3U>::EigenValues( VectorVariable<3U>& Ev ) const
 {
		return EigenValuesPositiveDefiniteSymmetricMatrix( Ev(0), Ev(1), Ev(2) );
	}


inline bool TensorVariable<3U>::EigenValues( std::vector<double64>& Ev ) const
	{
		return EigenValuesPositiveDefiniteSymmetricMatrix( Ev[0], Ev[1], Ev[2] );
	}




/**
    Initialiser for tensor variable, index is not used.
*/
inline void TensorVariable<3U>::Resize( size_t, double64 newValue )
  {
    data[0][0] = newValue;
    data[0][1] = newValue;
    data[1][0] = newValue;
    data[1][1] = newValue;
    data[0][2] = newValue;
    data[1][2] = newValue;
    data[2][0] = newValue;
    data[2][1] = newValue;
    data[2][2] = newValue;
  }



inline double64 TensorVariable<3U>::Trace() const
 {
    return data[0][0] + data[1][1] + data[2][2];   
 }



inline void TensorVariable<3U>::AssignToRow( size_t iRow, VectorVariable<3U>& vc )
{
  if ( iRow == 0U )
      flag[0U] = vc.Flag(0U);
  else if ( iRow == 1U )
      flag[1U] = vc.Flag(1U);
  else
      flag[2U] = vc.Flag(2U);

	data[iRow][0U] = vc[0U];
	data[iRow][1U] = vc[1U];
	data[iRow][2U] = vc[2U];
}


inline void TensorVariable<3U>::AssignToColumn( size_t iCol, VectorVariable<3U>& vc )
{
  if ( iCol == 0U )
      flag[0U] = vc.Flag(0U);
  else if ( iCol == 1U )
      flag[1U] = vc.Flag(1U);
  else
      flag[2U] = vc.Flag(2U);

	data[0U][iCol] = vc[0U];
	data[1U][iCol] = vc[1U];
	data[2U][iCol] = vc[2U];
}




inline VectorVariable<3U> TensorVariable<3U>::Row( size_t iRow ) const
{
	return VectorVariable<3U>( flag[iRow], flag[iRow], flag[iRow], 
	                           data[iRow][0U], data[iRow][1U], data[iRow][2U]);
}


inline VectorVariable<3U> TensorVariable<3U>::Column( size_t iCol ) const
{
	return VectorVariable<3U>( flag[iCol], flag[iCol], flag[iCol], 
	                           data[0U][iCol], data[1U][iCol], data[2U][iCol] );
}
 

inline bool  TensorVariable<3U>::operator==( const TensorVariable<3U>& ts ) const
 {
   return ( data == ts.data && flag == ts.flag );
 }
  

inline bool  TensorVariable<3U>::operator!=( const TensorVariable<3U>& t ) const
 {
     return !(*this == t);
 } 



inline bool  TensorVariable<3U>::operator<( const TensorVariable<3U>& t ) const
 {
     return (this->Determinant() < t.Determinant());
 } 

inline  bool TensorVariable<3U>::Out( FILE* fp ) const
  {
     fwrite( (void*)this, sizeof(TensorVariable<3U>), 1, fp );
     return true;
  }

inline  bool TensorVariable<3U>::In( FILE* fp )
  {
     fread( (void*)this, sizeof(TensorVariable<3U>), 1, fp );
     return true;
  }
  
  
/// fastest way to insert a tensor into an STL container; tensor only has flags for diagonal elements
inline const TensorVariable<2U>& makeTensor( VARIABLE_FLAG f1, VARIABLE_FLAG f2,
                                             double64 v11, double64 v12,
                                             double64 v21, double64 v22 )
 {
     return std::move(TensorVariable<2U>(f1,f2,v11,v12,v21,v22) );
 }
 

/// fastest way to insert a tensor into an STL container; tensor only has flags for diagonal elements
inline const TensorVariable<3U>& makeTensor( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3,
                                             double64 v11, double64 v12, double64 v13,
                                             double64 v21, double64 v22, double64 v23,
                                             double64 v31, double64 v32, double64 v33 )
 {
     return std::move(TensorVariable<3U>(f1,f2,f3,v11,v12,v13,v21,v22,v23,v31,v32,v33) );
 }
  
  

} // csmp 
 
#endif

