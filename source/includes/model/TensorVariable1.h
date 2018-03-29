#ifndef CSMP_TENSOR_VARIABLE1_H
#define CSMP_TENSOR_VARIABLE1_H

#include "ScalarVariable.h"
#include "VectorVariable1.h"

namespace csmp {

/// @brief 1D full specialization of tensor-type variable template (other tensors are declared in files #1, #2)
template<>
class TensorVariable<1U> {
  public:
    static constexpr VARIABLE_TYPE VariableType = TENSOR;

    TensorVariable();
    TensorVariable( const TensorVariable& );
    TensorVariable( VARIABLE_FLAG, double64 );
    TensorVariable( TensorVariable&& ) = default;
                                                      
    ~TensorVariable();
    
    double64&        operator()( size_t i, size_t j );
    const double64&  operator()( size_t i, size_t j ) const;
    void             Component( size_t, double64 );
    double64         Component( size_t i ) const;

    TensorVariable   operator+( double64 val ) const;
    TensorVariable   operator-( double64 val ) const;
    TensorVariable   operator*( double64 val ) const;
    TensorVariable   operator/( double64 val ) const;
    
    TensorVariable&  operator+=( double64 val );
    TensorVariable&  operator-=( double64 val );
    TensorVariable&  operator*=( double64 val );
    TensorVariable&  operator/=( double64 val );

    TensorVariable&  operator+=( const ScalarVariable& sc );
    TensorVariable&  operator-=( const ScalarVariable& sc );
    TensorVariable&  operator*=( const ScalarVariable& sc );
    TensorVariable&  operator/=( const ScalarVariable& sc );
    TensorVariable   operator+(  const TensorVariable& t ) const;
    TensorVariable   operator-(  const TensorVariable& t ) const;
    TensorVariable   operator*(  const TensorVariable& t ) const;

    VectorVariable<1U>  operator*( const VectorVariable<1U>& v ) const;

    Point<1U>  operator*( const Point<1U>& v ) const;

    TensorVariable   operator/(  const TensorVariable& t ) const;
    TensorVariable&  operator+=( const TensorVariable& t );
    TensorVariable&  operator-=( const TensorVariable& t );
    TensorVariable&  operator/=( const TensorVariable& t );
    TensorVariable&  operator*=( const TensorVariable& t );
    
    TensorVariable&  operator=( double64 val );
    TensorVariable&  operator=( const ScalarVariable& s );
    TensorVariable&  operator=( const VectorVariable<1U>& v );
    TensorVariable&  operator=( const TensorVariable& t );
    TensorVariable&  operator=( TensorVariable&& ) = default;
  
    bool             operator==( const TensorVariable& t ) const; 
    bool             operator!=( const TensorVariable& t ) const; 
    bool             operator<( const TensorVariable& t ) const; 
    
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;

    /// no resizing but new values
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    void             Identity();
    double64         MinElement() const;
    double64         MaxElement() const;
    double64         Determinant() const;
    double64         Trace() const;
    TensorVariable   Adjoint()     const;
    TensorVariable   Inverse()     const;
    TensorVariable   Transposed()  const;
    bool 		  	     EigenValues( VectorVariable<1U>& vecEigenvalues ) const;
    bool             Eigen( VectorVariable<1U>& vvEigenvalues, TensorVariable<1U>& tvEigenvectors, bool bNormalize ) const;
    bool             EigenNonSymmetric( VectorVariable<1U>& eigenVals, TensorVariable<1U>& eigenVecs ) const;
    void             AssignToRow( size_t, VectorVariable<1U>& vc );
    void             AssignToColumn( size_t, VectorVariable<1U>& vc );

	  VectorVariable<1U> Row( size_t ) const;
 	  VectorVariable<1U> Column( size_t ) const;

    void             In();
    void             Out() const;
    bool             In( FILE* fp );
    bool             Out( FILE* fp ) const;

  private:
    VARIABLE_FLAG  flag;
    double64       data;
};


VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts );

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts
 
} // csmp
 
#endif

