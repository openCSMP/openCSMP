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
    TensorVariable( VARIABLE_FLAG, double );
    
    double&        operator()( uint32_t i, uint32_t j );
    const double&  operator()( uint32_t i, uint32_t j ) const;
    void           Component( uint32_t, double );
    double         Component( uint32_t i ) const;

    TensorVariable   operator+( double val ) const;
    TensorVariable   operator-( double val ) const;
    TensorVariable   operator*( double val ) const;
    TensorVariable   operator/( double val ) const;
    
    TensorVariable&  operator+=( double val );
    TensorVariable&  operator-=( double val );
    TensorVariable&  operator*=( double val );
    TensorVariable&  operator/=( double val );

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
    
    TensorVariable&  operator=( double val );
    TensorVariable&  operator=( const ScalarVariable& );
    TensorVariable&  operator=( const VectorVariable<1U>& );
  
    bool             operator==( const TensorVariable& ) const;
    bool             operator!=( const TensorVariable& ) const;
    bool             operator<( const TensorVariable& ) const; 
    
    bool             IsWithinRange( double vmin, double vmax ) const;
    bool             Has_NaN_Values() const { return std::isnan(data); }

    VARIABLE_FLAG&   Flag( uint32_t i=0 );
    VARIABLE_FLAG    Flag( uint32_t i=0 ) const;
    uint32_t         Size() const;

    /// no resizing but new values
    void             Identity();
    double           MinElement() const;
    double           MaxElement() const;
    double           Determinant() const;
    double           Trace() const;
    TensorVariable   Adjoint()     const;
    TensorVariable   Inverse()     const;
    TensorVariable   Transposed()  const;
    bool 		  	     EigenValues( VectorVariable<1U>& vecEigenvalues ) const;
    bool             Eigen( VectorVariable<1U>& vvEigenvalues, TensorVariable<1U>& tvEigenvectors, bool bNormalize ) const;
    bool             EigenNonSymmetric( VectorVariable<1U>& eigenVals, TensorVariable<1U>& eigenVecs ) const;
    void             AssignToRow( uint32_t, VectorVariable<1U>& vc );
    void             AssignToColumn( uint32_t, VectorVariable<1U>& vc );

	  VectorVariable<1U> Row( uint32_t ) const;
 	  VectorVariable<1U> Column( uint32_t ) const;

    void             In();
    void             Out() const;
    bool             In( std::fstream& fp );
    bool             Out( std::fstream& fp ) const;

  private:
    VARIABLE_FLAG  flag;
    double         data;
};


VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts );
Point<1U>  operator*( const Point<1U>& vc, const TensorVariable<1U>& ts );

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts
 
} // csmp
 
#endif

