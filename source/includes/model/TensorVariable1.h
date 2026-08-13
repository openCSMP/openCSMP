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

    TensorVariable() noexcept;
    TensorVariable( VARIABLE_FLAG, double ) noexcept;
    
    double&        operator()( uint32_t i, uint32_t j ) noexcept;
    const double&  operator()( uint32_t i, uint32_t j ) const noexcept;
    void           Component( uint32_t, double ) noexcept;
    [[nodiscard]] double Component( uint32_t i ) const noexcept;

    TensorVariable   operator+( double val ) const noexcept;
    TensorVariable   operator-( double val ) const noexcept;
    TensorVariable   operator*( double val ) const noexcept;
    TensorVariable   operator/( double val ) const noexcept;
    
    TensorVariable&  operator+=( double val ) noexcept;
    TensorVariable&  operator-=( double val ) noexcept;
    TensorVariable&  operator*=( double val ) noexcept;
    TensorVariable&  operator/=( double val ) noexcept;

    TensorVariable&  operator+=( const ScalarVariable& sc ) noexcept;
    TensorVariable&  operator-=( const ScalarVariable& sc ) noexcept;
    TensorVariable&  operator*=( const ScalarVariable& sc ) noexcept;
    TensorVariable&  operator/=( const ScalarVariable& sc ) noexcept;
    TensorVariable   operator+(  const TensorVariable& t ) const noexcept;
    TensorVariable   operator-(  const TensorVariable& t ) const noexcept;
    TensorVariable   operator*(  const TensorVariable& t ) const noexcept;

    VectorVariable<1U>  operator*( const VectorVariable<1U>& v ) const noexcept;

    Point<1U>  operator*( const Point<1U>& v ) const noexcept;

    TensorVariable   operator/(  const TensorVariable& t ) const noexcept;
    TensorVariable&  operator+=( const TensorVariable& t ) noexcept;
    TensorVariable&  operator-=( const TensorVariable& t ) noexcept;
    TensorVariable&  operator/=( const TensorVariable& t ) noexcept;
    TensorVariable&  operator*=( const TensorVariable& t ) noexcept;
    
    TensorVariable&  operator=( double val ) noexcept;
    TensorVariable&  operator=( const ScalarVariable& ) noexcept;
    TensorVariable&  operator=( const VectorVariable<1U>& ) noexcept;
  
    bool             operator==( const TensorVariable& ) const noexcept;
    bool             operator!=( const TensorVariable& ) const noexcept;
    bool             operator<( const TensorVariable& ) const noexcept; 
    
    bool             IsWithinRange( double vmin, double vmax ) const noexcept;
    bool             Has_NaN_Values() const noexcept { return std::isnan(data); }

    VARIABLE_FLAG&   Flag( uint32_t i=0 ) noexcept;
    VARIABLE_FLAG    Flag( uint32_t i=0 ) const noexcept;

    static constexpr uint32_t Size() noexcept { return 1u; };

    /// no resizing but new values
    void             Identity() noexcept;
    double           MinElement() const noexcept;
    double           MaxElement() const noexcept;
    double           Determinant() const noexcept;
    double           Trace() const noexcept;
    TensorVariable   Adjoint()     const noexcept;
    TensorVariable   Inverse()     const noexcept;
    TensorVariable   Transposed()  const noexcept;
    bool 		  	     EigenValues( VectorVariable<1U>& vecEigenvalues ) const;
    bool             Eigen( VectorVariable<1U>& vvEigenvalues, TensorVariable<1U>& tvEigenvectors, bool normalize_Eigen_vectors ) const;
    bool             EigenWeaklyNonSymmetric( VectorVariable<1U>& eigenVals, TensorVariable<1U>& eigenVecs, double tolerance ) const;
    void             AssignToRow( uint32_t, VectorVariable<1U>& vc ) noexcept;
    void             AssignToColumn( uint32_t, VectorVariable<1U>& vc ) noexcept;

	  VectorVariable<1U> Row( uint32_t ) const noexcept;
 	  VectorVariable<1U> Column( uint32_t ) const noexcept;

    void             In();
    void             Out() const noexcept;
    bool             In( std::fstream& fp );
    bool             Out( std::fstream& fp ) const;

  private:
    VARIABLE_FLAG  flag;
    double         data;
};


VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts ) noexcept;
Point<1U>  operator*( const Point<1U>& vc, const TensorVariable<1U>& ts ) noexcept;

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts
 
} // csmp
 
#endif

