#ifndef CSMP_TENSOR_VARIABLE2_H
#define CSMP_TENSOR_VARIABLE2_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "VectorVariable2.h"

namespace csmp {

/// @brief 2D full specialization of tensor-type variable template (other tensors are declared in files #1, #)
template<>
class TensorVariable<2U> {
  public:
    static constexpr VARIABLE_TYPE VariableType = TENSOR;

    TensorVariable() noexcept;
    TensorVariable( VARIABLE_FLAG f, double val ) noexcept;
    TensorVariable( VARIABLE_FLAG f, double v11, double v12,
                    double v21, double v22 ) noexcept;
                                  
    TensorVariable( const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
                    const double&  v11, const double&  v12,
                    const double&  v21, const double&  v22 ) noexcept;
    
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

    TensorVariable&  operator+=( const ScalarVariable& ) noexcept;
    TensorVariable&  operator-=( const ScalarVariable& ) noexcept;
    TensorVariable&  operator*=( const ScalarVariable& ) noexcept;
    TensorVariable&  operator/=( const ScalarVariable& ) noexcept;
    TensorVariable   operator+(  const TensorVariable& ) const noexcept;
    TensorVariable   operator-(  const TensorVariable& ) const noexcept;
    TensorVariable   operator*(  const TensorVariable& ) const noexcept;

    VectorVariable<2U>  operator*( const VectorVariable<2U>& ) const noexcept;

    Point<2U>  operator*( const Point<2U>& ) const noexcept;

    TensorVariable   operator/(  const TensorVariable& ) const noexcept;
    TensorVariable&  operator+=( const TensorVariable& ) noexcept;
    TensorVariable&  operator-=( const TensorVariable& ) noexcept;
    TensorVariable&  operator/=( const TensorVariable& ) noexcept;
    TensorVariable&  operator*=( const TensorVariable& ) noexcept;
    
    TensorVariable&  operator=( double val ) noexcept;
    TensorVariable&  operator=( const ScalarVariable& ) noexcept;
    TensorVariable&  operator=( const VectorVariable<2U>& ) noexcept;
  
    bool             operator==( const TensorVariable& ) const noexcept;
    bool             operator!=( const TensorVariable& ) const noexcept;
    bool             operator<( const TensorVariable& ) const noexcept; 
    
    bool             IsWithinRange( double vmin, double vmax ) const noexcept;
    bool             Has_NaN_Values() const noexcept;
    VARIABLE_FLAG&   Flag( uint32_t i=0 ) noexcept;
    VARIABLE_FLAG    Flag( uint32_t i=0 ) const noexcept;
    static constexpr uint32_t Size() noexcept { return 4u; };
    void             Identity() noexcept;
    double           MinElement() const noexcept;
    double           MaxElement() const noexcept;
    double           Determinant() const noexcept;
    double           Trace() const noexcept;
    TensorVariable   Adjoint()     const noexcept;
    TensorVariable   Inverse()     const noexcept;
    TensorVariable   Transposed()  const noexcept;
    
    bool 				     Eigen( VectorVariable<2U>& vvEigenvalues, 
                            TensorVariable<2U>& tvEigenvectors, 
                            bool bNormalize ) const;
                            
    bool 				     EigenValues( VectorVariable<2U>& vvEigenvalues) const;
    bool             EigenValues( std::vector<double>& vecEigenvalues ) const;
    bool             EigenNonSymmetric( VectorVariable<2U>& eigenVals, TensorVariable<2U>& eigenVecs ) const;
 	  void				     DiagonalValues( double f_00, double f_11) noexcept;
 	  void				     DiagonalValues( const std::vector<double>& vecDiags ) noexcept;
 	  void				     DiagonalValues( const VectorVariable<2U>& vecDiags ) noexcept;
 	  void             AssignToRow( uint32_t i, VectorVariable<2U>& vc ) noexcept;
 	  void             AssignToColumn( uint32_t j, VectorVariable<2U>& vc ) noexcept;
	  VectorVariable<2U> Row( uint32_t iRow ) const noexcept;
 	  VectorVariable<2U> Column( uint32_t iCol ) const noexcept;
 		
    void             	In();
    void             	Out() const noexcept;
    bool              In( std::fstream& fp );
    bool              Out( std::fstream& fp ) const;

  private:
    std::array<VARIABLE_FLAG,2U>          flag;
    std::array<std::array<double,2U>,2U>  data;
};


VectorVariable<2U>  operator*( const VectorVariable<2U>& vc, const TensorVariable<2U>& ts ) noexcept;
Point<2U>           operator*( const Point<2U>& vc, const TensorVariable<2U>& ts ) noexcept;

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts

} // csmp 
 
#endif

