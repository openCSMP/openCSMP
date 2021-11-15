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

    TensorVariable();
    TensorVariable( const TensorVariable& );
    TensorVariable( TensorVariable&& ) = default;
    TensorVariable( VARIABLE_FLAG f, double val );
    TensorVariable( VARIABLE_FLAG f, double v11, double v12,
                    double v21, double v22 );
                                  
    TensorVariable( const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
                    const double&  v11, const double&  v12,
                    const double&  v21, const double&  v22 );
                                                      
    ~TensorVariable();
    
    double&        operator()( size_t i, size_t j );
    const double&  operator()( size_t i, size_t j ) const;
    void             Component( size_t, double );
    double         Component( size_t i ) const;

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
    TensorVariable   operator+(  const TensorVariable& ) const;
    TensorVariable   operator-(  const TensorVariable& ) const;
    TensorVariable   operator*(  const TensorVariable& ) const;

    VectorVariable<2U>  operator*( const VectorVariable<2U>& ) const;

    Point<2U>  operator*( const Point<2U>& ) const;

    TensorVariable   operator/(  const TensorVariable& ) const;
    TensorVariable&  operator+=( const TensorVariable& );
    TensorVariable&  operator-=( const TensorVariable& );
    TensorVariable&  operator/=( const TensorVariable& );
    TensorVariable&  operator*=( const TensorVariable& );
    
    TensorVariable&  operator=( double val );
    TensorVariable&  operator=( const ScalarVariable& );
    TensorVariable&  operator=( const VectorVariable<2U>& );
    TensorVariable&  operator=( const TensorVariable& );
    TensorVariable&  operator=( TensorVariable&& ) = default;
  
    bool             operator==( const TensorVariable& ) const;
    bool             operator!=( const TensorVariable& ) const;
    bool             operator<( const TensorVariable& ) const; 
    
    bool             IsWithinRange( double vmin, double vmax ) const;
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;
    void             Resize( size_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() );
    void             Identity();
    double         MinElement() const;
    double         MaxElement() const;
    double         Determinant() const;
    double         Trace() const;
    TensorVariable   Adjoint()     const;
    TensorVariable   Inverse()     const;
    TensorVariable   Transposed()  const;
    
    bool 				     Eigen( VectorVariable<2U>& vvEigenvalues, 
                            TensorVariable<2U>& tvEigenvectors, 
                            bool bNormalize ) const;
                            
    bool 				     EigenValues( VectorVariable<2U>& vvEigenvalues) const;
    bool             EigenValues( std::vector<double>& vecEigenvalues ) const;
    bool             EigenNonSymmetric( VectorVariable<2U>& eigenVals, TensorVariable<2U>& eigenVecs ) const;
 	  void				     DiagonalValues( double f_00, double f_11);
 	  void				     DiagonalValues( const std::vector<double>& vecDiags );
 	  void				     DiagonalValues( const VectorVariable<2U>& vecDiags );
 	  void             AssignToRow( size_t i, VectorVariable<2U>& vc );
 	  void             AssignToColumn( size_t j, VectorVariable<2U>& vc );
	  VectorVariable<2U> Row( size_t iRow ) const;
 	  VectorVariable<2U> Column( size_t iCol ) const;
 		
    void             	In();
    void             	Out() const;
    bool              In( std::fstream& fp );
    bool              Out( std::fstream& fp ) const;

  private:
    std::array<VARIABLE_FLAG,2U>            flag;
    std::array<std::array<double,2U>,2U>  data;
};


VectorVariable<2U>  operator*( const VectorVariable<2U>& vc, const TensorVariable<2U>& ts );
Point<2U>  operator*( const Point<2U>& vc, const TensorVariable<2U>& ts );

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts

} // csmp 
 
#endif

