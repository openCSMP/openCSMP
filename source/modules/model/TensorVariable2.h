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
    TensorVariable();
    TensorVariable( const TensorVariable& );
    TensorVariable( TensorVariable&& );
    TensorVariable( VARIABLE_FLAG f, double64 val );
    TensorVariable( VARIABLE_FLAG f, double64 v11, double64 v12,
                    double64 v21, double64 v22 );
                                  
    TensorVariable( const VARIABLE_FLAG& f11, const VARIABLE_FLAG& f22,
                    const double64&  v11, const double64&  v12,
                    const double64&  v21, const double64&  v22 );
                                                      
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

    TensorVariable&  operator+=( const ScalarVariable& );
    TensorVariable&  operator-=( const ScalarVariable& );
    TensorVariable&  operator*=( const ScalarVariable& );
    TensorVariable&  operator/=( const ScalarVariable& );
    TensorVariable   operator+(  const TensorVariable& ) const;
    TensorVariable   operator-(  const TensorVariable& ) const;
    TensorVariable   operator*(  const TensorVariable& ) const;

    VectorVariable<2U>  operator*( const VectorVariable<2U>& ) const;

    TensorVariable   operator/(  const TensorVariable& ) const;
    TensorVariable&  operator+=( const TensorVariable& );
    TensorVariable&  operator-=( const TensorVariable& );
    TensorVariable&  operator/=( const TensorVariable& );
    TensorVariable&  operator*=( const TensorVariable& );
    
    TensorVariable&  operator=( double64 val );
    TensorVariable&  operator=( const ScalarVariable& );
    TensorVariable&  operator=( const VectorVariable<2U>& );
    TensorVariable&  operator=( const TensorVariable& );
    TensorVariable&  operator=( TensorVariable&& );
  
    bool             operator==( const TensorVariable& ) const;
    bool             operator!=( const TensorVariable& ) const;
    bool             operator<( const TensorVariable& ) const; 
    
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    void             Identity();
    double64         MinElement() const;
    double64         MaxElement() const;
    double64         Determinant() const;
    double64         Trace() const;
    TensorVariable   Adjoint()     const;
    TensorVariable   Inverse()     const;
    TensorVariable   Transposed()  const;
    
    bool 				     Eigen( VectorVariable<2U>& vvEigenvalues, 
                            TensorVariable<2U>& tvEigenvectors, 
                            bool bNormalize ) const;
                            
    bool 				     EigenValues( VectorVariable<2U>& vvEigenvalues) const;
    bool             EigenValues( std::vector<double64>& vecEigenvalues ) const;
 	  void				     DiagonalValues( double64 f_00, double64 f_11);
 	  void				     DiagonalValues( const std::vector<double64>& vecDiags );
 	  void				     DiagonalValues( const VectorVariable<2U>& vecDiags );
 	  void             AssignToRow( size_t i, VectorVariable<2U>& vc );
 	  void             AssignToColumn( size_t j, VectorVariable<2U>& vc );
	  VectorVariable<2U> Row( size_t iRow ) const;
 	  VectorVariable<2U> Column( size_t iCol ) const;
 		
    void             	In();
    void             	Out() const;
    bool              In( FILE* fp );
    bool              Out( FILE* fp ) const;

  private:
    std::array<VARIABLE_FLAG,2U>            flag;
    std::array<std::array<double64,2U>,2U>  data;
};


VectorVariable<2U>  operator*( const VectorVariable<2U>& vc, const TensorVariable<2U>& ts );

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts

} // csmp 
 
#endif

