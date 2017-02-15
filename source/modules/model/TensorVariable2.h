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
    
    size_t           Components() const;
    double64         Value( size_t i, size_t j ) const;
    double64         Average() const;
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    void             Zero();
    void             Identity();
    double64         MinElement() const;
    double64         MaxElement() const;
    void             Fabs();
    void             Sqrt( bool from_absolute_value=false );
    void             Ln( bool from_absolute_value=false );
    void             Log10( bool from_absolute_value=false );
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
    void             	Out() const { Out(std::cout); }
    void             	Out(std::ostream& os) const;
    bool              In( FILE* fp );
    bool              Out( FILE* fp ) const;

  private:
    std::array<VARIABLE_FLAG,2U>            flag;
    std::array<std::array<double64,2U>,2U>  data;
};


VectorVariable<2U>  operator*( const VectorVariable<2U>& vc, const TensorVariable<2U>& ts );

// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts



inline double64& TensorVariable<2U>::operator()( size_t i, size_t j ) 
 {
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }




inline const double64& TensorVariable<2U>::operator()( size_t i, size_t j ) const
 {
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::operator(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j]; 
 }




inline size_t  TensorVariable<2U>::Components() const { return 4U; } 




inline void TensorVariable<2U>::Component( size_t i, double64 val )
 { 
    assert( i < Components() );
    // row by row
    if ( i == 0U )      data[0U][0U] = val;
    else if ( i == 1U ) data[0U][1U] = val;
    else if ( i == 2U ) data[1U][0U] = val;
    else                data[1U][1U] = val; // remaining case
 }



inline double64 TensorVariable<2U>::Component( size_t i ) const
 { 
    assert( i < Components() );
    // row by row
    if ( i == 0U ) return data[0U][0U];
    if ( i == 1U ) return data[0U][1U];
    if ( i == 2U ) return data[1U][0U];
    return data[1U][1U]; // remaining case
 }



inline double64 TensorVariable<2U>::Value( size_t i, size_t j ) const
 { 
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Value(): row access violation, i="<< i << std::endl;
         return data[0][0];
      }
    if ( j >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Value(): column access violation, j="<< j << std::endl;
         return data[0][0];
      }
#endif
    return data[i][j];
 }
 
 

inline VARIABLE_FLAG& TensorVariable<2U>::Flag( size_t i )      
 { 
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }



inline VARIABLE_FLAG  TensorVariable<2U>::Flag( size_t i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 2U ) { 
         std::cerr <<"\nTensorVariable<2U>::Flag(): diagonal access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }


inline size_t TensorVariable<2U>::Size() const
  {
    return 4U;
  }

inline void TensorVariable<2U>::Resize( size_t, double64 newValue )
  {
    data[0][0] = newValue;
    data[0][1] = newValue;
    data[1][0] = newValue;
    data[1][1] = newValue;
  }
 
                                   
 

inline TensorVariable<2U>::~TensorVariable() {}


// keep for storage of tensors in associative containers
inline bool  TensorVariable<2U>::operator==( const TensorVariable<2U>& ts ) const
 {
    return ( data == ts.data && flag == ts.flag );
 }


inline bool  TensorVariable<2U>::operator!=( const TensorVariable<2U>& t ) const
 {
     return !(*this == t);
 } 



inline bool  TensorVariable<2U>::operator<( const TensorVariable<2U>& t ) const
 {
     return (this->Determinant() < t.Determinant());
 } 


// -------
// METHODS
// -------



inline double64  TensorVariable<2U>::Average() const
 {
    return (data[0][0] + data[0][1] + data[1][0] + data[1][1]) / static_cast<double64>(4.);
 }




// re-tested: SKM 29-9-2001
// @test tested: O.K.

inline double64 TensorVariable<2U>::Determinant() const
 {
    return data[0][0]*data[1][1] - data[0][1]*data[1][0];
 }
 


inline double64 TensorVariable<2U>::Trace() const
 {
    return data[0][0]+data[1][1];
 }



inline void TensorVariable<2U>::AssignToRow( size_t iRow, VectorVariable<2U>& vc )
{
	if ( iRow == 0U )
	    flag[0U] = vc.Flag(0U);
	else
	    flag[1U] = vc.Flag(1U);
	data[iRow][0U] = vc[0U];
	data[iRow][1U] = vc[1U];
}


inline void TensorVariable<2U>::AssignToColumn( size_t iCol, VectorVariable<2U>& vc )
{
	if ( iCol == 0U )
	    flag[0U] = vc.Flag(0U);
	else
	    flag[1U] = vc.Flag(1U);
	data[0U][iCol] = vc[0U];
	data[1U][iCol] = vc[1U];
}




inline VectorVariable<2U> TensorVariable<2U>::Row( size_t iRow ) const
{
	return VectorVariable<2U>( flag[iRow], flag[iRow], 
	                           data[iRow][0U], data[iRow][1U] );
}


inline VectorVariable<2U> TensorVariable<2U>::Column( size_t iCol ) const
{
	return VectorVariable<2U>( flag[iCol], flag[iCol], 
	                           data[0U][iCol], data[1U][iCol] );
}


inline  bool TensorVariable<2>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(TensorVariable<2>), 1, fp );
  return true;
  }

inline  bool TensorVariable<2>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(TensorVariable<2>), 1, fp );
  return true;
  }


} // csmp 
 
#endif

