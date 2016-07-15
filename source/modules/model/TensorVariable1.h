#ifndef CSMP_TENSOR_VARIABLE1_H
#define CSMP_TENSOR_VARIABLE1_H

#include "ScalarVariable.h"
#include "VectorVariable1.h"

namespace csmp {

/// @brief 1D full specialization of tensor-type variable template (other tensors are declared in files #1, #2)
template<>
class TensorVariable<1U> {
  public:
    TensorVariable();
    TensorVariable( const TensorVariable& t );
    TensorVariable( VARIABLE_FLAG f, double64 val );
                                                      
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

    TensorVariable   operator/(  const TensorVariable& t ) const;
    TensorVariable&  operator+=( const TensorVariable& t );
    TensorVariable&  operator-=( const TensorVariable& t );
    TensorVariable&  operator/=( const TensorVariable& t );
    TensorVariable&  operator*=( const TensorVariable& t );
    
    TensorVariable&  operator=( double64 val );
    TensorVariable&  operator=( const ScalarVariable& s );
    TensorVariable&  operator=( const VectorVariable<1U>& v );
    TensorVariable&  operator=( const TensorVariable& t );
    
    bool             operator==( const TensorVariable& t ) const; 
    bool             operator!=( const TensorVariable& t ) const; 
    bool             operator<( const TensorVariable& t ) const; 
    
    size_t           Components() const;
    double64         Value( size_t i, size_t j ) const;
    double64         Average() const;
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;
    // no resizing but new values
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
    bool 		  	     EigenValues( VectorVariable<1U>& vecEigenvalues ) const;
    bool             Eigen( VectorVariable<1U>& vvEigenvalues, TensorVariable<1U>& tvEigenvectors, bool bNormalize ) const;
    
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



inline double64& TensorVariable<1U>::operator()( size_t, size_t ) 
 {
    return data; 
 }


inline const double64& TensorVariable<1U>::operator()( size_t, size_t ) const
 {
    return data; 
 }



inline double64 TensorVariable<1U>::Component( size_t ) const
 {
    return data; 
 }



inline void  TensorVariable<1U>::Component( size_t, double64 val ) 
 {
    data = val; 
 }



inline size_t  TensorVariable<1U>::Components() const { return 1U; } 



inline double64  TensorVariable<1U>::Value( size_t, size_t ) const 
 { 
    return data;
 }
 
 

inline VARIABLE_FLAG& TensorVariable<1U>::Flag( size_t )      
 { 
    return flag; 
 }



inline VARIABLE_FLAG  TensorVariable<1U>::Flag( size_t ) const 
 { 
    return flag; 
 }

inline size_t TensorVariable<1U>::Size() const
  {
    return 1U;
  }

inline void TensorVariable<1U>::Resize( size_t, double64 newValue )
  {
    data = newValue;
  }


inline TensorVariable<1U>::TensorVariable()
   : flag(ANY), data(std::numeric_limits<double64>::quiet_NaN())
  {
  }



inline TensorVariable<1U>&  TensorVariable<1U>::operator=( const TensorVariable<1U>& ts )
 {
    if ( &ts != this ) {
         flag = ts.flag;
         data = ts.data;
      }  
    return *this; 
 }




inline TensorVariable<1U>::TensorVariable( const TensorVariable<1U>& t )
 {
    *this = t;
 }




inline TensorVariable<1U>::TensorVariable( VARIABLE_FLAG f, double64 val )
 : flag(f), data(val)
 {
 } 
 
                                   
 

inline TensorVariable<1U>::~TensorVariable() {}




inline void  TensorVariable<1U>::Zero()
 {
    data = static_cast<double64>(0.0);
 }
 
 

// here the flag of the lefthand tensor-variable is sustained

inline TensorVariable<1U>  TensorVariable<1U>::operator+( const TensorVariable<1U>& t ) const
 {
    return TensorVariable<1U>( flag, t.data+data );
 }



inline TensorVariable<1U>  TensorVariable<1U>::operator-( const TensorVariable<1U>& t ) const
 {
      return TensorVariable<1U>( flag, data-t.data );
 }



inline TensorVariable<1U>  TensorVariable<1U>::operator+( double64 val ) const
 {
      return TensorVariable<1U>( flag, data+val );
 }
 
 

inline TensorVariable<1U>  TensorVariable<1U>::operator-( double64 val ) const
 {
      return TensorVariable<1U>( flag, data-val );
 }
 
 

inline TensorVariable<1U>  TensorVariable<1U>::operator*( double64 val ) const
 {
      return TensorVariable<1U>( flag, data*val );
 }
 
 

inline TensorVariable<1U>  TensorVariable<1U>::operator/( double64 val ) const
 {
      return TensorVariable<1U>( flag, data/val );
 }



// matrix vector multiplication: v = M * v

inline VectorVariable<1U>  TensorVariable<1U>::operator*( const VectorVariable<1U>& vc ) const
 {
    return VectorVariable<1U>( flag, data * vc.Value() );
 } 




inline TensorVariable<1U> TensorVariable<1U>::Adjoint() const
 {
      return TensorVariable( flag, 1. );
 }



inline TensorVariable<1U>  TensorVariable<1U>::operator*( const TensorVariable<1U>& ts ) const 
 {
   return TensorVariable<1U>( flag, data * ts.data );
 } 




inline TensorVariable<1U>&  TensorVariable<1U>::operator+=( const ScalarVariable& sc )
 {
    data += sc.Value();
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator-=( const ScalarVariable& sc )
 {
    data -= sc.Value();
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator*=( const ScalarVariable& sc )
 {
    data *= sc.Value();
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator/=( const ScalarVariable& sc )
 {
    data /= sc.Value();
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator+=( const TensorVariable<1U>& ts )
 {
    data += ts.data;
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator-=( const TensorVariable<1U>& ts )
 {
    data -= ts.data;
    return *this; 
 }


// element by element division

inline TensorVariable<1U>&  TensorVariable<1U>::operator/=( const TensorVariable<1U>& ts )
 {
    data /= ts.data;
    return *this; 
 }



// element by element division

inline TensorVariable<1U>  TensorVariable<1U>::operator/( const TensorVariable<1U>& ts ) const
 {
    return TensorVariable<1U>( flag, data / ts.data ); 
 }



// matrix multiplication

inline TensorVariable<1U>&  TensorVariable<1U>::operator*=( const TensorVariable<1U>& ts )
 {
    data *= ts.data;
    return *this; 
 }





inline TensorVariable<1U>&  TensorVariable<1U>::operator+=( double64 val )
 {
    data += val;
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator-=( double64 val )
 {
    data -= val;
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator*=( double64 val )
 {
    data *= val;
    return *this; 
 }



inline TensorVariable<1U>&  TensorVariable<1U>::operator/=( double64 val )
 {
    data /= val;
    return *this; 
 }


// --------------------
// ASSIGNMENT OPERATORS
// --------------------


inline TensorVariable<1U>&  TensorVariable<1U>::operator=( double64 val )
 {
    data = val;
    return *this; 
 }


// the flag is adopted from the scalar variable

inline TensorVariable<1U>&  TensorVariable<1U>::operator=( const ScalarVariable& sc )
 {
    flag = sc.Flag();
    data = sc.Value();
    return *this; 
 }



// writes vector into the diagonal of the zero'd tensor

inline TensorVariable<1U>&  TensorVariable<1U>::operator=( const VectorVariable<1U>& vc )
 {
    flag = vc.Flag();
    data = vc.Value();
    return *this; 
 }





inline bool  TensorVariable<1U>::operator==( const TensorVariable<1U>& ts ) const
 {
    return ( data == ts.data && flag == ts.flag );
 }

 
  

inline bool  TensorVariable<1U>::operator!=( const TensorVariable<1U>& t ) const
 {
     return !(*this == t);
 } 



inline bool  TensorVariable<1U>::operator<( const TensorVariable<1U>& t ) const
 {
     return (this->data < t.data);
 } 


// -------
// METHODS
// -------



inline double64  TensorVariable<1U>::Average() const
 {
    return data;
 }



inline void TensorVariable<1U>::Identity()
 {
    data = static_cast<double64>(1.0);
 }



inline TensorVariable<1U>  TensorVariable<1U>::Transposed() const
 {
    return TensorVariable( flag, data );
 }



inline double64 TensorVariable<1U>::Determinant() const
 {
    return data;
 }
 
inline double64 TensorVariable<1U>::Trace() const
 {
    return data;
 }

inline TensorVariable<1U> TensorVariable<1U>::Inverse() const
 {
    return TensorVariable( flag, 1. / data );
 }
 
 


inline  double64  TensorVariable<1U>::MinElement() const
 {
    return data;
 }
  
  


inline  double64  TensorVariable<1U>::MaxElement() const
 {
    return data;
 }
 
 


inline bool  TensorVariable<1U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data < vmin || data > vmax ) return false;
    return true;
 }
 
 


inline void  TensorVariable<1U>::Fabs() 
 { 
    data = std::fabs( data );
 }



inline void  TensorVariable<1U>::Sqrt( bool from_absolute_value ) 
 { 
    if ( from_absolute_value ) {
         data = std::sqrt( std::fabs(data) );
         return;
      }
    data = std::sqrt( data );
 }



inline void  TensorVariable<1U>::Ln( bool from_absolute_value ) 
 { 
    if ( from_absolute_value ) {
         data = std::log( std::fabs( data ) );
         return;
      }
    data = std::log( data );
 }



inline void  TensorVariable<1U>::Log10( bool from_absolute_value ) 
 { 
    if ( from_absolute_value ) {
         data = std::log10( std::fabs( data ) );
         return;
      }
    data = std::log10( data );
 }
 
 
// vector-matrix multiplication: v^T = (v^T * A)^T = A^T v  

inline VectorVariable<1U>  operator*( const VectorVariable<1U>& vc, const TensorVariable<1U>& ts )
 {
    return VectorVariable<1U>( vc.Flag(), ts.Value(0,0) * vc.Value() );
 }



inline void TensorVariable<1U>::AssignToRow( size_t, VectorVariable<1U>& vc )
{
	flag = vc.Flag(0U); 
	data = vc[0U];
}


inline void TensorVariable<1U>::AssignToColumn( size_t, VectorVariable<1U>& vc )
{
	flag = vc.Flag(0U);
	data = vc[0U];
}




inline VectorVariable<1U> TensorVariable<1U>::Row( size_t ) const
{
	return VectorVariable<1U>( flag, data );
}


inline VectorVariable<1U> TensorVariable<1U>::Column( size_t ) const
{
	return VectorVariable<1U>( flag, data );
}

inline  bool TensorVariable<1>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(TensorVariable<1>), 1, fp );
  return true;
  }

inline  bool TensorVariable<1>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(TensorVariable<1>), 1, fp );
  return true;
  }

 
} // csp 
 
#endif

