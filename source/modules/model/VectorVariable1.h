#ifndef CSMP_VECTOR_VARIABLE1_H
#define CSMP_VECTOR_VARIABLE1_H

#include "Point.h"
#include "ScalarVariable.h"
#include <cstdio>

namespace csmp {

/// @brief 1D full specialization of vector-type variable template (other vectors are declared in files #1, #2)
template<>
class VectorVariable<1U> {
  public:
    VectorVariable();
    VectorVariable( const VectorVariable& vc );
    VectorVariable( VARIABLE_FLAG f, double64 val );
    explicit VectorVariable( const std::vector<double64>& v );
    explicit VectorVariable( const csmp::Point<1U>& p );
    ~VectorVariable();
    
    // access of vector elements
    double64&        operator()( size_t i );
    double64         operator[]( size_t i ) const;
    void             Component( size_t, double64 );
    double64         Component( size_t i ) const;

    VectorVariable   operator+( double64 val ) const; 
    VectorVariable   operator-( double64 val ) const;
    VectorVariable   operator*( double64 val ) const; 
    VectorVariable   operator/( double64 val ) const;
    VectorVariable   operator^( double64 val ) const;

    VectorVariable&  operator+=( double64 val );
    VectorVariable&  operator-=( double64 val );
    VectorVariable&  operator*=( double64 val );
    VectorVariable&  operator/=( double64 val );
    
    VectorVariable&  operator+=( const ScalarVariable& sc );
    VectorVariable&  operator-=( const ScalarVariable& sc );
    VectorVariable&  operator*=( const ScalarVariable& sc );
    VectorVariable&  operator/=( const ScalarVariable& sc );
    
    VectorVariable   operator+(  const VectorVariable& v ) const;
    VectorVariable   operator-(  const VectorVariable& v ) const;
    VectorVariable   operator*(  const VectorVariable& v ) const; 
    VectorVariable   operator/(  const VectorVariable& v ) const;

    VectorVariable&  operator+=( const VectorVariable& v );
    VectorVariable&  operator-=( const VectorVariable& v );
    VectorVariable&  operator*=( const VectorVariable& v );
    VectorVariable&  operator/=( const VectorVariable& v );
    
    VectorVariable&  operator=( double64 val );
    VectorVariable&  operator=( const csmp::Point<1U>& p );
    VectorVariable&  operator=(  const ScalarVariable& s );
    VectorVariable&  operator=(  const VectorVariable& v );

    // extra operators
    bool             operator==( const VectorVariable& v ) const; 
    bool             operator!=( const VectorVariable& v ) const; 
    // compare length
    bool             operator<( const VectorVariable& v ) const; 
    double64         operator&( const VectorVariable& v ) const; // dot product = scalar product
    VectorVariable   operator%( const VectorVariable& v );       // cross product
    
    // Normal Methods
    size_t           Components() const;
    double64         Value( const size_t& i=0 ) const;
    VARIABLE_FLAG&   Flag( const size_t& i=0 );
    VARIABLE_FLAG    Flag( const size_t& i=0 ) const;
    size_t           Size() const;
  
    /// obviously no change to size, but only new value
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    double64         Length() const;
    double64         Average() const; // is needed, sigh !
    Point<1U>        P() const;
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    void             Zero();
    VectorVariable   Flip();
    void             Invert();
    void             Fabs();
    void             Sqrt();
    void             Ln();
    void             Log10();

    double64         DotProduct( const csmp::Point<1U>& p ) const;
    double64         DotProduct( const VectorVariable& v ) const;
    VectorVariable   CrossProduct( const csmp::Point<1U>& p ) const;
    VectorVariable   CrossProduct( const VectorVariable& v ) const;
    VectorVariable   ProjectOnto( const std::vector<double64>& v ) const;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const;

    void             In();
    void             Out() const;
    bool             In( FILE* fp );
    bool             Out( FILE* fp ) const;
    
    friend class TensorVariable<1U>;

  private:
    VARIABLE_FLAG flag;
    double64      data;
};

// copyright (c) 2001 by S.K. Matthai, S. Geiger & Stephen G. Roberts



inline VectorVariable<1U>::VectorVariable()
  : flag(ANY), data(std::numeric_limits<double64>::quiet_NaN())
  { 
  }
 


inline VectorVariable<1U>::VectorVariable( const VectorVariable& vc )
 : flag(vc.flag), data(vc.data)
 {
 }




inline VectorVariable<1U>::~VectorVariable() 
  {
  }




inline VectorVariable<1U>&  VectorVariable<1U>::operator=( const VectorVariable<1U>& v )
 {
    if ( &v != this ) {
        flag = v.flag;
        data = v.data;
       } 
    return *this; 
 }



inline double64& VectorVariable<1U>::operator()( size_t )       
  { 
     return data; 
  }



inline double64  VectorVariable<1U>::operator[]( size_t ) const 
 { 
    return data; 
 }



inline void  VectorVariable<1U>::Component( size_t, double64 val )
 { 
    data = val; 
 }



inline double64  VectorVariable<1U>::Component( size_t ) const 
 { 
    return data; 
 }

inline size_t VectorVariable<1U>::Size() const
  {
    return 1U;
  }

inline void VectorVariable<1U>::Resize( size_t, double64 newValue )
  {
    data = newValue;
  }


inline VectorVariable<1U>::VectorVariable( VARIABLE_FLAG f, double64 val )
 : flag(f), data(val)
 {
 }


 

inline VectorVariable<1U>::VectorVariable( const std::vector<double64>& v )
 : data(v[0]),flag(ANY)
 {
 }

inline VectorVariable<1U>::VectorVariable( const csmp::Point<1U>& p )
 : data(p[0]),flag(ANY)
 {
 }



inline VectorVariable<1U>  VectorVariable<1U>::operator+( const VectorVariable<1U>& v ) const
 {
    return VectorVariable( flag, data + v.data );
 }




inline VectorVariable<1U>  VectorVariable<1U>::operator-( const VectorVariable<1U>& v ) const
 {
    return VectorVariable( flag, data - v.data );
 }



inline VectorVariable<1U>  VectorVariable<1U>::operator*( const VectorVariable<1U>& v ) const 
 {
    return VectorVariable( flag, data * v.data );
 } 



inline VectorVariable<1U>  VectorVariable<1U>::operator/( const VectorVariable<1U>& v ) const
 {
    return VectorVariable( flag, data / v.data );
 }



inline VectorVariable<1U>  VectorVariable<1U>::operator+( double64 val ) const
 {
    return VectorVariable( flag, data + val );
 }



inline VectorVariable<1U>  VectorVariable<1U>::operator-( double64 val ) const
 {
    return VectorVariable( flag, data - val );
 }



inline VectorVariable<1U>  VectorVariable<1U>::operator*( double64 val ) const 
 {
    return VectorVariable( flag, data * val );
 } 



inline VectorVariable<1U>  VectorVariable<1U>::operator/( double64 val ) const
 {
    return VectorVariable( flag, data / val );
 }




inline VectorVariable<1U>  VectorVariable<1U>::operator^( double64 val ) const
 {
    return VectorVariable( flag, std::pow( data,val ) );
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator+=( double64 val )
 {
    data += val;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator-=( double64 val )
 {
    data -= val;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator*=( double64 val )
 {
    data *= val;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator/=( double64 val )
 {
    data /= val;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator+=( const ScalarVariable& sc )
 {
    data += sc.Value();
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator-=( const ScalarVariable& sc )
 {
    data -= sc.Value();
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator*=( const ScalarVariable& sc )
 {
    data *= sc.Value();
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator/=( const ScalarVariable& sc )
 {
    data /= sc.Value();
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator+=( const VectorVariable<1U>& v )
 {
    data += v.data;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator-=( const VectorVariable<1U>& v )
 {
    data -= v.data;
        
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator*=( const VectorVariable<1U>& v )
 {
    data *= v.data;
        
    return *this;  
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator/=( const VectorVariable<1U>& v )
 {
    data /= v.data;
        
    return *this; 
 }




// --------------------
// ASSIGNMENT OPERATORS
// --------------------
inline VectorVariable<1U>&  VectorVariable<1U>::operator=( double64 val )
 {
    data = val;
    return *this; 
 }




inline VectorVariable<1U>&  VectorVariable<1U>::operator=( const csmp::Point<1U>& p )
 {
    data = p[0];
    return *this; 
 }



inline VectorVariable<1U>&  VectorVariable<1U>::operator=( const ScalarVariable& sc )
 {
    flag = sc.Flag();
    data = sc.Value();
       
    return *this; 
 }



inline bool VectorVariable<1U>::operator==( const VectorVariable<1U>& v ) const
 {
    return( (v.flag==flag && v.data==data) );
 }

 

inline bool VectorVariable<1U>::operator!=( const VectorVariable<1U>& v ) const
 {
    return( v.flag!=flag || v.data!=data );
 } 



/// dot product

inline double64 VectorVariable<1U>::operator&( const VectorVariable<1U>& v ) const
 {
    return data*v.data;  
 } 



/// compare magnitude

inline bool VectorVariable<1U>::operator<( const VectorVariable<1U>& v ) const
 {
    return( (data<v.data) );
 } 


/// cross product ' % ' of two vectors 
inline VectorVariable<1U>  VectorVariable<1U>::operator%( const VectorVariable<1U>& )
 {
    // assuming that the third 1Uension has a zero coordinate
    return VectorVariable( flag, static_cast<double64>(0.) );
 } 




// -------
// METHODS
// -------


inline double64 VectorVariable<1U>::DotProduct( const csmp::Point<1U>& p ) const
 {
    return data * p[0]; 
 }

inline double64 VectorVariable<1U>::DotProduct( const VectorVariable& v ) const
 {
    return data * v[0];
 }


// not defined in 1D -> degenerate result = 0.

inline VectorVariable<1U> VectorVariable<1U>::CrossProduct( const csmp::Point<1U>& ) const
 {
     return VectorVariable<1U>( flag, 0. );
 }

inline VectorVariable<1U> VectorVariable<1U>::CrossProduct( const VectorVariable& ) const
 {
     return VectorVariable<1U>( flag, 0. );
 }


inline VectorVariable<1U>   VectorVariable<1U>::ProjectOnto( const std::vector<double64>& v ) const
 {
    return VectorVariable<1U>( flag, v[0] );
 }



inline VectorVariable<1U>   VectorVariable<1U>::ProjectOnto( const VectorVariable& v ) const
 {
    return VectorVariable<1U>( flag, v.data );
 }




inline size_t  VectorVariable<1U>::Components() const { return 1U; }



inline double64  VectorVariable<1U>::Value( const size_t& ) const 
 { 
    return data; 
 }


inline double64  VectorVariable<1U>::Average() const 
 { 
    return data; 
 }


inline VARIABLE_FLAG&  VectorVariable<1U>::Flag( const size_t& )       
 { 
    return flag; 
 }


inline VARIABLE_FLAG   VectorVariable<1U>::Flag( const size_t& ) const 
 { 
    return flag; 
 }


inline double64  VectorVariable<1U>::Length() const 
 {
    return std::fabs(data);
 }


inline void  VectorVariable<1U>::Zero()
 {
    data = static_cast<double64>(0.);
 }


inline VectorVariable<1U>  VectorVariable<1U>::Flip()
 {
    return VectorVariable<1U>( flag, -data );
 }



inline Point<1U>  VectorVariable<1U>::P() const
 {
    return csmp::Point<1U>(data);
 }




inline bool  VectorVariable<1U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data < vmin || data > vmax ) return false;
     
    return true;
 }




inline void  VectorVariable<1U>::Fabs() 
 { 
    data = std::fabs(data); 
 }




inline void  VectorVariable<1U>::Sqrt() 
 { 
     data = std::sqrt(data);
 }




inline void  VectorVariable<1U>::Ln() 
 { 
     data = std::log(data);
 }



inline void  VectorVariable<1U>::Log10() 
 { 
     data = std::log10(data); 
 }

inline  bool VectorVariable<1U>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(VectorVariable<1U>), 1, fp );
  return true;
  }

inline  bool VectorVariable<1U>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(VectorVariable<1U>), 1, fp );
  return true;
  }


} // end namespace csp

#endif












