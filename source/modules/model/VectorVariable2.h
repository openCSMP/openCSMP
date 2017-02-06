#ifndef CSMP_VECTOR_VARIABLE2_H
#define CSMP_VECTOR_VARIABLE2_H

#include "Point.h"
#include "ScalarVariable.h"

namespace csmp {

/// @brief 2D full specialization of vector-type variable template (other vectors are declared in files #1, #2)
template<>
class VectorVariable<2U> {
  public:
    VectorVariable();
    VectorVariable( const VectorVariable& );
    VectorVariable( VectorVariable&& );
    VectorVariable( VARIABLE_FLAG f, double64 val );
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double64 val1, double64 val2 ); // 2d
    explicit VectorVariable( const std::vector<double64>& );
    explicit VectorVariable( const csmp::Point<2U>& );
    ~VectorVariable();
    
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
    
    VectorVariable&  operator+=( const ScalarVariable& );
    VectorVariable&  operator-=( const ScalarVariable& );
    VectorVariable&  operator*=( const ScalarVariable& );
    VectorVariable&  operator/=( const ScalarVariable& );
    
    VectorVariable   operator+(  const VectorVariable& ) const;
    VectorVariable   operator-(  const VectorVariable& ) const;
    VectorVariable   operator*(  const VectorVariable& ) const;
    VectorVariable   operator/(  const VectorVariable& ) const;

    VectorVariable&  operator+=( const VectorVariable& );
    VectorVariable&  operator-=( const VectorVariable& );
    VectorVariable&  operator*=( const VectorVariable& );
    VectorVariable&  operator/=( const VectorVariable& );
    
    VectorVariable&  operator=( double64 );
    VectorVariable&  operator=( const Point<2U>& );
    VectorVariable&  operator=( const ScalarVariable& );
    VectorVariable&  operator=( const VectorVariable& );
    VectorVariable&  operator=( VectorVariable&& );

    // extra operators
    bool             operator==( const VectorVariable& ) const;
    bool             operator!=( const VectorVariable& ) const;
    // compare length
    bool             operator<(  const VectorVariable& ) const;
    double64         operator&(  const VectorVariable& ) const; // dot product = scalar product
    VectorVariable   operator%(  const VectorVariable& ) const; // cross product
    
    // Normal Methods
    size_t           Components() const;
    double64         Value( const size_t& i ) const;
    VARIABLE_FLAG&   Flag( const size_t& i=0 );
    VARIABLE_FLAG    Flag( const size_t& i=0 ) const;
    size_t           Size() const;
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    double64         Average() const;
    double64         Length() const;
    double64         AngleTo( const VectorVariable& v ) const;
    Point<2U>        P() const;
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    void             Zero();
    VectorVariable   Flip();
    VectorVariable   ProjectOnto( const std::vector<double64>& v ) const;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const;
    void             Invert();
    void             Fabs();
    void             Sqrt();
    void             Ln();
    void             Log10();
    
    void			 EuclideanNormalize();
    double64         DotProduct( const csmp::Point<2U>& p ) const;
    double64         DotProduct( const VectorVariable& v ) const;
    VectorVariable   CrossProduct( const csmp::Point<2U>& p ) const;
    VectorVariable   CrossProduct( const VectorVariable& v ) const;

    void             In();
    void             Out(std::ostream& is) const;
    bool             In( FILE* fp );
    bool             Out( FILE* fp ) const;

    friend class TensorVariable<2U>;

  private:
    std::array<VARIABLE_FLAG,2U> flag;
    std::array<double64,2U>      data;
};


// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts



inline VectorVariable<2U>::VectorVariable()
  : flag{{ANY,ANY}},
    data{{std::numeric_limits<double64>::quiet_NaN(),std::numeric_limits<double64>::quiet_NaN()}}
  { 
  }

  
inline VectorVariable<2U>::VectorVariable( const VectorVariable<2U>& v )
  : flag(v.flag),
    data(v.data)
 {
 }


inline VectorVariable<2U>::VectorVariable( VectorVariable<2U>&& v )
  : flag{v.flag},
    data{v.data}
 {
 }


inline VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f, double64 val )
  : flag{{f,f}},
    data{{val,val}}
 {
 }


inline VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double64 val1, double64 val2 )
  : flag{{f1,f2}},
    data{{val1,val2}}
 {
 }


inline VectorVariable<2U>::VectorVariable( const std::vector<double64>& v )
  : flag{{ANY,ANY}},
    data{{v[0],v[1]}}
 {
 }


inline VectorVariable<2U>::VectorVariable( const Point<2U>& p )
  : flag{{ANY,ANY}},
    data{{p[0],p[1]}}
 {
 }


inline VectorVariable<2U>::~VectorVariable() 
  {
  }


inline VectorVariable<2U>&  VectorVariable<2U>::operator=( const VectorVariable<2U>& v )
 {
    if ( &v != this ) {
          flag = v.flag;
          data = v.data;
       } 
    return *this; 
 }


inline VectorVariable<2U>&  VectorVariable<2U>::operator=( VectorVariable<2U>&& v )
 {
    if ( &v != this ) {
          flag = {v.flag};
          data = {v.data};
       } 
    return *this; 
 }


inline double64& VectorVariable<2U>::operator()( size_t i )       
  { 
     if ( i==0U ) return data[0];
     return              data[1]; 
  }



inline double64  VectorVariable<2U>::operator[]( size_t i ) const 
 { 
    if ( i==0U ) return data[0];
    return              data[1]; 
 }



inline void  VectorVariable<2U>::Component( size_t i, double64 val ) 
 { 
    if ( i==0U ) data[0] = val;
    else         data[1] = val; 
 }



inline double64  VectorVariable<2U>::Component( size_t i ) const 
 { 
    if ( i==0U ) return data[0];
    return              data[1]; 
 }


inline size_t VectorVariable<2U>::Size() const
  {
    return 2U;
  }

inline void VectorVariable<2U>::Resize( size_t, double64 newValue )
  {
    data[0] = newValue;
    data[1] = newValue;
  }

inline VectorVariable<2U>  VectorVariable<2U>::operator+( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] + v.data[0], data[1] + v.data[1] ));
 }




inline VectorVariable<2U>  VectorVariable<2U>::operator-( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] - v.data[0], data[1] - v.data[1] ));
 }



inline VectorVariable<2U>  VectorVariable<2U>::operator*( const VectorVariable<2U>& v ) const 
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] * v.data[0], data[1] * v.data[1] ));
 } 



inline VectorVariable<2U>  VectorVariable<2U>::operator/( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] / v.data[0], data[1] / v.data[1] ));
 }



inline VectorVariable<2U>  VectorVariable<2U>::operator+( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] + val, data[1] + val ));
 }



inline VectorVariable<2U>  VectorVariable<2U>::operator-( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] - val, data[1] - val ));
 }



inline VectorVariable<2U>  VectorVariable<2U>::operator*( double64 val ) const 
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] * val, data[1] * val ));
 } 



inline VectorVariable<2U>  VectorVariable<2U>::operator/( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] / val, data[1] / val ));
 }




inline VectorVariable<2U>  VectorVariable<2U>::operator^( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], std::pow( data[0],val ), std::pow( data[1], val ) ));
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator+=( double64 val )
 {
    data[0] += val;
    data[1] += val;
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator-=( double64 val )
 {
    data[0] -= val;
    data[1] -= val;
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator*=( double64 val )
 {
    data[0] *= val;
    data[1] *= val;
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator/=( double64 val )
 {
    data[0] /= val;
    data[1] /= val;
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator+=( const ScalarVariable& sc )
 {
    data[0] += sc.Value();
    data[1] += sc.Value();
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator-=( const ScalarVariable& sc )
 {
    data[0] -= sc.Value();
    data[1] -= sc.Value();
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator*=( const ScalarVariable& sc )
 {
    data[0] *= sc.Value();
    data[1] *= sc.Value();
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator/=( const ScalarVariable& sc )
 {
    data[0] /= sc.Value();
    data[1] /= sc.Value();
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator+=( const VectorVariable<2U>& v )
 {
    data[0] += v.data[0];
    data[1] += v.data[1];
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator-=( const VectorVariable<2U>& v )
 {
    data[0] -= v.data[0];
    data[1] -= v.data[1];
        
    return *this; 
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator*=( const VectorVariable<2U>& v )
 {
    data[0] *= v.data[0];
    data[1] *= v.data[1];
        
    return *this;  
 }




inline VectorVariable<2U>&  VectorVariable<2U>::operator/=( const VectorVariable<2U>& v )
 {
    data[0] /= v.data[0];
    data[1] /= v.data[1];
        
    return *this; 
 }




// --------------------
// ASSIGNMENT OPERATORS
// --------------------

inline VectorVariable<2U>&  VectorVariable<2U>::operator=( double64 val )
 {
    data[0] = val;
    data[1] = val;
        
    return *this; 
 }



inline VectorVariable<2U>&  VectorVariable<2U>::operator=( const csmp::Point<2U>& p )
 {
    data[0] = p[0];
    data[1] = p[1];
        
    return *this; 
 }



inline VectorVariable<2U>&  VectorVariable<2U>::operator=( const ScalarVariable& sc )
 {
    flag[0] = flag[1] = sc.Flag();
    data[0] = data[1] = sc.Value();
       
    return *this; 
 }


inline bool VectorVariable<2U>::operator==( const VectorVariable<2U>& v ) const
 {
    return( flag == v.flag && data == v.data );
 }


inline bool VectorVariable<2U>::operator!=( const VectorVariable<2U>& v ) const
 {
    return( flag != v.flag || data != v.data );
 }


// compare the length of two vectors

inline bool VectorVariable<2U>::operator<( const VectorVariable<2U>& v ) const
 {
    return (this->Length() < v.Length());
 } 



// dot product

inline double64 VectorVariable<2U>::operator&( const VectorVariable<2U>& v ) const
 {
    return data[0]*v.data[0] + data[1]*v.data[1];  
 } 




/// cross product ' % ' of two vectors 
inline VectorVariable<2U>  VectorVariable<2U>::operator%( const VectorVariable<2U>& v ) const
 {
    // assuming that the third dimension has a zero coordinate
    return std::move(VectorVariable( flag[0], flag[1], 0., data[0]*v.data[1] - v.data[0]*data[1] ));
 } 



// -------
// METHODS
// -------

/// L2 norm
inline void VectorVariable<2U>::EuclideanNormalize() 
 {
    const double64 fNorm(std::sqrt(data[0]*data[0] + data[1]*data[1]));
    
    if(fNorm == 0.) return; //added AP
    
    data[0] /= fNorm;
    data[1] /= fNorm;
 } 



inline double64 VectorVariable<2U>::DotProduct( const csmp::Point<2U>& p ) const
 {
    return data[0] * p[0] + data[1] * p[1]; 
 }

inline double64 VectorVariable<2U>::DotProduct( const VectorVariable& v ) const
 {
    return data[0] * v[0] + data[1] * v[1];
 }


inline VectorVariable<2U> VectorVariable<2U>::CrossProduct( const csmp::Point<2U>& p ) const
 {
     return std::move(VectorVariable<2U>( flag[0], flag[1], 0., data[0]*p[1] - data[1]*p[0] ));
 }

inline VectorVariable<2U> VectorVariable<2U>::CrossProduct( const VectorVariable& v ) const
 {
     return std::move(VectorVariable<2U>( flag[0], flag[1], 0., data[0]*v[1] - data[1]*v[0] ));
 }

inline size_t  VectorVariable<2U>::Components() const { return 2U; } 



inline double64  VectorVariable<2U>::Value( const size_t& i ) const
 { 
    if ( i==0U ) return data[0];
    return data[1]; 
 }


inline double64  VectorVariable<2U>::Average() const 
 {
    return (data[0]+data[1]) / 2.0;
 }


inline double64  VectorVariable<2U>::Length() const 
 {
    return std::sqrt( data[0]*data[0] + data[1]*data[1] );
 }



inline VARIABLE_FLAG&  VectorVariable<2U>::Flag( const size_t& i )       
 { 
    if ( i==0U ) return flag[0];
    return flag[1]; 
 }


inline VARIABLE_FLAG  VectorVariable<2U>::Flag( const size_t& i ) const 
 { 
    if ( i==0U ) return flag[0];
    return flag[1]; 
 }


inline void  VectorVariable<2U>::Zero()
 {
    data[0] = data[1] = static_cast<double64>(0.);
 }




inline Point<2U>  VectorVariable<2U>::P() const
 {
    return csmp::Point<2U>(data[0],data[1]);
 }



inline bool  VectorVariable<2U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data[0] < vmin || data[0] > vmax ) return false;
    if ( data[1] < vmin || data[1] > vmax ) return false;
     
    return true;
 }



inline void  VectorVariable<2U>::Fabs() 
 { 
    data[0] = std::fabs(data[0]); 
    data[1] = std::fabs(data[1]); 
 }


inline void  VectorVariable<2U>::Sqrt() 
 { 
     data[0] = std::sqrt(data[0]);
     data[1] = std::sqrt(data[1]);
 }




inline void  VectorVariable<2U>::Ln() 
 { 
     data[0] = std::log(data[0]);
     data[1] = std::log(data[1]);
 }



inline void  VectorVariable<2U>::Log10() 
 { 
     data[0] = std::log10(data[0]);
     data[1] = std::log10(data[1]);
 }


inline  bool VectorVariable<2U>::Out( FILE* fp ) const
  {
     fwrite( (void*)this, sizeof(VectorVariable<2U>), 1, fp );
     return true;
  }

inline  bool VectorVariable<2U>::In( FILE* fp )
  {
     fread( (void*)this, sizeof(VectorVariable<2U>), 1, fp );
     return true;
  }


} // end namespace csmp

#endif












