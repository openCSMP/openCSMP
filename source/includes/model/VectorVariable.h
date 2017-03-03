#ifndef CSMP_VECTOR_VARIABLE_H
#define CSMP_VECTOR_VARIABLE_H

#include "Point.h"
#include "ScalarVariable.h"
#include "VectorVariable1.h"
#include "VectorVariable2.h"

namespace csmp {


/**

@brief 3D full specialization of vector-type variable template (2D vectors are declared in files #1, #2)

@author S.K. Matthai
@author Stephen G. Roberts
@author S. Geiger
@date 2001

3D full specialization of vector-type variable template (2D vectors are in files #1, #2).

Implements basic vector variable type that consists of dim=spatial
dimensions entries (options are 2 x 1 and 3 x 1 for 2D and 3D models,
respectively) plus as many variable flags as vector
elements. The flags are of type enum VARIABLE_FLAG.  
VectorVariable is defined as a C++ template for the 2D and 3D
cases (In 1D you should use ScalarVariable to represent vector type
properties).  

NOTE: The definitions of the 2D constructors are in the
file VectorVariable2.h.  

 
@section motivation Motivation

A plain STL vector cannot be used to implement this datatype since there
is a need to flag each vector entry so that it can receive a special
treatment in a finite element computation, i.e. as a fixed (Dirichlet)
boundary value or the like. Also there is a need to make vector operations
as efficient as possible. This has been done by unrolling all loops.
 
 
@section applicability Applicability

Any kind of computations that involve Scalar-, Vector-, or TensorVariable
objects in CSMP. For standard element-by-element calculations, operators
are provided. This includes special operators like '^' which will 
calculate the power of each vector element for the user-specified 
value. Slightly odd is the use of the operators & and % for the scalar
and the cross products of vectors, respectively (there are no better
options, try maybe templatized enums, suggestions are welcome). 

All vector operations are already unrolled to give best performance.
Temporaries are eliminated as far as is possible and all *=, /=, +=,
-=, dot and cross product, and Length() operators / member functions
have no overhead due to assigments of flags. Hence there is nothing 
to be gained from additional handcoding of such operations.
 
 
@section examples Application Examples

Operations of the following types are possible: 
 
@code
VectorVariable<3U>  v1(PLAIN,PLAIN,PLAIN,1.,1.,1.);
cout <<"\main: Length of "<< v1 <<" is "<< v1.Length() << endl;

VectorVariable<3U>  v2;
v2 = sqrt(3.);
VectorVariable<3U>  v3 = v2 & v1; // dot product

v3.Out();
@endcode

*/
template<>
class VectorVariable<3U> {
  public:
<<<<<<< HEAD
    VectorVariable();                             ///< default constructor
    ~VectorVariable();                            ///< destructor
    VectorVariable( const VectorVariable& );      ///< copy constructor
    VectorVariable( VectorVariable&& ) = default; ///< move constructor
=======
    VectorVariable();                        ///< default constructor
    ~VectorVariable();                       ///< destructor
    VectorVariable( const VectorVariable& ); ///< copy constructor
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
    /// sets all elements to fl, val
    VectorVariable( VARIABLE_FLAG, double64 );
  
    /// constructor for 3D version
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3, 
                    double64 val1, double64 val2, double64 val3 );
     
    /// initialize with an STL vector
    explicit VectorVariable( const std::vector<double64>& );
  
    /// initialize with a Point
    explicit VectorVariable( const csmp::Point<3U>& );

    // access of vector elements
    double64&        operator()( size_t );
    const double64&  operator()( size_t ) const;
    double64         operator[]( size_t ) const;
    void             Component( size_t, double64 );
    double64         Component( size_t ) const;
    
    // assigments
    VectorVariable&  operator=( double64 );
    VectorVariable&  operator=( const Point<3U>& );
    VectorVariable&  operator=( const ScalarVariable& );
    VectorVariable&  operator=( const VectorVariable& );
    VectorVariable&  operator=( VectorVariable&& ) = default;

    // standard operators
    VectorVariable   operator+( double64 ) const;
    VectorVariable   operator-( double64 ) const;
    VectorVariable   operator*( double64 ) const;
    VectorVariable   operator/( double64 ) const;
  
    /// squares all elements of the vector
    VectorVariable   operator^( double64 ) const;

    // element-by-element operations
    VectorVariable   operator+(  const VectorVariable& ) const;
    VectorVariable   operator-(  const VectorVariable& ) const;
    VectorVariable   operator*(  const VectorVariable& ) const;
    VectorVariable   operator/(  const VectorVariable& ) const;

    /// vec1 += vec2  enables shorthand for  vec1 = vec1 + vec2
    VectorVariable&  operator+=( double64 );
    VectorVariable&  operator-=( double64 );
    VectorVariable&  operator*=( double64 );
    VectorVariable&  operator/=( double64 );
    
    VectorVariable&  operator+=( const ScalarVariable& );
    VectorVariable&  operator-=( const ScalarVariable& );
    VectorVariable&  operator*=( const ScalarVariable& );
    VectorVariable&  operator/=( const ScalarVariable& );
    
    VectorVariable&  operator+=( const VectorVariable& );
    VectorVariable&  operator-=( const VectorVariable& );
    VectorVariable&  operator*=( const VectorVariable& );
    VectorVariable&  operator/=( const VectorVariable& );
    
    /// comparison of flags and values
    bool             operator==( const VectorVariable& ) const;
    bool             operator!=( const VectorVariable& ) const;

    /// comparison by length (to allow ordering in containers)
    bool             operator<(  const VectorVariable& ) const;
  
    /// dot product = scalar product
    double64         operator&(  const VectorVariable& ) const;
  
    /// cross product (vector perpendicular to input vectors
    VectorVariable   operator%(  const VectorVariable& ) const;
    
    /// returns spatial i-th dimension
    double64        Length() const;
    double64        AngleTo( const VectorVariable& v ) const;
  
    /// returns csmp::Point initialised with vector values; @note name avoids GNU clash
    Point<3U>       P() const;
  
    /// checks vector length against the value range supplied as arguments
    bool            IsWithinRange( double64 vmin, double64 vmax ) const;
    VARIABLE_FLAG   Flag( size_t i=0 ) const;
    VARIABLE_FLAG&  Flag( size_t i=0 );
  
    /// for the PropertyStorage
    size_t          Size() const;
  
    /// sets all values to newValue; @todo SKM (1) deprecate
    void            Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
  
    // projections
    double64        DotProduct( const csmp::Point<3U>& ) const;
    double64        DotProduct( const VectorVariable& ) const;
    VectorVariable  CrossProduct( const csmp::Point<3U>& ) const;
    VectorVariable  CrossProduct( const VectorVariable& ) const;
    VectorVariable  ProjectOnto( const std::vector<double64>& ) const;
    VectorVariable  ProjectOnto( const VectorVariable& ) const;

    // assignment and modification
  
    /// Multiplies by negative unity vector
    VectorVariable  Flip();
  
    /// reverts the sequence of entries
    void            Invert();

    /// normalizes the vector variable by its length given the unit normal @todo replace with more generic function
    void            EuclideanNormalize();
    
    // IO
    void            In();
    void            Out() const { Out(std::cout); }
    void            Out(std::ostream& is) const;
    bool            In( FILE* );
    bool            Out( FILE* ) const;

    friend class TensorVariable<3U>;

  private:
    std::array<VARIABLE_FLAG,3U>  flag;  ///< flags that specify the treatment of the values in computations
    std::array<double64,3U>       data;  ///< values of the vector components
};


VectorVariable<1U> makeVector( VARIABLE_FLAG, double64 );
VectorVariable<2U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG, double64, double64 );
VectorVariable<3U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG, double64, double64, double64 );
VectorVariable<3U> makeVector( const std::array<VARIABLE_FLAG,3U>&, const std::array<double64,3U>& );
VectorVariable<3U> makeVector( const std::vector<VARIABLE_FLAG>&, const std::vector<double64>& );


template<size_t dim>
std::ostream&  operator<<( std::ostream&, const VectorVariable<dim>& );

// 1D

Point<1U> operator+( const Point<1U>&, const VectorVariable<1U>& );


Point<1U> operator-( const Point<1U>&, const VectorVariable<1U>& );


Point<1U> operator*( const Point<1U>&, const VectorVariable<1U>& );


Point<1U> operator/( const Point<1U>&, const VectorVariable<1U>& );

// 2D

Point<2U> operator+( const Point<2U>&, const VectorVariable<2U>& );


Point<2U> operator-( const Point<2U>&, const VectorVariable<2U>& );


Point<2U> operator*( const Point<2U>&, const VectorVariable<2U>& );


Point<2U> operator/( const Point<2U>&, const VectorVariable<2U>& );

// 3D

Point<3U> operator+( const Point<3U>&, const VectorVariable<3U>& );


Point<3U> operator-( const Point<3U>&, const VectorVariable<3U>& );


Point<3U> operator*( const Point<3U>&, const VectorVariable<3U>& );


Point<3U> operator/( const Point<3U>&, const VectorVariable<3U>& );


<<<<<<< HEAD
=======
// 1D specializations
inline const VectorVariable<1U>& makeVector( VARIABLE_FLAG fx, double64 vx )
 {
    return std::move(VectorVariable<1U>(fx,vx));
 }
 
// 2D specializations
inline const VectorVariable<2U>& makeVector( VARIABLE_FLAG fx, VARIABLE_FLAG fy, double64 vx, double64 vy )
 {
    return std::move(VectorVariable<2U>(fx,fy,vx,vy));
 }


// 3D specializations

inline const VectorVariable<3U>& makeVector( VARIABLE_FLAG fx, VARIABLE_FLAG fy, VARIABLE_FLAG fz, double64 vx, double64 vy, double64 vz )
 {
    return std::move(VectorVariable<3U>(fx,fy,fz,vx,vy,vz));
 }

inline const VectorVariable<3U>& makeVector( const std::array<VARIABLE_FLAG,3U>& flags, const std::array<double64,3U>& vals )
 {
    return std::move(VectorVariable<3U>(flags[0],flags[1],flags[2],vals[0],vals[1],vals[2]));
 }

inline const VectorVariable<3U>& makeVector( const std::vector<VARIABLE_FLAG>& flags, const std::vector<double64>& vals )
 {
    assert( flags.size() == 3U );
    assert( vals.size() == 3U );
    return std::move(VectorVariable<3U>(flags[0],flags[1],flags[2],vals[0],vals[1],vals[2]));
 }



inline VectorVariable<3U>::VectorVariable( const VectorVariable<3U>& v )
 : flag(v.flag),
   data(v.data)
 {
 }



inline double64& VectorVariable<3U>::operator()( size_t i )       
  { 
#ifdef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nVectorVariable<3U>::operator(): vector access violation, i="<< i << std::endl;
         return data[0];
      }
#endif
     return data[i]; 
  }



inline double64  VectorVariable<3U>::operator[]( size_t i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nVectorVariable<3U>::operator[]: vector access violation, i="<< i << std::endl;
         return data[0];
      }
#endif
    return data[i]; 
 }



inline void  VectorVariable<3U>::Component( size_t i, double64 val ) 
 { 
    assert( i < 3U );
    data[i] = val; 
 }




inline double64  VectorVariable<3U>::Component( size_t i ) const 
 { 
    assert( i < 3U );
    return data[i]; 
 }



    
    // Normal Methods

inline double64  VectorVariable<3U>::Value( const size_t& i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) {
         std::cerr <<"\nVectorVariable<3U>::Value(): vector access violation, i="<< i << std::endl;
         return data[0];
      }
#endif
    return data[i]; 
 }


inline VARIABLE_FLAG&  VectorVariable<3U>::Flag( const size_t& i )       
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nVectorVariable<3U>::Flag(): access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }


inline VARIABLE_FLAG   VectorVariable<3U>::Flag( const size_t& i ) const 
 { 
#ifndef NDEBUG 
    if ( i >= 3U ) { 
         std::cerr <<"\nVectorVariable<3U>::Flag(): access violation, i="<< i << std::endl;
         return flag[0];
      }
#endif
    return flag[i]; 
 }

inline size_t VectorVariable<3U>::Size() const
  {
    return 3U;
  }

inline void VectorVariable<3U>::Resize( size_t, double64 newValue )
  {
    data[0] = newValue;
    data[1] = newValue;
    data[2] = newValue;
  }





inline VectorVariable<3U>::~VectorVariable() 
  {
  }



inline size_t VectorVariable<3U>::Components() const { return 3U; }



inline VectorVariable<3U>::VectorVariable( const std::vector<double64>& v )
 : flag{{ANY,ANY,ANY}},
   data{{v[0],v[1],v[2]}}
 {
 }

inline VectorVariable<3U>::VectorVariable( const csmp::Point<3U>& p )
 : flag{{ANY,ANY,ANY}},
   data{{p[0],p[1],p[2]}}
 {
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator+=( double64 val )
 {
    data[0] += val;
    data[1] += val;
    data[2] += val;
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator-=( double64 val )
 {
    data[0] -= val;
    data[1] -= val;
    data[2] -= val;
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator*=( double64 val )
 {
    data[0] *= val;
    data[1] *= val;
    data[2] *= val;
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator/=( double64 val )
 {
    data[0] /= val;
    data[1] /= val;
    data[2] /= val;
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator+=( const ScalarVariable& sc )
 {
    data[0] += sc.Value();
    data[1] += sc.Value();
    data[2] += sc.Value();
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator-=( const ScalarVariable& sc )
 {
    data[0] -= sc.Value();
    data[1] -= sc.Value();
    data[2] -= sc.Value();
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator*=( const ScalarVariable& sc )
 {
    data[0] *= sc.Value();
    data[1] *= sc.Value();
    data[2] *= sc.Value();
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator/=( const ScalarVariable& sc )
 {
    data[0] /= sc.Value();
    data[1] /= sc.Value();
    data[2] /= sc.Value();
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator+=( const VectorVariable<3U>& v )
 {
    data[0] += v.data[0];
    data[1] += v.data[1];
    data[2] += v.data[2];
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator-=( const VectorVariable<3U>& v )
 {
    data[0] -= v.data[0];
    data[1] -= v.data[1];
    data[2] -= v.data[2];
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator*=( const VectorVariable<3U>& v )
 {
    data[0] *= v.data[0];
    data[1] *= v.data[1];
    data[2] *= v.data[2];
        
    return *this;  
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator/=( const VectorVariable<3U>& v )
 {
    data[0] /= v.data[0];
    data[1] /= v.data[1];
    data[2] /= v.data[2];
        
    return *this; 
 }




// --------------------
// ASSIGNMENT OPERATORS
// --------------------

inline VectorVariable<3U>&  VectorVariable<3U>::operator=( double64 val )
 {
    data[0] = val;
    data[1] = val;
    data[2] = val;
        
    return *this; 
 }




inline VectorVariable<3U>&  VectorVariable<3U>::operator=( const ScalarVariable& sc )
 {
    flag[0] = flag[1] = flag[2] = sc.Flag();
    data[0] = data[1] = data[2] = sc.Value();
       
    return *this; 
 }
 


inline VectorVariable<3U>&  VectorVariable<3U>::operator=( const csmp::Point<3U>& pt )
 {
    data[0] = pt[0];
    data[1] = pt[1];
    data[2] = pt[2];
       
    return *this; 
 }

// using the comparitor of the standard array
inline bool VectorVariable<3U>::operator==( const VectorVariable<3U>& v ) const
 {
    return( flag == v.flag && data == v.data );
 }
 

inline bool VectorVariable<3U>::operator!=( const VectorVariable<3U>& v ) const
 {
    return( flag != v.flag || data != v.data );
 }


/// compare the length of two vectors
inline bool VectorVariable<3U>::operator<( const VectorVariable<3U>& v ) const
 {
    return (this->Length() < v.Length());
 } 



/// dot product

inline double64 VectorVariable<3U>::operator&( const VectorVariable<3U>& v ) const
 {
    return data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]; 
 } 



// -------
// METHODS
// -------

/// normalize L2

inline void VectorVariable<3U>::EuclideanNormalize() 
 {
    const double64 fNorm(std::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2]));
    
    if(fNorm == 0.) return; //added AP
    
    data[0] /= fNorm;
    data[1] /= fNorm;
    data[2] /= fNorm;
 } 


inline double64 VectorVariable<3U>::DotProduct( const csmp::Point<3U>& p ) const
 {
    return data[0]*p[0] + data[1]*p[1] + data[2]*p[2]; 
 }

inline double64 VectorVariable<3U>::DotProduct( const VectorVariable& v ) const
 {
    return data[0]*v[0] + data[1]*v[1] + data[2]*v[2];
 }

inline VectorVariable<3U> VectorVariable<3U>::CrossProduct( const csmp::Point<3U>& p ) const
 {
     return VectorVariable<3U>(    flag[0], flag[1], flag[2],
                                   data[1]*p[2] - data[2]*p[1],
                                   data[2]*p[0] - data[0]*p[2],
                                   data[0]*p[1] - data[1]*p[0] );
 }

inline VectorVariable<3U> VectorVariable<3U>::CrossProduct( const VectorVariable& v ) const
 {
     return VectorVariable<3U>(    flag[0], flag[1], flag[2],
                                   data[1]*v[2] - data[2]*v[1],
                                   data[2]*v[0] - data[0]*v[2],
                                   data[0]*v[1] - data[1]*v[0] );
 }

inline double64  VectorVariable<3U>::Average() const 
 {
   return (data[0]+data[1]+data[2]) / static_cast<double64>(3.);
 }


inline double64  VectorVariable<3U>::Length() const 
 {
    return std::sqrt( data[0]*data[0] + data[1]*data[1] + data[2]*data[2] );
 }


inline void  VectorVariable<3U>::Zero()
 {
    data[0] = data[1] = data[2] = static_cast<double64>(0.);
 }



inline Point<3U>  VectorVariable<3U>::P() const
 {
    return csmp::Point<3U>(data[0],data[1],data[2]);
 }



inline bool  VectorVariable<3U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data[0] < vmin || data[0] > vmax ) return false;
    if ( data[1] < vmin || data[1] > vmax ) return false;
    if ( data[2] < vmin || data[2] > vmax ) return false;
     
    return true;
 }




inline void  VectorVariable<3U>::Fabs() 
 { 
    data[0] = std::fabs(data[0]); 
    data[1] = std::fabs(data[1]); 
    data[2] = std::fabs(data[2]); 
 }

inline  bool VectorVariable<3U>::Out( FILE* fp ) const
  {
  fwrite( (void*)this, sizeof(VectorVariable<3U>), 1, fp );
  return true;
  }

inline  bool VectorVariable<3U>::In( FILE* fp )
  {
  fread( (void*)this, sizeof(VectorVariable<3U>), 1, fp );
  return true;
  }

/// operators with Points
// 1D

inline Point<1U> operator+( const Point<1U>& p, const VectorVariable<1U>& vc )
 {
    return Point<1U>(p[0] + vc[0]);
 }


inline Point<1U> operator-( const Point<1U>& p, const VectorVariable<1U>& vc )
 {
    return Point<1U>(p[0] - vc[0]);
 }


inline Point<1U> operator*( const Point<1U>& p, const VectorVariable<1U>& vc )
 {
    return Point<1U>(p[0] * vc[0]);
 }


inline Point<1U> operator/( const Point<1U>& p, const VectorVariable<1U>& vc )
 {
    return Point<1U>(p[0] / vc[0]);
 }

// 2D

inline Point<2U> operator+( const Point<2U>& p, const VectorVariable<2U>& vc )
 {
    return Point<2U>(p[0] + vc[0], p[1] + vc[1]);
 }


inline Point<2U> operator-( const Point<2U>& p, const VectorVariable<2U>& vc )
 {
    return Point<2U>(p[0] - vc[0], p[1] - vc[1]);
 }


inline Point<2U> operator*( const Point<2U>& p, const VectorVariable<2U>& vc )
 {
    return Point<2U>(p[0] * vc[0], p[1] * vc[1]);
 }


inline Point<2U> operator/( const Point<2U>& p, const VectorVariable<2U>& vc )
 {
    return Point<2U>(p[0] / vc[0], p[1] / vc[1]);
 }

// 3D

inline Point<3U> operator+( const Point<3U>& p, const VectorVariable<3U>& vc )
 {
    return Point<3U>(p[0] + vc[0], p[1] + vc[1], p[2] + vc[2]);
 }


inline Point<3U> operator-( const Point<3U>& p, const VectorVariable<3U>& vc )
 {
    return Point<3U>(p[0] - vc[0], p[1] - vc[1], p[2] - vc[2]);
 }


inline Point<3U> operator*( const Point<3U>& p, const VectorVariable<3U>& vc )
 {
    return Point<3U>(p[0] * vc[0], p[1] * vc[1], p[2] * vc[2]);
 }


inline Point<3U> operator/( const Point<3U>& p, const VectorVariable<3U>& vc )
 {
    return Point<3U>(p[0] / vc[0], p[1] / vc[1], p[2] / vc[2]);
 }



>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
} // end namespace csmp

#endif












