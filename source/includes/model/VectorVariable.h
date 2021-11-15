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
    static constexpr VARIABLE_TYPE VariableType = VECTOR;

    VectorVariable();                             ///< default constructor
    ~VectorVariable();                            ///< destructor
    VectorVariable( const VectorVariable& );      ///< copy constructor
    VectorVariable( VectorVariable&& ) = default; ///< move constructor
    /// sets all elements to fl, val
    VectorVariable( VARIABLE_FLAG, double );
  
    /// constructor for 3D version
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, VARIABLE_FLAG f3, 
                    double val1, double val2, double val3 );
     
    /// initialize with an STL vector
    explicit VectorVariable( const std::vector<double>& );
  
    /// initialize with a Point
    explicit VectorVariable( const csmp::Point<3U>& );

    // access of vector elements
    double&        operator()( size_t );
    const double&  operator()( size_t ) const;
    double         operator[]( size_t ) const;
    void             Component( size_t, double );
    double         Component( size_t ) const;
    
    // assigments
    VectorVariable&  operator=( double );
    VectorVariable&  operator=( const Point<3U>& );
    VectorVariable&  operator=( const ScalarVariable& );
    VectorVariable&  operator=( const VectorVariable& );
    VectorVariable&  operator=( VectorVariable&& ) = default;

    // standard operators
    VectorVariable   operator+( double ) const;
    VectorVariable   operator-( double ) const;
    VectorVariable   operator*( double ) const;
    VectorVariable   operator/( double ) const;
  
    /// squares all elements of the vector
    VectorVariable   operator^( double ) const;

    // element-by-element operations
    VectorVariable   operator+(  const VectorVariable& ) const;
    VectorVariable   operator-(  const VectorVariable& ) const;
    VectorVariable   operator*(  const VectorVariable& ) const;
    VectorVariable   operator/(  const VectorVariable& ) const;

    /// vec1 += vec2  enables shorthand for  vec1 = vec1 + vec2
    VectorVariable&  operator+=( double );
    VectorVariable&  operator-=( double );
    VectorVariable&  operator*=( double );
    VectorVariable&  operator/=( double );
    
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
    double         operator&(  const VectorVariable& ) const;
  
    /// cross product (vector perpendicular to input vectors
    VectorVariable   operator%(  const VectorVariable& ) const;
    
    /// returns spatial i-th dimension
    double        Length() const;
    double        AngleTo( const VectorVariable& v ) const;
  
    /// returns csmp::Point initialised with vector values; @note name avoids GNU clash
    Point<3U>       P() const;
  
    /// checks vector length against the value range supplied as arguments
    bool            IsWithinRange( double vmin, double vmax ) const;
    VARIABLE_FLAG   Flag( size_t i=0 ) const;
    VARIABLE_FLAG&  Flag( size_t i=0 );
  
    /// for the PropertyStorage
    size_t          Size() const;
  
    /// sets all values to newValue; @todo SKM (1) deprecate
    void            Resize( size_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() );
  
    // projections
    double        DotProduct( const csmp::Point<3U>& ) const;
    double        DotProduct( const VectorVariable& ) const;
    VectorVariable  CrossProduct( const csmp::Point<3U>& ) const;
    VectorVariable  CrossProduct( const VectorVariable& ) const;
    VectorVariable  ProjectOnto( const std::vector<double>& ) const;
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
    void            Out() const;
    bool            In( std::fstream& );
    bool            Out( std::fstream& ) const;

    friend class TensorVariable<3U>;

  private:
    std::array<VARIABLE_FLAG,3U>  flag;  ///< flags that specify the treatment of the values in computations
    std::array<double,3U>       data;  ///< values of the vector components
};


VectorVariable<1U> makeVector( VARIABLE_FLAG, double );
VectorVariable<2U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG, double, double );
VectorVariable<3U> makeVector( VARIABLE_FLAG, VARIABLE_FLAG, VARIABLE_FLAG, double, double, double );
VectorVariable<3U> makeVector( const std::array<VARIABLE_FLAG,3U>&, const std::array<double,3U>& );
VectorVariable<3U> makeVector( const std::vector<VARIABLE_FLAG>&, const std::vector<double>& );


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


} // end namespace csmp

#endif












