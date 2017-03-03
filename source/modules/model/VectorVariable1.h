#ifndef CSMP_VECTOR_VARIABLE1_H
#define CSMP_VECTOR_VARIABLE1_H

#include "Point.h"
#include "ScalarVariable.h"

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
    const double64&  operator()( size_t i ) const;
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
  
    // Normal Methods
    VARIABLE_FLAG&   Flag( size_t i=0 );
    VARIABLE_FLAG    Flag( size_t i=0 ) const;
    size_t           Size() const;
  
    void             Resize( size_t newSize, double64 newValue = std::numeric_limits<double64>::quiet_NaN() );
    double64         Length() const;
    Point<1U>        P() const;
    bool             IsWithinRange( double64 vmin, double64 vmax ) const;
    VectorVariable   Flip();
    void             Invert();

    double64         DotProduct( const csmp::Point<1U>& p ) const;
    double64         DotProduct( const VectorVariable& v ) const;
    VectorVariable   CrossProduct( const csmp::Point<1U>& p ) const;
    VectorVariable   CrossProduct( const VectorVariable& v ) const;
    VectorVariable   ProjectOnto( const std::vector<double64>& v ) const;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const;

    void             In();
    void             Out() const { Out(std::cout); }
    void             Out(std::ostream& os) const;
    bool             In( FILE* fp );
    bool             Out( FILE* fp ) const;
    
    friend class TensorVariable<1U>;

  private:
    VARIABLE_FLAG flag;
    double64      data;
};

// copyright (c) 2001 by S.K. Matthai, S. Geiger & Stephen G. Roberts




} // end namespace csmp

#endif












