#ifndef CSMP_VECTOR_VARIABLE2_H
#define CSMP_VECTOR_VARIABLE2_H

#include "Point.h"
#include "ScalarVariable.h"

namespace csmp {

/// @brief 2D full specialization of vector-type variable template (other vectors are declared in files #1, #2)
template<>
class VectorVariable<2U> {
  public:
    static constexpr VARIABLE_TYPE VariableType = VECTOR;

    VectorVariable();
    VectorVariable( const VectorVariable& );
    VectorVariable( VectorVariable&& ) = default;
    VectorVariable( VARIABLE_FLAG f, double val );
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double val1, double val2 );
    explicit VectorVariable( const std::vector<double>& );
    explicit VectorVariable( const csmp::Point<2U>& );
    ~VectorVariable();
    
    double&        operator()( size_t i );
    const double&  operator()( size_t i ) const;
    double         operator[]( size_t i ) const;
    void             Component( size_t, double );
    double         Component( size_t i ) const;

    VectorVariable   operator+( double val ) const; 
    VectorVariable   operator-( double val ) const;
    VectorVariable   operator*( double val ) const; 
    VectorVariable   operator/( double val ) const;
    VectorVariable   operator^( double val ) const;

    VectorVariable&  operator+=( double val );
    VectorVariable&  operator-=( double val );
    VectorVariable&  operator*=( double val );
    VectorVariable&  operator/=( double val );
    
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
    
    VectorVariable&  operator=( double );
    VectorVariable&  operator=( const Point<2U>& );
    VectorVariable&  operator=( const ScalarVariable& );
    VectorVariable&  operator=( const VectorVariable& );
    VectorVariable&  operator=( VectorVariable&& ) = default;

    // extra operators
    bool             operator==( const VectorVariable& ) const;
    bool             operator!=( const VectorVariable& ) const;
    // compare length
    bool             operator<(  const VectorVariable& ) const;
  
    // Normal Methods
    VARIABLE_FLAG&   Flag( const size_t& i=0 );
    VARIABLE_FLAG    Flag( const size_t& i=0 ) const;
    size_t           Size() const;
    void             Resize( size_t newSize, double newValue = std::numeric_limits<double>::quiet_NaN() );
    double         Length() const;
    double         AngleTo( const VectorVariable& v ) const;
    Point<2U>        P() const;
    bool             IsWithinRange( double vmin, double vmax ) const;
    VectorVariable   Flip();
    VectorVariable   ProjectOnto( const std::vector<double>& v ) const;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const;
    void             Invert();
  
    void			       EuclideanNormalize();
    double         DotProduct( const csmp::Point<2U>& p ) const;
    double         DotProduct( const VectorVariable& v ) const;
    VectorVariable   CrossProduct( const csmp::Point<2U>& p ) const;
    VectorVariable   CrossProduct( const VectorVariable& v ) const;

    void             In();
    void             Out() const;
    bool             In( std::fstream& fp );
    bool             Out( std::fstream& fp ) const;

    friend class TensorVariable<2U>;

  private:
    std::array<VARIABLE_FLAG,2U> flag;
    std::array<double,2U>      data;
};


// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts

} // end namespace csmp

#endif












