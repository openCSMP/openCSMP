#ifndef CSMP_VECTOR_VARIABLE1_H
#define CSMP_VECTOR_VARIABLE1_H

#include "Point.h"
#include "ScalarVariable.h"

namespace csmp {

/// @brief 1D full specialization of vector-type variable template (other vectors are declared in files #1, #2)
template<>
class VectorVariable<1U> {
  public:
    static constexpr VARIABLE_TYPE VariableType = VECTOR;

    VectorVariable();

    VectorVariable( VARIABLE_FLAG f, double val );
    explicit VectorVariable( const std::vector<double>& );
    explicit VectorVariable( const csmp::Point<1U>& );
    
    // access of vector elements
    double&        operator()( uint32_t i );
    const double&  operator()( uint32_t i ) const;
    double         operator[]( uint32_t i ) const;
    void           Component( uint32_t, double );
    double         Component( uint32_t i ) const;

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
    
    VectorVariable&  operator=( double val );
    VectorVariable&  operator=( const csmp::Point<1U>& );
    VectorVariable&  operator=(  const ScalarVariable& );

    // extra operators
    bool             operator==( const VectorVariable& ) const;
    bool             operator!=( const VectorVariable& ) const;
    // compare length
    bool             operator<( const VectorVariable& ) const;
  
    // Normal Methods
    VARIABLE_FLAG&   Flag( uint32_t i=0 );
    VARIABLE_FLAG    Flag( uint32_t i=0 ) const;
    uint32_t         Size() const;
  
    double           Length() const;
    Point<1U>        P() const;
    bool             IsWithinRange( double vmin, double vmax ) const;
    bool             Has_NaN_Values() const { return std::isnan(data); }

    VectorVariable   Flip();
    void             Invert();

    double           DotProduct( const csmp::Point<1U>& p ) const;
    double           DotProduct( const VectorVariable& v ) const;
    VectorVariable   CrossProduct( const csmp::Point<1U>& p ) const;
    VectorVariable   CrossProduct( const VectorVariable& v ) const;
    VectorVariable   ProjectOnto( const std::vector<double>& v ) const;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const;
  
    /// normalize length to 1
    void             EuclideanNormalize() { data = 1.; }

    void             In();
    void             Out() const;
    bool             In( std::fstream& fp );
    bool             Out( std::fstream& fp ) const;
    
    friend class TensorVariable<1U>;

  private:
    VARIABLE_FLAG flag;
    double        data;
};

// copyright (c) 2001 by S.K. Matthai, S. Geiger & Stephen G. Roberts




} // end namespace csmp

#endif












