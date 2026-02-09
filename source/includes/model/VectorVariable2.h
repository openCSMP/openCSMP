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

    VectorVariable() noexcept;

    VectorVariable( VARIABLE_FLAG f, double val ) noexcept;
    VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double val1, double val2 ) noexcept;
    explicit VectorVariable( const std::vector<double>& ) noexcept;
    explicit VectorVariable( const csmp::Point<2U>& ) noexcept;
    
    double&        operator()( uint32_t i ) noexcept;
    const double&  operator()( uint32_t i ) const noexcept;
    double         operator[]( uint32_t i ) const noexcept;
    void           Component( uint32_t, double ) noexcept;
    double         Component( uint32_t i ) const noexcept;

    VectorVariable   operator^( double val ) const noexcept;

    VectorVariable&  operator+=( double val ) noexcept;
    VectorVariable&  operator-=( double val ) noexcept;
    VectorVariable&  operator*=( double val ) noexcept;
    VectorVariable&  operator/=( double val ) noexcept;
    
    VectorVariable&  operator+=( const ScalarVariable& ) noexcept;
    VectorVariable&  operator-=( const ScalarVariable& ) noexcept;
    VectorVariable&  operator*=( const ScalarVariable& ) noexcept;
    VectorVariable&  operator/=( const ScalarVariable& ) noexcept;
    
    VectorVariable   operator+(  const VectorVariable& ) const noexcept;
    VectorVariable   operator-(  const VectorVariable& ) const noexcept;
    VectorVariable   operator*(  const VectorVariable& ) const noexcept;
    VectorVariable   operator/(  const VectorVariable& ) const noexcept;

    VectorVariable&  operator+=( const VectorVariable& ) noexcept;
    VectorVariable&  operator-=( const VectorVariable& ) noexcept;
    VectorVariable&  operator*=( const VectorVariable& ) noexcept;
    VectorVariable&  operator/=( const VectorVariable& ) noexcept;
    
    VectorVariable&  operator=( double ) noexcept;
    VectorVariable&  operator=( const Point<2U>& ) noexcept;
    VectorVariable&  operator=( const ScalarVariable& ) noexcept;

    // extra operators
    bool             operator==( const VectorVariable& ) const noexcept;
    bool             operator!=( const VectorVariable& ) const noexcept;
    // compare length
    bool             operator<(  const VectorVariable& ) const noexcept;
  
    // Normal Methods
    VARIABLE_FLAG&   Flag( uint32_t i=0 ) noexcept;
    VARIABLE_FLAG    Flag( uint32_t i=0 ) const noexcept;

    static constexpr uint32_t Size() noexcept { return 2u; };

    double           Length() const noexcept;
    double           AngleTo( const VectorVariable& v ) const noexcept;
    Point<2U>        P() const noexcept;
    bool             IsWithinRange( double vmin, double vmax ) const noexcept;
    bool             Has_NaN_Values() const noexcept;
    VectorVariable   Flip() noexcept;
    VectorVariable   ProjectOnto( const std::vector<double>& v ) const noexcept;
    VectorVariable   ProjectOnto( const VectorVariable& v ) const noexcept;
    void             Invert() noexcept;
  
    void			       EuclideanNormalize() noexcept;
    double           DotProduct( const csmp::Point<2U>& p ) const noexcept;
    double           DotProduct( const VectorVariable& v ) const noexcept;
    VectorVariable   CrossProduct( const csmp::Point<2U>& p ) const noexcept;
    VectorVariable   CrossProduct( const VectorVariable& v ) const noexcept;

    void             In();
    void             Out() const noexcept;
    bool             In( std::fstream& fp );
    bool             Out( std::fstream& fp ) const;

    friend class TensorVariable<2U>;

  private:
    std::array<VARIABLE_FLAG,2U> flag;
    std::array<double,2U>        data;
};


// copyright (c) 2001 by S.K. Matthaei, S. Geiger & Stephen G. Roberts

} // end namespace csmp

#endif












