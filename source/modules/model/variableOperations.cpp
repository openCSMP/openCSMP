#include "variableOperations.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
double dotProduct( const vector<double>& v1, const vector<double>& v2 )
  {
      assert( !v1.empty() );
      assert( !v2.empty() );
      assert( v1.size() == v2.size() );
      double value(0.);
      for ( auto i=0U; i<dim; ++i ) value += v1[i] * v2[i];
      return value;
  }

template double dotProduct<1U>( const vector<double>&, const vector<double>& );
template double dotProduct<2U>( const vector<double>&, const vector<double>& );
template double dotProduct<3U>( const vector<double>&, const vector<double>& );




// addition
template<>
Point<1U> operator+(const Point<1U>& p, const VectorVariable<1U>& v ) {
     return Point<1U>(p[0]+v[0]);
  }

template<>
Point<2U> operator+(const Point<2U>& p, const VectorVariable<2U>& v ) {
     return Point<2U>(p[0]+v[0],p[1]+v[1]);
  }

template<>
Point<3U> operator+(const Point<3U>& p, const VectorVariable<3U>& v ) {
     return Point<3U>(p[0]+v[0],p[1]+v[1],p[2]+v[2]);
  }

// subtraction
template<>
Point<1U> operator-(const Point<1U>& p, const VectorVariable<1U>& v ) {
     return Point<1U>(p[0]-v[0]);
  }

template<>
Point<2U> operator-(const Point<2U>& p, const VectorVariable<2U>& v ) {
     return Point<2U>(p[0]-v[0],p[1]-v[1]);
  }

template<>
Point<3U> operator-(const Point<3U>& p, const VectorVariable<3U>& v ) {
     return Point<3U>(p[0]-v[0],p[1]-v[1],p[2]-v[2]);
  }


// multiplication
template<>
Point<1U> operator*(const Point<1U>& p, const VectorVariable<1U>& v ) {
     return Point<1U>(p[0]*v[0]);
  }

template<>
Point<2U> operator*(const Point<2U>& p, const VectorVariable<2U>& v ) {
     return Point<2U>(p[0]*v[0],p[1]*v[1]);
  }

template<>
Point<3U> operator*(const Point<3U>& p, const VectorVariable<3U>& v ) {
     return Point<3U>(p[0]*v[0],p[1]*v[1],p[2]*v[2]);
  }


// division result = p / v
template<>
Point<1U> operator/(const Point<1U>& p, const VectorVariable<1U>& v ) {
     return Point<1U>(p[0]/v[0]);
  }

template<>
Point<2U> operator/(const Point<2U>& p, const VectorVariable<2U>& v ) {
     return Point<2U>(p[0]/v[0],p[1]/v[1]);
  }

template<>
Point<3U> operator/(const Point<3U>& p, const VectorVariable<3U>& v ) {
     return Point<3U>(p[0]/v[0],p[1]/v[1],p[2]/v[2]);
  }


// dot products: result = v . p
template<>
double dotProduct( const VectorVariable<1>& v,  const csmp::Point<1>& p ) {
   return std::move(v[0]*p[0]);
}

template<>
double dotProduct( const VectorVariable<2>& v,  const csmp::Point<2>& p ) {
   return std::move(v[0]*p[0] + v[1]*p[1]);
}

template<>
double dotProduct( const VectorVariable<3>& v,  const csmp::Point<3>& p ) {
   return std::move(v[0]*p[0] + v[1]*p[1] + v[2]*p[2]);
}


// dot products vector variables: result = v1 . v2
template<>
double dotProduct( const VectorVariable<1>& v1, const VectorVariable<1>& v2 ) {
   return std::move(v1[0]*v2[0]);
}

template<>
double dotProduct( const VectorVariable<2>& v1, const VectorVariable<2>& v2 ) {
   return std::move(v1[0]*v2[0] + v1[1]*v2[1]);
}

template<>
double dotProduct( const VectorVariable<3>& v1, const VectorVariable<3>& v2 ) {
   return std::move(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}



VectorVariable<1U> crossProduct(  const VectorVariable<1U>& v,  const csmp::Point<1U>& )
{
    return VectorVariable<1U>(v.Flag(),0.);
}

VectorVariable<2U> crossProduct(  const VectorVariable<2U>& v1,  const csmp::Point<2U>& p )
{
    return VectorVariable<2U>(v1.Flag(0),v1.Flag(1),0.0,v1[0]*p[1]-v1[1]*p[0]);
}

VectorVariable<3U> crossProduct(  const VectorVariable<3U>& v1,  const csmp::Point<3U>& p ) {

    return VectorVariable<3U>( v1.Flag(0), v1.Flag(1), v1.Flag(2),
                                        (v1[1]*p[2]-v1[2]*p[1]),
                                        (v1[2]*p[0] - v1[0]*p[2]),
                                        (v1[0]*p[1] - v1[1]*p[0]));
}


VectorVariable<1U> crossProduct(  const VectorVariable<1U>& v1,  const VectorVariable<1U>& v2 ) {
    return VectorVariable<1U>(v1.Flag(),0.0);
}

VectorVariable<2U> crossProduct(  const VectorVariable<2U>& v1,  const VectorVariable<2U>& v2 ) {
    return VectorVariable<2U>(v1.Flag(0),v1.Flag(1),0.0,v1[0]*v2[1]-v1[1]*v2[0]);
}

VectorVariable<3U> crossProduct(  const VectorVariable<3U>& v1,  const VectorVariable<3U>& v2 ) {

    return VectorVariable<3U>( v1.Flag(0), v1.Flag(1), v1.Flag(2),
                                         (v1[1]*v2[2]-v1[2]*v2[1]),(v1[2]*v2[0] - v1[0]*v2[2]),
                                         (v1[0]*v2[1] - v1[1]*v2[0]));
}

/**
    @note the flags being kept belong to the vector (legacy).
    perhaps, because of the multiplication order, the tensor ones should be kept.
    Nevertheless, I will keep this since the TensorVariable_Test uses this convention

    @author Julian Mindel 
    @date 20.04.2016
*/
template<uint32_t dim>
VectorVariable<dim> multiplyTensorByVector( const TensorVariable<dim>& ts, const VectorVariable<dim>& vc )
{
    VectorVariable<dim> resultvec(ANY,0.0);

    for (auto i = 0U ; i < dim; ++i)
        for (auto j = 0U ; j < dim; ++j)
            resultvec(i)+=ts(i,j) * vc[j];

    for (auto i = 0U ; i < dim; ++i)
        resultvec.Flag(i)=ts.Flag(i);
  
    return std::move(resultvec);
}

/**
    this essentially carries out the same task as MultiplyTensorByVector, but has a meaningful name
    emphasizing the fact that only vectors that are horizontal/transposed can be multiplied by tensors.
*/
template<uint32_t dim>
VectorVariable<dim> multiplyHorizontalVectorByTensor( const VectorVariable<dim>& vc, const TensorVariable<dim>& ts )
{
    VectorVariable<dim> resultvec(ANY,0.0);

    for (auto i = 0U ; i < dim; ++i)
        for (auto j = 0U ; j < dim; ++j)
            resultvec(i)+=vc[j] * ts(j,i);

    for (auto i = 0U ; i < dim; ++i)
        resultvec.Flag(i)=vc.Flag(i);

    return std::move(resultvec);
}

template<uint32_t dim>
TensorVariable<dim> multiplyTensorByTensor( const TensorVariable<dim>& ts, const TensorVariable<dim>& ts2 )
{
    TensorVariable<dim> resulttensor(ANY,0.0);

    for (auto i = 0U ; i < dim; ++i)
        for (auto j = 0U ; j < dim; ++j)
            for (auto k = 0U ; k < dim; ++k)
                resulttensor(i,j)+=ts(i,k) * ts2(k,j);
    for (auto i = 0U ; i < dim; ++i)
        resulttensor.Flag(i)=ts.Flag(i);

    return std::move(resulttensor);
}



template<uint32_t>
bool isDiagonalTensor( const TensorVariable<1U>& ts ) { return true; }

/// 2D version
template<uint32_t>
bool isDiagonalTensor( const TensorVariable<2U>& ts )
{
  // if the off-diagonal elements are numerically zero
  if ( ts( 0, 1 ) <= numeric_limits<double>::epsilon() &&
       ts( 1, 0 ) <= numeric_limits<double>::epsilon() ) return true;
  return false;
}

/// 3D version
template<uint32_t>
bool isDiagonalTensor( const TensorVariable<3U>& ts )
{
  // if the off-diagonal elements are numerically zero
  if ( ts( 0, 1 ) <= numeric_limits<double>::epsilon() &&
       ts( 0, 2 ) <= numeric_limits<double>::epsilon() &&
       ts( 1, 2 ) <= numeric_limits<double>::epsilon() &&
       ts( 1, 0 ) <= numeric_limits<double>::epsilon() &&
       ts( 2, 0 ) <= numeric_limits<double>::epsilon() &&
       ts( 2, 1 ) <= numeric_limits<double>::epsilon() ) return true;
  return false;
}

// this wants to be a lambda function in the next method
/// recovering and sorting to find minimum and maximum Eigen values
template<uint32_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double& tmin, double& tmax )
{
  VectorVariable<dim>  evals;

  // doing simple case first
  if ( isDiagonalTensor( ts ) ) for ( auto i{0U}; i<dim; ++i ) evals( i ) = ts( i, i );
  else ts.EigenValues( evals );
  std::set<double> min_max;
  for ( auto i{0U}; i<dim; i++ ) min_max.insert( evals[i] );
  tmin = (*min_max.begin());
  tmax = (*min_max.rbegin());
}




/// return angle in degrees
/// this function was taken from an old one called AngleTo, and templetized
template<uint32_t dim>
double  angleBetween( const VectorVariable<dim>& v1, const VectorVariable<dim>& v2 )
{
    double ab, a_dot_b;

    // a b
    // ---
    ab      = dotProduct(v1,v2);
    // |a| . |b|
    // ---------
    a_dot_b = v1.Length()*v2.Length();

    // a b
    // ---
    double cos_angle = ab / a_dot_b;

    // if zero intercept
    if ( cos_angle == 0.0 ) return  90.0;
    // if outside of range of 'acos' function
    if ( cos_angle >  1.0 ) return   0.0;
    if ( cos_angle < -1.0 ) return 180.0;

    return (180.0/3.14159265358979324) * std::acos(cos_angle);
}



template double  angleBetween( const VectorVariable<1U> & v1, const VectorVariable<1U>& v2 );
template double  angleBetween( const VectorVariable<2U> & v1, const VectorVariable<2U>& v2 );
template double  angleBetween( const VectorVariable<3U> & v1, const VectorVariable<3U>& v2 );


// TODO: this is what operator overloading should be for; revise!
template VectorVariable<1U> multiplyTensorByVector( const TensorVariable<1U>&, const VectorVariable<1U>& ) ;
template VectorVariable<2U> multiplyTensorByVector( const TensorVariable<2U>&, const VectorVariable<2U>& ) ;
template VectorVariable<3U> multiplyTensorByVector( const TensorVariable<3U>&, const VectorVariable<3U>& ) ;

template VectorVariable<1U> multiplyHorizontalVectorByTensor( const VectorVariable<1U>&, const TensorVariable<1U>&  ) ;
template VectorVariable<2U> multiplyHorizontalVectorByTensor( const VectorVariable<2U>&, const TensorVariable<2U>&  ) ;
template VectorVariable<3U> multiplyHorizontalVectorByTensor( const VectorVariable<3U>&, const TensorVariable<3U>&  ) ;

template TensorVariable<1U> multiplyTensorByTensor( const TensorVariable<1U>&, const TensorVariable<1U>& ) ;
template TensorVariable<2U> multiplyTensorByTensor( const TensorVariable<2U>&, const TensorVariable<2U>& ) ;
template TensorVariable<3U> multiplyTensorByTensor( const TensorVariable<3U>&, const TensorVariable<3U>& ) ;

} // end namespace csmp

