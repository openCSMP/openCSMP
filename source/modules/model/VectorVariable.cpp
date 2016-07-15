#include "VectorVariable.h"

using namespace std;

namespace csmp {


VectorVariable<3U>::VectorVariable()
  : flag{{ANY,ANY,ANY}},
    data{{std::numeric_limits<double64>::quiet_NaN(),std::numeric_limits<double64>::quiet_NaN(),std::numeric_limits<double64>::quiet_NaN()}}
  { 
  }
  




VectorVariable<3U>::VectorVariable( VARIABLE_FLAG f, double64 val )
  : flag{{f,f,f}},
    data{{val,val,val}}
 {
 }




VectorVariable<3U>::VectorVariable( VARIABLE_FLAG f1,   VARIABLE_FLAG f2, VARIABLE_FLAG f3, 
                                    double64  val1, double64  val2, double64  val3 )
  : flag{{f1,f2,f3}},
    data{{val1,val2,val3}}
 {
 }





VectorVariable<3U>&  VectorVariable<3U>::operator=( const VectorVariable<3U>& v )
 {
    if ( &v != this ) {
          flag = v.flag;
          data = v.data;
       } 
    return *this; 
 }


VectorVariable<3U>&  VectorVariable<3U>::operator=( VectorVariable<3U>&& v )
 {
    if ( &v != this ) {
          flag = {v.flag};
          data = {v.data};
       } 
    return *this; 
 }




VectorVariable<3U>  VectorVariable<3U>::operator+( const VectorVariable<3U>& v ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2], 
                           data[0] + v.data[0], data[1] + v.data[1], data[2] + v.data[2] );
 }




VectorVariable<3U>  VectorVariable<3U>::operator-( const VectorVariable<3U>& v ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] - v.data[0], data[1] - v.data[1], data[2] - v.data[2] );
 }



VectorVariable<3U>  VectorVariable<3U>::operator*( const VectorVariable<3U>& v ) const 
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] * v.data[0], data[1] * v.data[1], data[2] * v.data[2] );
 } 



VectorVariable<3U>  VectorVariable<3U>::operator/( const VectorVariable<3U>& v ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] / v.data[0], data[1] / v.data[1], data[2] / v.data[2] );
 }



VectorVariable<3U>  VectorVariable<3U>::operator+( double64 val ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] + val, data[1] + val, data[2] + val );
 }



VectorVariable<3U>  VectorVariable<3U>::operator-( double64 val ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] - val, data[1] - val, data[2] - val );
 }



VectorVariable<3U>  VectorVariable<3U>::operator*( double64 val ) const 
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] * val, data[1] * val, data[2] * val );
 } 



VectorVariable<3U>  VectorVariable<3U>::operator/( double64 val ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[0] / val, data[1] / val, data[2] / val );
 }




VectorVariable<3U>  VectorVariable<3U>::operator^( double64 val ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           std::pow( data[0],val ), std::pow( data[1], val ), std::pow( data[2], val ) );
 }




/// cross product ' % ' of two vectors
VectorVariable<3U>  VectorVariable<3U>::operator%( const VectorVariable<3U>& v ) const
 {
    return VectorVariable( flag[0], flag[1], flag[2],
                           data[1]*v.data[2] - v.data[1]*data[2],
                         -(data[0]*v.data[2] - v.data[0]*data[2]),
                           data[0]*v.data[1] - v.data[0]*data[1] );
 } 




/// return angle in degrees
double64  VectorVariable<3U>::AngleTo( const VectorVariable<3U>& v ) const 
 {
    double64 ab      = data[0]*v.data[0]+data[1]*v.data[1]+data[2]*v.data[2];
    double64 a_dot_b = std::sqrt( (data[0]*data[0]+data[1]*data[1]+data[2]*data[2]) * 
                        (v.data[0]*v.data[0]+v.data[1]*v.data[1]+v.data[2]*v.data[2]) );
    // a b
    // ---
    double64 cos_angle = ab / a_dot_b;

    // if zero intercept
    if ( cos_angle == 0.0 ) return  90.0;
    // if outside of range of 'acos' function
    if ( cos_angle >  1.0 ) return   0.0;
    if ( cos_angle < -1.0 ) return 180.0;
        
    return (static_cast<double64>(180.)/static_cast<double64>(3.14159265358979324)) * std::acos(cos_angle);
 }    



VectorVariable<3U>  VectorVariable<3U>::Flip()
 {
    VectorVariable<3U>  temp;
    
    temp.flag[0] = flag[2];
    temp.flag[2] = flag[0];
    temp.flag[1] = flag[1];
    temp.data[0] = data[2];
    temp.data[2] = data[0];
    temp.data[1] = data[1];

    return temp; 
 }


/// Multiplies by negative unity vector
void  VectorVariable<3U>::Invert()
  {
    data[0] *= -1.;
    data[1] *= -1.;
    data[2] *= -1.;
  }




/// projects this vector variable onto the supplied vector v and returns the projection vector

VectorVariable<3U>  VectorVariable<3U>::ProjectOnto( const std::vector<double64>& v ) const
 {
    double64 ratio((data[0]*v[0] + data[1]*v[1] + data[2]*v[2]) / (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));

    return VectorVariable<3U>( flag[0], flag[1], flag[2], v[0]*ratio, v[1]*ratio, v[2]*ratio );
 }



VectorVariable<3U>  VectorVariable<3U>::ProjectOnto( const VectorVariable<3U>& v ) const
 {
    double64 ratio((data[0]*v.data[0] + data[1]*v.data[1] + data[2]*v.data[2]) / 
                    (v.data[0]*v.data[0] + v.data[1]*v.data[1] + v.data[2]*v.data[2]));

    return VectorVariable<3U>( flag[0], flag[1], flag[2],
                               v.data[0]*ratio, v.data[1]*ratio, v.data[2]*ratio );
 }




void  VectorVariable<3U>::Sqrt() 
 { 
     data[0] = std::sqrt(data[0]);
     data[1] = std::sqrt(data[1]);
     data[2] = std::sqrt(data[2]);
 }




void  VectorVariable<3U>::Ln() 
 { 
     data[0] = std::log(data[0]);
     data[1] = std::log(data[1]);
     data[2] = std::log(data[2]);
 }



void  VectorVariable<3U>::Log10() 
 { 
     data[0] = std::log10(data[0]);
     data[1] = std::log10(data[1]);
     data[2] = std::log10(data[2]); 
 }





template<size_t dim>
ostream&  operator<<( ostream& stream, const VectorVariable<dim>& o )
 {
     for ( size_t i=0U; i<dim; i++ )
       stream << o[i] <<" ("<< parseStatus(o.Flag(i)) <<") ";
       
     return stream;
 }



void  VectorVariable<3U>::In()
 {
     string  status;
     
     for ( size_t i=0; i<3U; i++ ) {
          if ( i == 0 )      cout <<"\nEnter status for x-component of variable: ";
          else if ( i == 1 ) cout <<"\nEnter status for y-component of variable: ";
          else               cout <<"\nEnter status for z-component of variable: ";
          cout.flush();
          cin >> status;
          flag[i] = parseStatus( status.c_str() );
       }
     
     cout <<"\nEnter 3 vector elements: ";
     cout.flush();
     for ( size_t i=0; i<3U; i++ ) cin >> data[i];

 } // end In



void  VectorVariable<3U>::Out() const 
 {
     cout <<"\nStatus: "<< endl; 
     for ( size_t i=0; i<3U; i++ ) 
       cout << parseStatus(flag[i]) <<"\t\t";
     cout << endl;
     for ( size_t i=0; i<3U; i++ ) cout << data[i] <<"\t\t";
     cout << endl;

 } // end Out


template ostream&  operator<< <1U>( ostream& stream, const VectorVariable<1U>& o );
template ostream&  operator<< <2U>( ostream& stream, const VectorVariable<2U>& o );
template ostream&  operator<< <3U>( ostream& stream, const VectorVariable<3U>& o );

} // end namespace csmp

