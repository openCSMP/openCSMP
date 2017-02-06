#include "VectorVariable2.h"

using namespace std;

namespace csmp {


/** return angle in degrees

*/
double64  VectorVariable<2U>::AngleTo( const VectorVariable<2U>& v ) const
 {
    double64 ab, a_dot_b;
 
     // a b
     // ---
     ab      = data[0]*v.data[0] + data[1]*v.data[1];
     // |a| . |b|
     // ---------
     a_dot_b = std::sqrt( (data[0]*data[0]+data[1]*data[1]) * 
                         (v.data[0]*v.data[0]+v.data[1]*v.data[1]) );
      
    // a b
    // ---
    double64 cos_angle = ab / a_dot_b;

    // if zero intercept
    if ( cos_angle == 0.0 ) return  90.0;
    // if outside of range of 'acos' function
    if ( cos_angle >  1.0 ) return   0.0;
    if ( cos_angle < -1.0 ) return 180.0;
        
    return (180.0/3.14159265358979324) * std::acos(cos_angle);
 }    



VectorVariable<2U>  VectorVariable<2U>::Flip()
 {
    VectorVariable<2U>  temp;
    
    temp.flag[0] = flag[1];
    temp.flag[1] = flag[0];
    temp.data[0] = data[1];
    temp.data[1] = data[0];

    return std::move(temp);
 }

/// Multiplies by negative unity vector
void  VectorVariable<2U>::Invert()
  {
    data[0] *= -1.;
    data[1] *= -1.;
  }

VectorVariable<2U>  VectorVariable<2U>::ProjectOnto( const std::vector<double64>& v ) const
 {
    double64 ratio((data[0]*v[0] + data[1]*v[1]) / (v[0]*v[0] + v[1]*v[1]));

    return std::move(VectorVariable<2U>( flag[0], flag[1],
                               v[0]*ratio, v[1]*ratio ));
 }



VectorVariable<2U>  VectorVariable<2U>::ProjectOnto( const VectorVariable<2U>& v ) const
 {
    double64 ratio((data[0]*v.data[0] + data[1]*v.data[1]) / (v.data[0]*v.data[0] + v.data[1]*v.data[1]));

    return std::move(VectorVariable<2U>( flag[0], flag[1],
                               v.data[0]*ratio, v.data[1]*ratio ));
 }



void  VectorVariable<2U>::In()
 {
     string  status;
     
     cout.flush();
     for ( size_t i=0; i<2U; i++ )
       {
          if ( i == 0 ) cout <<"\nEnter status for x-component of variable: ";
          else          cout <<"\nEnter status for y-component of variable: ";
          cout.flush();
          cin >> status;
          flag[i] = parseStatus( status.c_str() );
       }
     
     cout <<"\nEnter x=0 and y=1 vector variable elements: ";
     cout.flush();
     for ( size_t i=0; i<2U; i++ ) cin >> data[i];

 } // end In




void  VectorVariable<2U>::Out(std::ostream& os) const 
 {
     string  status;
     
     os <<"\nStatus: "<< endl; 
     for ( size_t i=0; i<2U; i++ ) 
       os << (status = parseStatus(flag[i])) <<"\t\t";
     os << endl;
     for ( size_t i=0; i<2U; i++ ) os << data[i] <<"\t\t";
     os << endl;

 } // end Out


} // end namespace csmp

