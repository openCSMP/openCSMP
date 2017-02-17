#include "VectorVariable2.h"

using namespace std;

namespace csmp {

VectorVariable<2U>::VectorVariable()
  : flag{{ANY,ANY}},
    data{{std::numeric_limits<double64>::quiet_NaN(),std::numeric_limits<double64>::quiet_NaN()}}
  { 
  }

  
VectorVariable<2U>::VectorVariable( const VectorVariable<2U>& v )
  : flag(v.flag),
    data(v.data)
 {
 }


VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f, double64 val )
  : flag{{f,f}},
    data{{val,val}}
 {
 }


VectorVariable<2U>::VectorVariable( VARIABLE_FLAG f1, VARIABLE_FLAG f2, double64 val1, double64 val2 )
  : flag{{f1,f2}},
    data{{val1,val2}}
 {
 }


VectorVariable<2U>::VectorVariable( const std::vector<double64>& v )
  : flag{{ANY,ANY}},
    data{{v[0],v[1]}}
 {
 }


VectorVariable<2U>::VectorVariable( const Point<2U>& p )
  : flag{{ANY,ANY}},
    data{{p[0],p[1]}}
 {
 }


VectorVariable<2U>::~VectorVariable() 
  {
  }


VectorVariable<2U>&  VectorVariable<2U>::operator=( const VectorVariable<2U>& v )
 {
    if ( &v != this ) {
          flag = v.flag;
          data = v.data;
       } 
    return *this; 
 }


double64& VectorVariable<2U>::operator()( size_t i )       
  { 
     if ( i==0U ) return data[0];
     return              data[1]; 
  }


const double64& VectorVariable<2U>::operator()( size_t i ) const 
  { 
     if ( i==0U ) return data[0];
     return              data[1]; 
  }



double64  VectorVariable<2U>::operator[]( size_t i ) const 
 { 
    if ( i==0U ) return data[0];
    return              data[1]; 
 }



void  VectorVariable<2U>::Component( size_t i, double64 val ) 
 { 
    if ( i==0U ) data[0] = val;
    else         data[1] = val; 
 }



double64  VectorVariable<2U>::Component( size_t i ) const 
 { 
    if ( i==0U ) return data[0];
    return              data[1]; 
 }


size_t VectorVariable<2U>::Size() const
  {
    return 2U;
  }

void VectorVariable<2U>::Resize( size_t, double64 newValue )
  {
    data[0] = newValue;
    data[1] = newValue;
  }

VectorVariable<2U>  VectorVariable<2U>::operator+( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] + v.data[0], data[1] + v.data[1] ));
 }




VectorVariable<2U>  VectorVariable<2U>::operator-( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] - v.data[0], data[1] - v.data[1] ));
 }



VectorVariable<2U>  VectorVariable<2U>::operator*( const VectorVariable<2U>& v ) const 
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] * v.data[0], data[1] * v.data[1] ));
 } 



VectorVariable<2U>  VectorVariable<2U>::operator/( const VectorVariable<2U>& v ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] / v.data[0], data[1] / v.data[1] ));
 }



VectorVariable<2U>  VectorVariable<2U>::operator+( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] + val, data[1] + val ));
 }



VectorVariable<2U>  VectorVariable<2U>::operator-( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] - val, data[1] - val ));
 }



VectorVariable<2U>  VectorVariable<2U>::operator*( double64 val ) const 
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] * val, data[1] * val ));
 } 



VectorVariable<2U>  VectorVariable<2U>::operator/( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], data[0] / val, data[1] / val ));
 }




VectorVariable<2U>  VectorVariable<2U>::operator^( double64 val ) const
 {
    return std::move(VectorVariable( flag[0], flag[1], std::pow( data[0],val ), std::pow( data[1], val ) ));
 }




VectorVariable<2U>&  VectorVariable<2U>::operator+=( double64 val )
 {
    data[0] += val;
    data[1] += val;
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator-=( double64 val )
 {
    data[0] -= val;
    data[1] -= val;
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator*=( double64 val )
 {
    data[0] *= val;
    data[1] *= val;
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator/=( double64 val )
 {
    data[0] /= val;
    data[1] /= val;
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator+=( const ScalarVariable& sc )
 {
    data[0] += sc();
    data[1] += sc();
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator-=( const ScalarVariable& sc )
 {
    data[0] -= sc();
    data[1] -= sc();
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator*=( const ScalarVariable& sc )
 {
    data[0] *= sc();
    data[1] *= sc();
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator/=( const ScalarVariable& sc )
 {
    data[0] /= sc();
    data[1] /= sc();
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator+=( const VectorVariable<2U>& v )
 {
    data[0] += v.data[0];
    data[1] += v.data[1];
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator-=( const VectorVariable<2U>& v )
 {
    data[0] -= v.data[0];
    data[1] -= v.data[1];
        
    return *this; 
 }




VectorVariable<2U>&  VectorVariable<2U>::operator*=( const VectorVariable<2U>& v )
 {
    data[0] *= v.data[0];
    data[1] *= v.data[1];
        
    return *this;  
 }




VectorVariable<2U>&  VectorVariable<2U>::operator/=( const VectorVariable<2U>& v )
 {
    data[0] /= v.data[0];
    data[1] /= v.data[1];
        
    return *this; 
 }




// --------------------
// ASSIGNMENT OPERATORS
// --------------------

VectorVariable<2U>&  VectorVariable<2U>::operator=( double64 val )
 {
    data[0] = val;
    data[1] = val;
        
    return *this; 
 }



VectorVariable<2U>&  VectorVariable<2U>::operator=( const csmp::Point<2U>& p )
 {
    data[0] = p[0];
    data[1] = p[1];
        
    return *this; 
 }



VectorVariable<2U>&  VectorVariable<2U>::operator=( const ScalarVariable& sc )
 {
    flag[0] = flag[1] = sc.Flag();
    data[0] = data[1] = sc();
       
    return *this; 
 }


bool VectorVariable<2U>::operator==( const VectorVariable<2U>& v ) const
 {
    return( flag == v.flag && data == v.data );
 }


bool VectorVariable<2U>::operator!=( const VectorVariable<2U>& v ) const
 {
    return( flag != v.flag || data != v.data );
 }


// compare the length of two vectors

bool VectorVariable<2U>::operator<( const VectorVariable<2U>& v ) const
 {
    return (this < &v);
 } 





// -------
// METHODS
// -------

/// L2 norm
void VectorVariable<2U>::EuclideanNormalize() 
 {
    const double64 fNorm(std::sqrt(data[0]*data[0] + data[1]*data[1]));
    
    if(fNorm == 0.) return; //added AP
    
    data[0] /= fNorm;
    data[1] /= fNorm;
 } 



double64 VectorVariable<2U>::DotProduct( const csmp::Point<2U>& p ) const
 {
    return data[0] * p[0] + data[1] * p[1]; 
 }

double64 VectorVariable<2U>::DotProduct( const VectorVariable& v ) const
 {
    return data[0] * v[0] + data[1] * v[1];
 }


VectorVariable<2U> VectorVariable<2U>::CrossProduct( const csmp::Point<2U>& p ) const
 {
     return std::move(VectorVariable<2U>( flag[0], flag[1], 0., data[0]*p[1] - data[1]*p[0] ));
 }

VectorVariable<2U> VectorVariable<2U>::CrossProduct( const VectorVariable& v ) const
 {
     return std::move(VectorVariable<2U>( flag[0], flag[1], 0., data[0]*v[1] - data[1]*v[0] ));
 }


double64  VectorVariable<2U>::Length() const 
 {
    return std::sqrt( data[0]*data[0] + data[1]*data[1] );
 }



VARIABLE_FLAG&  VectorVariable<2U>::Flag( const size_t& i )       
 { 
    if ( i==0U ) return flag[0];
    return flag[1]; 
 }


VARIABLE_FLAG  VectorVariable<2U>::Flag( const size_t& i ) const 
 { 
    if ( i==0U ) return flag[0];
    return flag[1]; 
 }


Point<2U>  VectorVariable<2U>::P() const
 {
    return csmp::Point<2U>(data[0],data[1]);
 }



bool  VectorVariable<2U>::IsWithinRange( double64 vmin, double64 vmax ) const
 {
    if ( data[0] < vmin || data[0] > vmax ) return false;
    if ( data[1] < vmin || data[1] > vmax ) return false;
     
    return true;
 }


 bool VectorVariable<2U>::Out( FILE* fp ) const
  {
     fwrite( (void*)this, sizeof(VectorVariable<2U>), 1, fp );
     return true;
  }

 bool VectorVariable<2U>::In( FILE* fp )
  {
     fread( (void*)this, sizeof(VectorVariable<2U>), 1, fp );
     return true;
  }

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




void  VectorVariable<2U>::Out() const 
 {
     string  status;
     
     cout <<"\nStatus: "<< endl; 
     for ( size_t i=0; i<2U; i++ ) 
       cout << (status = parseStatus(flag[i])) <<"\t\t";
     cout << endl;
     for ( size_t i=0; i<2U; i++ ) cout << data[i] <<"\t\t";
     cout << endl;

 } // end Out


} // end namespace csmp

