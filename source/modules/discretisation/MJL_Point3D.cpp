#include "MJL_Point3D.h"

using namespace std;

namespace mjl {

bool Point3D::CoincidesWithWithinTolerance( const Point3D& p,
                                            double tolerance ) const
 {
    if ( (x_-tolerance <= p.x_ && x_+tolerance >= p.x_) &&
         (y_-tolerance <= p.y_ && y_+tolerance >= p.y_) &&
         (z_-tolerance <= p.x_ && z_+tolerance >= p.z_) ) return true;
         
    return false;
 }





void Point3D::Out() const
 {
    cout <<"\nMJL_Point3D::Out: "<< x_ <<"  "<< y_ <<"  "<< z_ << endl;
 }
 
} // end mjl

