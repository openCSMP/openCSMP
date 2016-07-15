#include "MJL_Triangle3D.h"
#include "MJL_Edge3D.h"

using namespace std;

namespace mjl {

Triangle3D::Triangle3D( const Point3D& v0, 
                        const Point3D& v1, 
                        const Point3D& v2, 
                        long  id )
  : id_(id), mark_(0),
    v0_(v0), v1_(v1), v2_(v2),
    boundingBox_( Point3D( min3(v0.x_, v1.x_, v2.x_),
                           min3(v0.y_, v1.y_, v2.y_),
                           min3(v0.z_, v1.z_, v2.z_) ),
                  Point3D( max3(v0.x_, v1.x_, v2.x_),
                           max3(v0.y_, v1.y_, v2.y_),
                           max3(v0.z_, v1.z_, v2.z_) ) )
 {
    n_ = crossProduct(v1-v0, v2-v0);
    n_ = (1. / n_.Length()) * n_;
 }



Triangle3D& Triangle3D::operator=( const Triangle3D& t )
 {
    if ( &t != this ) {
         v0_          = t.v0_;
         v1_          = t.v1_;
         v2_          = t.v2_;
         boundingBox_ = t.boundingBox_;
         n_           = t.n_;
         id_          = t.id_; 
         mark_        = t.mark_;
      }
    return *this; 
 }



void Triangle3D::Set( const Point3D& v0, 
                      const Point3D& v1, 
                      const Point3D& v2, 
                      long id )
 {
    id_   = id;
    mark_ = 0;
    v0_   = v0;
    v1_   = v1;
    v2_   = v2;
    boundingBox_.org_.x_  = min3(v0.x_, v1.x_, v2.x_);
    boundingBox_.org_.y_  = min3(v0.y_, v1.y_, v2.y_);
    boundingBox_.org_.z_  = min3(v0.z_, v1.z_, v2.z_);
    boundingBox_.dest_.x_ = max3(v0.x_, v1.x_, v2.x_);
    boundingBox_.dest_.y_ = max3(v0.y_, v1.y_, v2.y_);
    boundingBox_.dest_.z_ = max3(v0.z_, v1.z_, v2.z_);
    n_ = crossProduct(v1-v0, v2-v0);
    n_ = (1.0 / n_.Length()) * n_;
 }




/**
 
Identifies the facing direction of the father triangle, by comparing the
orientation of its normal to the principle coordinate directions. In CSP,
these directions are defined in a righthand coordinate system: 
 
         TOP
       
        y |
    3D:   | 
          |
          o-------> RIGHT
         /       x
     z  / 
     
    FRONT

  

@return the integer value definition from 'MJL_Geometry.h' for the
result of the comparison. If none of the principal axes match, the 
method will return the value MJL3D_TILTED. 

@section implementation Implementation

The comparison will evaluate to one of the principal facing directions if 
the face normal is parallel to-, or inclined by less that 45o to the normal
representing that direction. 

@section application Application

Identify to which model boundary a face of a tetrahedron in an input mesh 
belongs. In CSP, for instance, this comparison will be carried out in
MeshInterface::ReadTSolid(). 

@section messages Messages 

If the determined facing direction does not match any of the principal
axis, a message will report this. 

*/
int  Triangle3D::FacingDirection() const
  {
    static const double cut_off = 45.0;
    // normals corresponding to the possible directions
    static Edge3D  front ( Point3D(0.0,0.0,0.0), Point3D(0.0,0.0,1.0) );
    static Edge3D  back  ( Point3D(0.0,0.0,0.0), Point3D(0.0,0.0,-1.0) );
    static Edge3D  top   ( Point3D(0.0,0.0,0.0), Point3D(0.0,1.0,0.0) );
    static Edge3D  bottom( Point3D(0.0,0.0,0.0), Point3D(0.0,-1.0,0.0) );
    static Edge3D  left  ( Point3D(0.0,0.0,0.0), Point3D(-1.0,0.0,0.0) );
    static Edge3D  right ( Point3D(0.0,0.0,0.0), Point3D(1.0,0.0,0.0) );
    static Edge3D  tnormal;
    
    UnitNormal( tnormal );
 
    if ( bottom.AngleTo(tnormal) <= cut_off ) return MJL3D_BOTTOM;
    if ( front.AngleTo(tnormal)  <= cut_off ) return MJL3D_FRONT;
    if ( right.AngleTo(tnormal)  <= cut_off ) return MJL3D_RIGHT;
    if ( back.AngleTo(tnormal)   <= cut_off ) return MJL3D_BACK;
    if ( left.AngleTo(tnormal)   <= cut_off ) return MJL3D_LEFT;
    if ( top.AngleTo(tnormal)    <= cut_off ) return MJL3D_TOP;

    cout <<"\nTriangle3D::FacingDirection: ";
    cout <<"facing direction of triangle could not be identified..."<< endl;
    return MJL3D_TILTED;
  }
  
  
  
void  Triangle3D::Out() const
 {
    cout <<"\nTriangle3D::Out:"<< endl;
    cout <<"\nVertices: ";
    v0_.Out();
    v1_.Out();
    v2_.Out();
    cout <<"\nEdge marking diagonal of bounding box:";
    boundingBox_.Out();
    cout <<"\nCentral point (n): ";
    n_.Out();
    cout << endl;
 }
  
  
}  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
