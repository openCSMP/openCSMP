#ifndef MJL_TRIANGLE_3D_H
#define MJL_TRIANGLE_3D_H

#include "MJL_Point3D.h"
#include "MJL_Edge3D.h"

namespace mjl {

#define min3(A,B,C) ((A)<(B) ? ((A)<(C)?(A):(C)) : ((B)<(C)?(B):(C)))

#define max3(A,B,C) ((A)>(B) ? ((A)>(C)?(A):(C)) : ((B)>(C)?(B):(C)))

#define MJL3D_BOTTOM -4
#define MJL3D_FRONT  -6
#define MJL3D_RIGHT  -3
#define MJL3D_BACK   -7
#define MJL3D_LEFT   -2
#define MJL3D_TOP    -5
#define MJL3D_TILTED -1

#define EPSILON1 1e-12

enum { POSITIVE, NEGATIVE, ON };

  
class Triangle3D {
  public:  
    Triangle3D()  {};
    Triangle3D( const Point3D& v0, const Point3D& v1, const Point3D& v2, long id );
    Triangle3D( const Triangle3D& t ) { *this = t; };
    ~Triangle3D() {};
    Triangle3D& operator=( const Triangle3D& t );
    Point3D     operator[]( int i ) const;
    Edge3D      boundingBox()       const { return boundingBox_; };
    Point3D     n()                 const { return n_; };
    void        UnitNormal( Edge3D& e ) const; 
    void        UnitNormal( std::vector<double>& vec ) const;
    int         Classify( const Point3D& p ) const;
    int         FacingDirection() const;
    void        Set( const Point3D& v0, const Point3D& v1, const Point3D& v2, long id );
    void        Out() const { Out(std::cout); }
    void        Out(std::ostream& os) const;

    long     id_;
    int      mark_;

  private:
    Point3D  v0_, v1_, v2_, n_;
    Edge3D   boundingBox_;
};



inline Point3D  Triangle3D::operator[]( int i ) const
 { 
    if ( i == 0 ) return v0_;
    if ( i == 1 ) return v1_;
    return v2_;
 }


inline int  Triangle3D::Classify( const Point3D& p ) const
   {
        Point3D v  = p - v0_;
        double     len = v.Length();
        if ( len == 0. ) return ON;
        v = (1. / len) * v;
        double d = v.DotProduct(n_);
        if ( d > EPSILON1)       return POSITIVE;
        else if ( d < -EPSILON1) return NEGATIVE;
        return ON;
   }
    

inline void  Triangle3D::UnitNormal( Edge3D& e ) const
  {
     e.org_  = Point3D(0.,0.,0.);
     e.dest_ = n_;
  }


inline void  Triangle3D::UnitNormal( std::vector<double>& vec ) const
  {
     vec.resize(3u);
     vec[0] = n_.X();
     vec[1] = n_.Y();
     vec[2] = n_.Z();
  }

}

#endif















