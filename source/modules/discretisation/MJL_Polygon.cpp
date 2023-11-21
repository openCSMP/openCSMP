#include "MJL_Polygon.h"

using namespace std;

namespace mjl {


// points in list should be unique and in clockwise order
Polygon::Polygon( list<mjl::Point>::const_iterator first, 
                  list<mjl::Point>::const_iterator last )
 {
    v_ = new Vertex( *first++ );
 
    while ( first != last ) {
         v_ = v_->Insert( new Vertex( *first++ ) );
         size_++;
      }
    v_ = v_->Cw();
 }  



/// advances the argument polygon until its back to
/// where it was initially
Polygon&  Polygon::operator=( const Polygon& p )
 {
    if ( &p != this ) {
         size_ = p.size_;
         if ( size_ == 0 ) v_ = 0;
         else {
              v_ = new Vertex( p.Point() );
              for ( unsigned int i=1U; i<size_; i++ ) {
                   p.Advance( CLOCKWISE );
                   v_ = v_->Insert( new Vertex( p.Point() ) );
                }
           }
         p.Advance( CLOCKWISE );
         v_ = v_->Cw();
      }
    return *this; 
 }



Polygon::~Polygon()
 {
    if ( v_ ) {
         Vertex* w = v_->Cw();
         
         while( v_ != w ) {
              delete w->Remove();
              w = v_->Cw();
           }
         delete v_;
      }
 }


void Polygon::Erase()
 {
    if ( v_ ) {
         Vertex* w = v_->Cw();
         
         while( v_ != w ) {
              delete w->Remove();
              w = v_->Cw();
           }
         delete v_;
      }
      
    size_=0;
    v_=0;
 }
 
 

void Polygon::Resize()
 { 
    if ( v_ == 0 ) size_ = 0;
    else {
         Vertex* v = v_->Cw();
         for ( size_=1U; v!=v_; ++size_ )
           v = v->Cw();
      }
 }




void Polygon::BoundingRectangle( mjl::Point& cnr_min, mjl::Point& cnr_max ) const
 {
    cnr_min = cnr_max = Point();
 
    for ( unsigned int i=1u; i<=size_; i++, Advance(CLOCKWISE) ) {
         mjl::Point a = v_->Point();
         if ( a < cnr_min ) cnr_min = a;
         if ( a > cnr_max ) cnr_max = a;
      }
 }





void Polygon::Out() const
 {
    cout <<"\nPolygon::Out: "<< size_ <<" vertices (in clockwise order): ";
    
    for ( unsigned int i=1u; i<=size_; i++, Advance(CLOCKWISE) ) {
         mjl::Point a = v_->Point();
         cout <<"\n\t"<< a.X() <<"  "<< a.Y();
      }
    cout << endl;
 }



/// writes to file with extension .txt
/// point x,y,z are komma separated
void Polygon::Out( const char* file_txt ) const
 {
    string    file_name(file_txt); file_name +=".txt";
    ofstream  ofs( file_name.c_str() );
    
    if ( !ofs.is_open() ) {
         cout <<"\nPolygon::Out: can't open output file: "<< file_name << endl;
         return;
      }

    ofs <<"Polygon::Out: "<< size_ <<" vertices (in clockwise order): ";
    ofs <<"\n"<< size_;
    
    for ( unsigned int i=1u; i<=size_; i++, Advance(CLOCKWISE) ) {
         mjl::Point a = v_->Point();
         ofs <<"\n"<< a.X() <<","<< a.Y() <<",0.";
      }
    ofs << endl;
    
    cout <<"\nPolygon::Out: '"<< file_name <<"' written successfully."<< endl;
 }





void Polygon::OutputCoordinatesTo( std::list<mjl::Point>& points ) const
 {
    if ( !points.empty() ) points.erase( points.begin(), points.end() );

    for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) )
      points.push_back( v_->Point() );
 }



void Polygon::Revert()
 {
    list<mjl::Point>  points;
    
    for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) )
      points.push_front( v_->Point() );
      
    Erase();
    
    list<mjl::Point>::const_iterator it=points.begin();

    v_ = new Vertex( (*it++) );
 
    while ( it != points.end() )
      v_ = v_->Insert( new Vertex( (*it++) ) );
    v_ = v_->Cw();
    
    size_ = points.size();
 }



Vertex*  Polygon::LeastVertex( int (*cmp)( const mjl::Point& a, const mjl::Point& b ) ) // p.87
 {
    Vertex* bestV = V();
    
    Advance(CLOCKWISE);
    for ( int i=1; i<static_cast<int>(size_); i++, Advance(CLOCKWISE) )
      if ( (*cmp)( V()->Point(), bestV->Point() ) > 0 ) bestV = V();
    
    SetV(bestV);
    
    return bestV;
    
 } // end LeastVertex




void  Polygon::CenterOfGravity( double& x, double& y ) const
 {
     x = y = 0.;
 
     for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) ) {
          x += v_->Point()[0];
          y += v_->Point()[1];
       }
     x /= static_cast<double>(size_);
     y /= static_cast<double>(size_);

 } // end CenterOfGravity



double Polygon::Perimeter() const
 {
     double perimeter(0.);

     for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) ) {
          mjl::Edge segment = Edge();
          perimeter += segment.Length();
       }
       
    return perimeter;
 }



void  Polygon::Scale( double factor )
 {
     for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) )
       static_cast<mjl::Point>(*v_) *= factor;
 }
 
 
void  Polygon::Move( double dx, double dy )
 {
     for ( unsigned int i=0U; i<size_; i++, Advance(CLOCKWISE) ) {
          static_cast<mjl::Point>(*v_)(0) += dx;
          static_cast<mjl::Point>(*v_)(1) += dy;
       }
 }




/**
 
Loops around polygons and tests whether convex or concave angles between
Edges along the polygon boundary are smaller than 'tolerated_bound'.
If so, method returns false and prints the angle.  
 */
bool Polygon::CheckAngles( double min_angle_permitted ) const
 {
    bool  check(true);
    
    // always go one step back so that all vertices are tested
    for ( unsigned int i=0; i<size_; Advance(CLOCKWISE), i++ )
      {
         mjl::Point a = Ccw()->Point();
         mjl::Point p = Point(); // on perimeter of polygon between a and b
         mjl::Point b = Cw()->Point();
         // edge is oriented clockwise (see MJL definition p. 120, fig. 5.12)
         mjl::Edge e(a,b);
         double angle = signedAngle(p,e); 
         
         // testing the angle
         if ( fabs(angle) < min_angle_permitted ) {
              cout <<"\nMJL_Polygon::CheckAngles: angle is too small: "<< angle;
              cout <<" vs. "<< min_angle_permitted <<" at point:";
              p.Out();
              check = false;
           }
      }
      
    return check;
 
 } // end CheckAngles



mjl::Point  Polygon::InsidePointSA() // uses signed angle > and returns edge midpoint 
 {
    Vertex* curr = V();

    // always go one step back so that all vertices are tested
    for ( unsigned int i=0; i<size_; Advance(CLOCKWISE), i++ )
      {
         mjl::Point a = Ccw()->Point();
         mjl::Point p = Point(); // on perimeter of polygon between a and b
         mjl::Point b = Cw()->Point();
         // edge is oriented clockwise (see MJL definition p. 120, fig. 5.12)
         mjl::Edge e(a,b);
         double angle = signedAngle(p,e); 
         
         // if point is convex the midpoint of the respective edge is returned
         if ( angle > 0. and angle < 180. ) {
              SetV(curr); // restore original current vertex
              return e.MidPoint();
           }
      }
      
    cout <<"\nPolygon::InsidePointSA: unable to find point inside polygon."<< endl;  
    return mjl::Point(std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN());  
      
 } // end InsidePointSA 


} // end namespace csmp
