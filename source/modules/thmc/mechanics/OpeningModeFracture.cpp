#include "OpeningModeFracture.h"
#include "CSMP_mathUtilities.h"
#include "MJL_Edge.h"
#include "MJL_Polygon.h"

using namespace std;

namespace csmp {

OpeningModeFracture::OpeningModeFracture()          // 0.1 mm
   : length(15.5), a(15.5/2.), nu(0.23), E(4.5e10), MINIMUM_APERTURE(1.0e-4)
 {
    cout <<"\nOpeningModeFracture::default: "; 
    cout <<"mean laboratory data for granite, D.D. Pollard 2/27/95, chapter 10A, p. 46.";
    cout <<"\n    Poisson's ratio:   0.23 ";
    cout <<"\n    Young's modulus:  45.0 GPa, "<< endl;
 }



OpeningModeFracture::OpeningModeFracture( double64 l, double64 n, double64 E )
   : length(l), a(l/2.), nu(n), E(E), MINIMUM_APERTURE(1.0e-4)
 {
 }


/// copy-constructor
OpeningModeFracture::OpeningModeFracture( const OpeningModeFracture& of )
 : nu(of.nu), E(of.E),
   a(of.a), x(of.x), pos(of.pos), length(of.length),
   MINIMUM_APERTURE(of.MINIMUM_APERTURE)
 {
 }


/**
    NB: also setting fracture half-length a !
*/
void  OpeningModeFracture::Length( double64 len )
 {
    length = len;
    a      = len / 2.;
 }


double64  OpeningModeFracture::Length() const
 {
    return length;
 }



/// normal stress is negative if tensile
double64  OpeningModeFracture::CenterAperture( double64 pf, double64 syy ) const
 {
    assert( pf > 0.);
    if ( syy > pf ) {
         std::cerr <<"\nOpeningModeFracture::CenterAperture: This fracture would be closed ";
         std::cerr <<"since syy > pf. If length < 5 m an aperture of 1 mm is assigned (else 1 cm)..." << std::endl;
         return MINIMUM_APERTURE;
      }
    return (2. * ((-syy)+pf)*(1.0-nu)*std::sqrt( (a*a) )) / E;
 }


double64  OpeningModeFracture::Aperture( double64 pf, double64 syy, double64 cx ) const
 {
    assert( pf > 0.);
    if ( syy > pf ) {
         std::cerr <<"\nOpeningModeFracture::Aperture: This fracture would be closed ";
         std::cerr <<"since syy > pf. If length < 5 m aperture of 1 mm is assigned (else 1 cm)..." << std::endl;
         return MINIMUM_APERTURE;
      }
    return (2. * ((-syy)+pf)*(1.0-nu)*std::sqrt( (a*a - cx*cx) )) / E;
 }


double64  OpeningModeFracture::Transmissivity( double64 pf, double64 syy, double64 visc ) const
  {
     assert( pf > 0.);
     // Renshaw 95, JGR 100:B2, 24,629-24,636, Tf in equation 13, p. 24,631
     double64 ap  = CenterAperture( pf, syy );
     ap *= ap;
     ap *= ap;
     return ap / (12. * visc);
  }




/// joint volume (assuming unit thickness)
double64  OpeningModeFracture::Mode_I_Volume( double64 p, double64 syy ) const
 {
    assert( p > 0.);
    if ( syy - p > 0. ) return MINIMUM_APERTURE * length;
    return (1./4. * (1 - nu) * (syy + p) * (length*length) * PI) / E;
   
 } // end Mode_I_Volume


/// volume of a joint in shear (assuming unit thickness)
double64  OpeningModeFracture::Mode_II_Volume( double64 sxy ) const
 {
    assert( sxy >= 0.);
    if ( sxy < numeric_limits<double64>::epsilon() ) return MINIMUM_APERTURE * length;
    return (1./4. * sxy * (1 - nu) * (length*length) * PI) / E;

 } // end Mode_II_Volume


/// maximum volume in response to tensile and shear stresses
double64  OpeningModeFracture::MaximumVolume( double64 pf, double64 syy, double64 sxy ) const
 {
    assert( pf > 0.);
    assert( sxy >= 0.);
    if ( sxy < numeric_limits<double64>::epsilon() && syy > pf ) return MINIMUM_APERTURE * length;
    return std::max( Mode_I_Volume( pf, syy ), Mode_II_Volume( sxy ) );
   
 } // end MaximumVolume


/**
   dilatation of a square block of rock with fracture intensity (f/m), and dimensions = frac length^2
*/
double64  OpeningModeFracture::Dilatation( double64 frac_intensity, double64 pf, double64 syy, double64 sxy ) const
 {
    assert( pf > 0.);
    assert( sxy >= 0.);
    assert( frac_intensity >= 0. );
    if ( sxy < numeric_limits<double64>::epsilon() && syy > pf ) return 0.;
    //      n fractures
    return (frac_intensity * length) * MaximumVolume( pf, syy, sxy ) / (length * length);
   
 } // end Dilatation








// tested: O.K.
void OpeningModeFracture::SixteenPointConvexHull( double64 pf,
                                                  double64 syy,
                                                  const mjl::Point& b, // left
                                                  const mjl::Point& c, // right
                                                  list<mjl::Point>& chain )
 {
    mjl::Point point;
    mjl::Edge  line( b, c );
    
    // getting the length of the fracture given by 'a' and 'b'
    // O.K.
    length = sqrt( (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) );
    a      = length / 2.0;
           
    // getting apertures at x0, x1 and x2 
    // O.K.      
    double64 a0 = Aperture( pf, syy, 0.0 );
    double64 a1 = Aperture( pf, syy, a - a*1.0e-1 );
    double64 a2 = Aperture( pf, syy, a - a*1.0e-2 );
    double64 a3 = Aperture( pf, syy, a / 2.0 );
    
    // checking apertures
    if ( a0 < MINIMUM_APERTURE ) a0 = MINIMUM_APERTURE;
    if ( a1 < MINIMUM_APERTURE ) a1 = MINIMUM_APERTURE;
    if ( a2 < MINIMUM_APERTURE ) a2 = MINIMUM_APERTURE;
    if ( a3 < MINIMUM_APERTURE ) a3 = MINIMUM_APERTURE;
 
    // checking supplied chain of points
    if ( !chain.empty() )
      chain.erase( chain.begin(), chain.end() );

    // Along the fracture find the two points p and q that
    // are separated by a distance corresponding to the 
    // fracture aperture and that have, as a center between
    // them, the point along the fracture from which this
    // aperture is to be constructed.
    
    // 1. find the 2 points for the first fracture opening
    double64 t_p = (a * 1.0e-2) / length - (a2/2.0) / length, 
           t_q = (a * 1.0e-2) / length + (a2/2.0) / length;    
         
   // 2. build new edge for those points  
   mjl::Edge  oe1( line.Point(t_q), line.Point(t_p) );
   mjl::Point p2, p3, p4, p5, p6, p7, p8, p10, p11, p12,
             p13, p14, p15, p16;
     
   // rotate edge 90o counter clockwise around its center
   oe1.Rot();
   p2  = oe1.Origin();
   p16 = oe1.Destination();
    
   // second point on segment 
   t_p = (a * 1.0e-1) / length  -  (a1/2.0) / length; 
   t_q = (a * 1.0e-1) / length  +  (a1/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p3  = oe1.Origin();
   p15 = oe1.Destination();

   // third point on segment 
   t_p = (a / 2.0) / length  -  (a3/2.0) / length; 
   t_q = (a / 2.0) / length  +  (a3/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p4  = oe1.Origin();
   p14 = oe1.Destination();

   // fracture center 
   t_p = a / length  -  (a0/2.0) / length; 
   t_q = a / length  +  (a0/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p5  = oe1.Origin();
   p13 = oe1.Destination();

   // fourth point on segment 
   t_p = (a + a / 2.0) / length  -  (a3/2.0) / length; 
   t_q = (a + a / 2.0) / length  +  (a3/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p6  = oe1.Origin();
   p12 = oe1.Destination();

   // third point on segment 
   t_p = (length - a * 1.0e-1) / length  -  (a1/2.0) / length; 
   t_q = (length - a * 1.0e-1) / length  +  (a1/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p7  = oe1.Origin();
   p11 = oe1.Destination();

   // fourth point on segment 
   t_p = (length - a * 1.0e-2) / length  -  (a2/2.0) / length; 
   t_q = (length - a * 1.0e-2) / length  +  (a2/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p8 = oe1.Origin();
   p10 = oe1.Destination();

   // adding points to chain in counter-clockwise order
   chain.push_back( b );
   chain.push_back( p2 );
   chain.push_back( p3 );
   chain.push_back( p4 );
   chain.push_back( p5 );
   chain.push_back( p6 );
   chain.push_back( p7 );
   chain.push_back( p8 );
   chain.push_back( c );
   chain.push_back( p10 );
   chain.push_back( p11 );
   chain.push_back( p12 );
   chain.push_back( p13 );
   chain.push_back( p14 );
   chain.push_back( p15 );
   chain.push_back( p16 );
   
 } // end SixteenPointConvexHull
 





// tested: O.K. 
void OpeningModeFracture::BluntTenPointHull( double64 pf,
                                             double64 syy,
                                             const mjl::Point& b, // left
                                             const mjl::Point& c, // right
                                             list<mjl::Point>& chain )
 {
    mjl::Point point;
    mjl::Edge  line( b, c );
    
    // make sure that counter-clockwise building occurs
    if ( b > c ) line.Flip();
    
    // getting the length of the fracture given by 'a' and 'b'
    // O.K.
    length = sqrt( (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) );
    a      = length / 2.0;
           
    // getting apertures at x0, x1 and x2 
    // O.K.      
    double64 a0 = Aperture( pf, syy, 0.0 );
    double64 a1 = Aperture( pf, syy, a - a*0.1 ); // is moved to fracture tip
    double64 a3 = Aperture( pf, syy, a / 2.0 );

    // checking apertures
    if ( a0 < MINIMUM_APERTURE ) a0 = MINIMUM_APERTURE;
    if ( a1 < MINIMUM_APERTURE ) a1 = MINIMUM_APERTURE;
    if ( a3 < MINIMUM_APERTURE ) a3 = MINIMUM_APERTURE;
 
    // checking supplied chain of points
    if ( !chain.empty() )
      chain.erase( chain.begin(), chain.end() );

    // Along the fracture find the two points p and q that
    // are separated by a distance corresponding to the 
    // fracture aperture and that have, as a center between
    // them, the point along the fracture from which this
    // aperture is to be constructed.
    mjl::Point p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
    
    // 1. find the 2 points for the blunt end of fracture
    //    exaggerating its aperture
    double64 t_p = 0.0  -  (a1/2.0) / length, 
           t_q = 0.0  +  (a1/2.0) / length;    
         
   // 2. build new edge for those points  
   mjl::Edge  oe1( line.Point(t_q), line.Point(t_p) );
     
   // rotate edge 90o counter clockwise around its center
   oe1.Rot();
   p1  = oe1.Origin();
   p10 = oe1.Destination();
    
   // second segment along fracture
   t_p = (a / 2.0) / length  -  (a3/2.0) / length; 
   t_q = (a / 2.0) / length  +  (a3/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p2  = oe1.Origin();
   p9  = oe1.Destination();

   // fracture center 
   t_p = a / length  -  (a0/2.0) / length; 
   t_q = a / length  +  (a0/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p3  = oe1.Origin();
   p8  = oe1.Destination();

   // fourth point on segment 
   t_p = (a + a / 2.0) / length  -  (a3/2.0) / length; 
   t_q = (a + a / 2.0) / length  +  (a3/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p4  = oe1.Origin();
   p7  = oe1.Destination();

   // blunt second fracture tip 
   t_p = 1.0  -  (a1/2.0) / length; 
   t_q = 1.0  +  (a1/2.0) / length;    
   oe1.Set( line.Point(t_q), line.Point(t_p) );
   oe1.Rot();
   p5  = oe1.Origin();
   p6  = oe1.Destination();

   // adding points to chain in counter-clockwise order
   chain.push_back( p1 );
   chain.push_back( p2 );
   chain.push_back( p3 );
   chain.push_back( p4 );
   chain.push_back( p5 );
   chain.push_back( p6 );
   chain.push_back( p7 );
   chain.push_back( p8 );
   chain.push_back( p9 );
   chain.push_back( p10 );
   
 } // end BluntTenPointHull
 








void OpeningModeFracture::RectangleHull( double64 pf, double64 syy,
                                         mjl::Point& b,
                                         mjl::Point& c,
                                         list<mjl::Point>& chain )
 {
    mjl::Edge  line( b, c );

    // make sure that counter-clockwise building occurs
    if ( b > c ) line.Flip();

    double64 aperture = CenterAperture( pf, syy );

    // checking apertures
    if ( aperture < MINIMUM_APERTURE ) aperture = MINIMUM_APERTURE;

    // 1. find the 2 points for the end edge of the fracture
    length = sqrt( (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) );
    double64 t_p = -(aperture/2.0) / length, 
           t_q =  (aperture/2.0) / length;  
         
    // 2. build new edge for those points  
    mjl::Edge  oe1( line.Point(t_q), line.Point(t_p) );
     
    // rotate edge 90o counter clockwise around its center
    oe1.Rot();
    mjl::Point p1 = oe1.Origin();
    mjl::Point p2 = oe1.Destination();
    
    // second point on segment 
    t_p = 1.0 - (aperture/2.0) / length; 
    t_q = 1.0 + (aperture/2.0) / length;    
    oe1.Set( line.Point(t_q), line.Point(t_p) );
    oe1.Rot();
    mjl::Point p4 = oe1.Origin();
    mjl::Point p3 = oe1.Destination();
    
    // checking supplied chain of points
    if ( !chain.empty() )
      chain.erase( chain.begin(), chain.end() );

    // adding points to chain in counter-clockwise order
    chain.push_back( p1 );
    chain.push_back( p2 );
    chain.push_back( p3 );
    chain.push_back( p4 );
    
 } // end RectangleHull                                                                                







void OpeningModeFracture::RectangleHull( double64 fixed_aperture,
                                         mjl::Point& b,
                                         mjl::Point& c,
                                         list<mjl::Point>& chain )
 {
    mjl::Edge  line( b, c );

    // make sure that counter-clockwise building occurs
    if ( b > c ) line.Flip();

    double64 aperture = fixed_aperture;

    // 1. find the 2 points for the end edge of the fracture
    length = sqrt( (c[0]-b[0])*(c[0]-b[0]) + (c[1]-b[1])*(c[1]-b[1]) );
    double64 t_p = -(aperture/2.0) / length, 
           t_q =  (aperture/2.0) / length;  
         
    // 2. build new edge for those points  
    mjl::Edge  oe1( line.Point(t_q), line.Point(t_p) );
     
    // rotate edge 90o counter clockwise around its center
    oe1.Rot();
    mjl::Point p1 = oe1.Origin();
    mjl::Point p2 = oe1.Destination();
    
    // second point on segment 
    t_p = 1.0 - (aperture/2.0) / length; 
    t_q = 1.0 + (aperture/2.0) / length;    
    oe1.Set( line.Point(t_q), line.Point(t_p) );
    oe1.Rot();
    mjl::Point p4 = oe1.Origin();
    mjl::Point p3 = oe1.Destination();
    
    // checking supplied chain of points
    if ( !chain.empty() )
      chain.erase( chain.begin(), chain.end() );

    // adding points to chain in counter-clockwise order
    chain.push_back( p1 );
    chain.push_back( p2 );
    chain.push_back( p3 );
    chain.push_back( p4 );
    
 } // end RectangleHull                                                                                







// TEST PROGRAMS
// =============

/*
    // making a rectangular fracture
    OpeningModeFracture frac10;
    mjl::Point           x(2.0,5.0), y(8.0,4.0);
    list<mjl::Point>     chain_x;

    frac10.RectangleHull( 1.0e+9, 2.0e+5, x, y, chain_x );
    mjl::Polygon poly10( chain_x.begin(), chain_x.end() );
    poly10.Out("rectangle");
*/


/*
   // BUILDING TWO FRACTURES AND INTERSECTING THEM
   // --------------------------------------------

   OpeningModeFracture frac1, frac2;
   mjl::Point           a(2.0,1.0), b(8.0,1.0);
   list<mjl::Point>     chain;
    
   frac1.SixteenPointConvexHull( 1.0e+9, 2.0e+5, a, b, chain );
   mjl::Polygon poly1( chain.begin(), chain.end() );
   poly1.Out();
 
   a.Set(6.0,1.1);
   b.Set(12.0,1.1);  
   frac2.SixteenPointConvexHull( 1.0e+9, 2.0e+5, a, b, chain );
   mjl::Polygon poly2( chain.begin(), chain.end() );
   poly2.Out(); 
   
   mjl::Polygon poly3;
   
   poly3.BuildIntersectionPolygon( poly1, poly2 );
   poly3.Out();
   
   mjl::Point p = poly3.Point();
   
   poly1.RemoveSharedPoints( poly3 );
   poly2.RemoveSharedPoints( poly3 );
  
   combinePolygonsTo( poly1, poly2, chain );
  
   mjl::Polygon poly4;
   
   poly4.BuildStarShapedAround( p, chain.begin(), chain.end() );
   
   // removing point in kernel
   poly4.Remove( p );
   poly4.Out();


*/


} // end namespace csmp


