// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "OpeningModeFracture.h"
#include "CSMP_mathUtilities.h"

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



OpeningModeFracture::OpeningModeFracture( double l, double n, double E )
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
void  OpeningModeFracture::Length( double len )
 {
    length = len;
    a      = len / 2.;
 }


double  OpeningModeFracture::Length() const
 {
    return length;
 }



/// normal stress is negative if tensile
double  OpeningModeFracture::CenterAperture( double pf, double syy ) const
 {
    assert( pf > 0.);
    if ( syy > pf ) {
         std::cerr <<"\nOpeningModeFracture::CenterAperture: This fracture would be closed ";
         std::cerr <<"since syy > pf. If length < 5 m an aperture of 1 mm is assigned (else 1 cm)..." << std::endl;
         return MINIMUM_APERTURE;
      }
    return (2. * ((-syy)+pf)*(1.0-nu)*std::sqrt( (a*a) )) / E;
 }


double  OpeningModeFracture::Aperture( double pf, double syy, double cx ) const
 {
    assert( pf > 0.);
    if ( syy > pf ) {
         std::cerr <<"\nOpeningModeFracture::Aperture: This fracture would be closed ";
         std::cerr <<"since syy > pf. If length < 5 m aperture of 1 mm is assigned (else 1 cm)..." << std::endl;
         return MINIMUM_APERTURE;
      }
    return (2. * ((-syy)+pf)*(1.0-nu)*std::sqrt( (a*a - cx*cx) )) / E;
 }


double  OpeningModeFracture::Transmissivity( double pf, double syy, double visc ) const
  {
     assert( pf > 0.);
     // Renshaw 95, JGR 100:B2, 24,629-24,636, Tf in equation 13, p. 24,631
     double ap  = CenterAperture( pf, syy );
     ap *= ap;
     ap *= ap;
     return ap / (12. * visc);
  }




/// joint volume (assuming unit thickness)
double  OpeningModeFracture::Mode_I_Volume( double p, double syy ) const
 {
    assert( p > 0.);
    if ( syy - p > 0. ) return MINIMUM_APERTURE * length;
    return (1./4. * (1 - nu) * (syy + p) * (length*length) * PI) / E;
   
 } // end Mode_I_Volume


/// volume of a joint in shear (assuming unit thickness)
double  OpeningModeFracture::Mode_II_Volume( double sxy ) const
 {
    assert( sxy >= 0.);
    if ( sxy < numeric_limits<double>::epsilon() ) return MINIMUM_APERTURE * length;
    return (1./4. * sxy * (1 - nu) * (length*length) * PI) / E;

 } // end Mode_II_Volume


/// maximum volume in response to tensile and shear stresses
double  OpeningModeFracture::MaximumVolume( double pf, double syy, double sxy ) const
 {
    assert( pf > 0.);
    assert( sxy >= 0.);
    if ( sxy < numeric_limits<double>::epsilon() && syy > pf ) return MINIMUM_APERTURE * length;
    return std::max( Mode_I_Volume( pf, syy ), Mode_II_Volume( sxy ) );
   
 } // end MaximumVolume


/**
   dilatation of a square block of rock with fracture intensity (f/m), and dimensions = frac length^2
*/
double  OpeningModeFracture::Dilatation( double frac_intensity, double pf, double syy, double sxy ) const
 {
    assert( pf > 0.);
    assert( sxy >= 0.);
    assert( frac_intensity >= 0. );
    if ( sxy < numeric_limits<double>::epsilon() && syy > pf ) return 0.;
    //      n fractures
    return (frac_intensity * length) * MaximumVolume( pf, syy, sxy ) / (length * length);
   
 } // end Dilatation








void OpeningModeFracture::SixteenPointConvexHull( double pf,
                                                  double syy,
                                                  const Point<2>& b, // left tip
                                                  const Point<2>& c, // right tip
                                                  list<Point<2>>& chain )
{
    // 1. Basic Geometry
    Point<2> vec = c - b;
    double L = distance(b, c);
    double a_half = L / 2.0;

    // Unit vector along fracture and the perpendicular (normal) vector
    Point<2> u = vec / L;
    Point<2> n(-u[1], u[0]); // 90 deg CCW rotation

    // 2. Sample apertures at specific locations along the fracture
    // Locations: near tip (1%), slightly further (10%), quarter point (50% of half-length)
    double a0 = std::max(Aperture(pf, syy, 0.0), MINIMUM_APERTURE);            // Center
    double a1 = std::max(Aperture(pf, syy, a_half * 0.9), MINIMUM_APERTURE);   // 10% from tip
    double a2 = std::max(Aperture(pf, syy, a_half * 0.99), MINIMUM_APERTURE);  // 1% from tip
    double a3 = std::max(Aperture(pf, syy, a_half * 0.5), MINIMUM_APERTURE);   // Mid-wing

    chain.clear();

    // 3. Define the sampling positions (t from 0 to 1) and their associated apertures
    // We follow the CCW order: Start at tip B, go along top, reach C, return along bottom.
    struct Sample { double t; double half_ap; };
    vector<Sample> top_samples = {
        { (a_half * 0.01) / L, a2 / 2.0 },
        { (a_half * 0.1) / L,  a1 / 2.0 },
        { (a_half * 0.5) / L,  a3 / 2.0 },
        { a_half / L,          a0 / 2.0 }, // Center
        { (a_half + a_half * 0.5) / L, a3 / 2.0 },
        { (L - a_half * 0.1) / L,      a1 / 2.0 },
        { (L - a_half * 0.01) / L,     a2 / 2.0 }
    };

    // 4. Build the Chain
    chain.push_back(b); // Start at Left Tip

    // Add Top Points (Offset by +n)
    for (const auto& s : top_samples) {
        chain.push_back(b + (u * (s.t * L)) + (n * s.half_ap));
    }

    chain.push_back(c); // Right Tip

    // Add Bottom Points (Offset by -n, in reverse order of t)
    for (auto it = top_samples.rbegin(); it != top_samples.rend(); ++it) {
        chain.push_back(b + (u * (it->t * L)) - (n * it->half_ap));
    }
} 





void OpeningModeFracture::BluntTenPointHull( double pf,
                                             double syy,
                                             const Point<2>& b, // left tip
                                             const Point<2>& c, // right tip
                                             list<Point<2>>& chain )
{
    // 1. Basic Geometry
    // Handle the b > c flip logic using simple swap if needed
    Point<2> start = b;
    Point<2> end = c;
    if (start[0] > end[0] || (start[0] == end[0] && start[1] > end[1])) {
        std::swap(start, end);
    }

    Point<2> vec = end - start;
    double L = distance(start, end);
    double a_half = L / 2.0;

    // Unit vector along fracture and the perpendicular normal (90 deg CCW)
    Point<2> u = vec / L;
    Point<2> n(-u[1], u[0]); 

    // 2. Calculate Apertures
    double a0 = std::max(Aperture(pf, syy, 0.0), MINIMUM_APERTURE);            // Center
    double a1 = std::max(Aperture(pf, syy, a_half * 0.9), MINIMUM_APERTURE);   // Near tips
    double a3 = std::max(Aperture(pf, syy, a_half * 0.5), MINIMUM_APERTURE);   // Mid-wing

    chain.clear();

    // 3. Define sampling positions along the segment (t = 0.0 to 1.0)
    // Note: p1 and p5 are at the actual tips (0.0 and 1.0)
    struct Sample { double t; double half_ap; };
    vector<Sample> samples = {
        { 0.0, a1 / 2.0 },    // Left blunt tip
        { 0.25, a3 / 2.0 },   // Mid-left
        { 0.5, a0 / 2.0 },    // Center
        { 0.75, a3 / 2.0 },   // Mid-right
        { 1.0, a1 / 2.0 }     // Right blunt tip
    };

    // 4. Build the Chain (Counter-Clockwise)
    // First pass: Top side (Positive normal offset)
    for (const auto& s : samples) {
        chain.push_back(start + (u * (s.t * L)) + (n * s.half_ap));
    }

    // Second pass: Bottom side (Negative normal offset, reverse order)
    for (auto it = samples.rbegin(); it != samples.rend(); ++it) {
        chain.push_back(start + (u * (it->t * L)) - (n * it->half_ap));
    }
}
 








void OpeningModeFracture::RectangleHull( double pf, 
                                         double syy,
                                         Point<2>& b,
                                         Point<2>& c,
                                         list<Point<2>>& chain )
 {
    // 1. Basic Geometry & Sorting
    Point<2> start = b;
    Point<2> end = c;
    
    // Ensure consistent CCW winding by sorting points
    if (start[0] > end[0] || (start[0] == end[0] && start[1] > end[1])) {
        std::swap(start, end);
    }

    Point<2> vec = end - start;
    double L = distance(start, end);

    // Unit vector along fracture (u) and perpendicular normal (n)
    Point<2> u = vec / L;
    Point<2> n(-u[1], u[0]); // 90 deg CCW rotation

    // 2. Aperture handling
    double aperture = CenterAperture(pf, syy);
    if (aperture < MINIMUM_APERTURE) aperture = MINIMUM_APERTURE;
    double half_ap = aperture / 2.0;

    // 3. Clear and rebuild the chain
    chain.clear();

    // The four corners of the rectangle:
    // p1: Bottom-Left (offset back from start and down)
    // p2: Top-Left    (offset back from start and up)
    // p3: Top-Right   (offset forward from end and up)
    // p4: Bottom-Right(offset forward from end and down)
    
    // Using your original logic where the rectangle slightly overshoots the tips
    // by (half_ap) to encapsulate the tips.
    
    chain.push_back(start - (u * half_ap) - (n * half_ap)); // p1
    chain.push_back(start - (u * half_ap) + (n * half_ap)); // p2
    chain.push_back(end   + (u * half_ap) + (n * half_ap)); // p3
    chain.push_back(end   + (u * half_ap) - (n * half_ap)); // p4
    
 } // end RectangleHull






void OpeningModeFracture::RectangleHull( double fixed_aperture,
                                         Point<2>& b,
                                         Point<2>& c,
                                         list<Point<2>>& chain )
{
    // 1. Basic Geometry and Sorting
    Point<2> start = b;
    Point<2> end = c;
    
    // Maintain consistent winding order regardless of input point order
    if (start[0] > end[0] || (start[0] == end[0] && start[1] > end[1])) {
        std::swap(start, end);
    }

    Point<2> vec = end - start;
    double L = distance(start, end);

    // Safeguard for zero-length fractures
    if (L < 1e-12) return;

    // Unit vector along fracture (u) and perpendicular normal (n)
    Point<2> u = vec / L;
    Point<2> n(-u[1], u[0]); // 90 deg CCW rotation

    double half_ap = fixed_aperture / 2.0;

    // 2. Clear existing chain
    chain.clear();

    // 3. Calculate 4 corners
    // The original logic extends the rectangle beyond the tips by half_ap
    // to create a box that fully "contains" the segment nodes.
    
    // p1: Back-Bottom
    chain.push_back(start - (u * half_ap) - (n * half_ap)); 
    // p2: Back-Top
    chain.push_back(start - (u * half_ap) + (n * half_ap)); 
    // p3: Front-Top
    chain.push_back(end   + (u * half_ap) + (n * half_ap)); 
    // p4: Front-Bottom
    chain.push_back(end   + (u * half_ap) - (n * half_ap)); 
    
} // end RectangleHull






// TEST PROGRAMS
// =============

/*
    // making a rectangular fracture
    OpeningModeFracture frac10;
    Point<2>           x(2.0,5.0), y(8.0,4.0);
    list<Point<2>>     chain_x;

    frac10.RectangleHull( 1.0e+9, 2.0e+5, x, y, chain_x );
    mjl::Polygon poly10( chain_x.begin(), chain_x.end() );
    poly10.Out("rectangle");
*/


/*
   // BUILDING TWO FRACTURES AND INTERSECTING THEM
   // --------------------------------------------

   OpeningModeFracture frac1, frac2;
   Point<2>           a(2.0,1.0), b(8.0,1.0);
   list<Point<2>>     chain;
    
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
   
   Point<2> p = poly3.Point();
   
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


