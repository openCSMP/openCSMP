#include "QuadrilateralFacet.h"

using namespace std;

namespace csmp {

// Private default constructor, should never be used

QuadrilateralFacet::QuadrilateralFacet()
{
} // end constructor


QuadrilateralFacet::QuadrilateralFacet( const Point<3U>& p0,
                                        const Point<3U>& p1,
                                        const Point<3U>& p2,
                                        const Point<3U>& p3 )
: pt0_(p0),
  pt1_(p1),
  pt2_(p2),
  pt3_(p3)
{
} // end constructor



QuadrilateralFacet::~QuadrilateralFacet()
{
}


double64 QuadrilateralFacet::Area() const
 {
    return facetArea4( pt0_, pt1_, pt2_, pt3_ );
 }



Point<3U> QuadrilateralFacet::NormalAtCenter() const
 {
    return normalAtFacetCenter( pt0_, pt1_, pt2_, pt3_ );
 }

/**

This computes the global shape function derivative matrix at the
center of the quadrilateral facet from the Jacobian and the node
coordinates. With this result the normal
is constructed as the crossproduct of the derivative vectors 
and subsequently scaled to unit length. 
*/
Point<3U>  normalAtFacetCenter( const Point<3U>& pt0, 
                                const Point<3U>& pt1, 
                                const Point<3U>& pt2, 
                                const Point<3U>& pt3 )
{
  // x - coordinate
  double64 jac00(0.25 * (-pt0[0] + pt1[0] + pt2[0] - pt3[0])); // dnr
  double64 jac10(0.25 * (-pt0[0] - pt1[0] + pt2[0] + pt3[0])); // dns
  
  // y - coordinate
  double64 jac01(0.25 * (-pt0[1] + pt1[1] + pt2[1] - pt3[1]));
  double64 jac11(0.25 * (-pt0[1] - pt1[1] + pt2[1] + pt3[1]));

  // z - coordinate
  double64 jac02(0.25 * (-pt0[2] + pt1[2] + pt2[2] - pt3[2]));
  double64 jac12(0.25 * (-pt0[2] - pt1[2] + pt2[2] + pt3[2]));
   
  Point<3U> vecNormal( jac01*jac12-jac02*jac11,
                       jac02*jac10-jac00*jac12,
                       jac00*jac11-jac01*jac10 );

  vecNormal.NormalizeLengthTo(1.); 

  return vecNormal;
  
} // end normalAtFacetCenter

/// normal is zero as it points into the coordinate direction that does not exist
Point<2U>  normalAtFacetCenter( const Point<2U>& pt0,
                                const Point<2U>& pt1,
                                const Point<2U>& pt2,
                                const Point<2U>& pt3 )
{
    return Point<2U>();
}

/// normal is zero as it points into the coordinate direction that does not exist
Point<1U>  normalAtFacetCenter( const Point<1U>& pt0,
                                const Point<1U>& pt1,
                                const Point<1U>& pt2,
                                const Point<1U>& pt3 )
{
    return Point<1U>();
}

/// one-point quadrature
double64  facetArea1( const Point<3U>& pt0, 
                      const Point<3U>& pt1,
                      const Point<3U>& pt2,
                      const Point<3U>& pt3 )
{
  // x - coordinate
  double64 jac00(0.25 * (-pt0[0] + pt1[0] + pt2[0] - pt3[0])); // dnr
  double64 jac10(0.25 * (-pt0[0] - pt1[0] + pt2[0] + pt3[0])); // dns
  
  // y - coordinate
  double64 jac01(0.25 * (-pt0[1] + pt1[1] + pt2[1] - pt3[1]));
  double64 jac11(0.25 * (-pt0[1] - pt1[1] + pt2[1] + pt3[1]));

  // z - coordinate
  double64 jac02(0.25 * (-pt0[2] + pt1[2] + pt2[2] - pt3[2]));
  double64 jac12(0.25 * (-pt0[2] - pt1[2] + pt2[2] + pt3[2]));
   
  double64 efg0(jac00 * jac00 + jac01 * jac01 + jac02 * jac02);
  double64 efg1(jac00 * jac10 + jac01 * jac11 + jac02 * jac12);
  double64 efg2(jac10 * jac10 + jac11 * jac11 + jac12 * jac12);
  
  return 4. * sqrt(efg0 * efg2 - efg1 * efg1);
  
} // end facetArea1

double64  facetArea1( const Point<2U>& pt0,
                      const Point<2U>& pt1,
                      const Point<2U>& pt2,
                      const Point<2U>& pt3 )
{
    return 1.;
}

double64  facetArea1( const Point<1U>& pt0,
                      const Point<1U>& pt1,
                      const Point<1U>& pt2,
                      const Point<1U>& pt3 )
{
    return 1.;
}

/**

Computes the area of the facet using 4-point numerical integration
with symmetric integration point locations about sqrt(1/3).  

@section arguments Input Arguments 

The corner points of the quadrilateral facet.  

@return Returns the facet area.  

@section implementation Implementation

Local shape function derivatives dnr and dns at integration points: 

IP 1, dnr, dns:
-0.105662432702594  0.105662432702594 0.394337567297406 -0.394337567297406
-0.105662432702594 -0.394337567297406 0.394337567297406  0.105662432702594

IP 2, dnr, dns:
-0.105662432702594  0.105662432702594 0.394337567297406 -0.394337567297406
-0.394337567297406 -0.105662432702594 0.105662432702594  0.394337567297406

IP 3, dnr, dns:
-0.394337567297406  0.394337567297406 0.105662432702594 -0.105662432702594
-0.394337567297406 -0.105662432702594 0.105662432702594  0.394337567297406

IP 4, dnr, dns:
-0.394337567297406  0.394337567297406 0.105662432702594 -0.105662432702594
-0.105662432702594 -0.394337567297406 0.394337567297406  0.105662432702594


@section application Application

Finite volume computations on volumetric finite elements.  
*/
double64  facetArea4( const Point<3U>& p0, 
                      const Point<3U>& p1, 
                      const Point<3U>& p2, 
                      const Point<3U>& p3 )
{
  // integration point 1
  // -------------------
  // Jacobian column 1           dnr1 * x1 +               dnr2 * x2 
  double64 jac00(-0.105662432702594*p0[0] + 0.105662432702594*p1[0] + 0.394337567297406*p2[0] - 0.394337567297406*p3[0]); // dnr
  double64 jac10(-0.105662432702594*p0[0] - 0.394337567297406*p1[0] + 0.394337567297406*p2[0] + 0.105662432702594*p3[0]); // dns
  
  // column 2
  double64 jac01(-0.105662432702594*p0[1] + 0.105662432702594*p1[1] + 0.394337567297406*p2[1] - 0.394337567297406*p3[1]);
  double64 jac11(-0.105662432702594*p0[1] - 0.394337567297406*p1[1] + 0.394337567297406*p2[1] + 0.105662432702594*p3[1]);

  // column 3
  double64 jac02(-0.105662432702594*p0[2] + 0.105662432702594*p1[2] + 0.394337567297406*p2[2] - 0.394337567297406*p3[2]);
  double64 jac12(-0.105662432702594*p0[2] - 0.394337567297406*p1[2] + 0.394337567297406*p2[2] + 0.105662432702594*p3[2]);
   
  // compute J' 1 
  double64 efg0(jac00 * jac00 + jac01 * jac01 + jac02 * jac02);
  double64 efg1(jac00 * jac10 + jac01 * jac11 + jac02 * jac12);
  double64 efg2(jac10 * jac10 + jac11 * jac11 + jac12 * jac12);
  double64 detJ(sqrt(efg0 * efg2 - efg1 * efg1));


  // integration point 2
  // -------------------
  // x - coordinate
  jac00 = -0.105662432702594*p0[0] + 0.105662432702594*p1[0] + 0.394337567297406*p2[0] - 0.394337567297406*p3[0]; // dnr
  jac10 = -0.394337567297406*p0[0] - 0.105662432702594*p1[0] + 0.105662432702594*p2[0] + 0.394337567297406*p3[0]; // dns
  
  // y - coordinate
  jac01 = -0.105662432702594*p0[1] + 0.105662432702594*p1[1] + 0.394337567297406*p2[1] - 0.394337567297406*p3[1];
  jac11 = -0.394337567297406*p0[1] - 0.105662432702594*p1[1] + 0.105662432702594*p2[1] + 0.394337567297406*p3[1];

  // z - coordinate
  jac02 = -0.105662432702594*p0[2] + 0.105662432702594*p1[2] + 0.394337567297406*p2[2] - 0.394337567297406*p3[2];
  jac12 = -0.394337567297406*p0[2] - 0.105662432702594*p1[2] + 0.105662432702594*p2[2] + 0.394337567297406*p3[2];

  // compute J' 2 
  efg0 = jac00 * jac00 + jac01 * jac01 + jac02 * jac02;
  efg1 = jac00 * jac10 + jac01 * jac11 + jac02 * jac12;
  efg2 = jac10 * jac10 + jac11 * jac11 + jac12 * jac12;
  detJ += sqrt(efg0 * efg2 - efg1 * efg1);


  // integration point 3
  // -------------------
  // x - coordinate
  jac00 = -0.394337567297406*p0[0] + 0.394337567297406*p1[0] + 0.105662432702594*p2[0] - 0.105662432702594*p3[0]; // dnr
  jac10 = -0.394337567297406*p0[0] - 0.105662432702594*p1[0] + 0.105662432702594*p2[0] + 0.394337567297406*p3[0]; // dns
  
  // y - coordinate
  jac01 = -0.394337567297406*p0[1] + 0.394337567297406*p1[1] + 0.105662432702594*p2[1] - 0.105662432702594*p3[1];
  jac11 = -0.394337567297406*p0[1] - 0.105662432702594*p1[1] + 0.105662432702594*p2[1] + 0.394337567297406*p3[1];

  // z - coordinate
  jac02 = -0.394337567297406*p0[2] + 0.394337567297406*p1[2] + 0.105662432702594*p2[2] - 0.105662432702594*p3[2];
  jac12 = -0.394337567297406*p0[2] - 0.105662432702594*p1[2] + 0.105662432702594*p2[2] + 0.394337567297406*p3[2];

  // compute J' 3 
  efg0 = jac00 * jac00 + jac01 * jac01 + jac02 * jac02;
  efg1 = jac00 * jac10 + jac01 * jac11 + jac02 * jac12;
  efg2 = jac10 * jac10 + jac11 * jac11 + jac12 * jac12;
  detJ += sqrt(efg0 * efg2 - efg1 * efg1);


  // integration point 4
  // -------------------
  // x - coordinate
  jac00 = -0.394337567297406*p0[0] + 0.394337567297406*p1[0] + 0.105662432702594*p2[0] - 0.105662432702594*p3[0]; 
  jac10 = -0.105662432702594*p0[0] - 0.394337567297406*p1[0] + 0.394337567297406*p2[0] + 0.105662432702594*p3[0]; 
  
  // y - coordinate
  jac01 = -0.394337567297406*p0[1] + 0.394337567297406*p1[1] + 0.105662432702594*p2[1] - 0.105662432702594*p3[1];
  jac11 = -0.105662432702594*p0[1] - 0.394337567297406*p1[1] + 0.394337567297406*p2[1] + 0.105662432702594*p3[1];

  // z - coordinate
  jac02 = -0.394337567297406*p0[2] + 0.394337567297406*p1[2] + 0.105662432702594*p2[2] - 0.105662432702594*p3[2];
  jac12 = -0.105662432702594*p0[2] - 0.394337567297406*p1[2] + 0.394337567297406*p2[2] + 0.105662432702594*p3[2];

  // compute J' 4 
  efg0 = jac00 * jac00 + jac01 * jac01 + jac02 * jac02;
  efg1 = jac00 * jac10 + jac01 * jac11 + jac02 * jac12;
  efg2 = jac10 * jac10 + jac11 * jac11 + jac12 * jac12;
  detJ += sqrt(efg0 * efg2 - efg1 * efg1);

  if ( detJ <= 0. ) {
     std::cout <<"\nfacetArea: Erroneous determinant of 2-3D Jacobian matrix: ";
     std::cout << detJ << endl;
     cout.flush();
     throw std::range_error("facetArea");
  }
   
  return detJ;
  
} // end facetArea


double64  facetArea4( const Point<2U>& pt0,
                      const Point<2U>& pt1,
                      const Point<2U>& pt2,
                      const Point<2U>& pt3 )
{
    return 1.;
}

double64  facetArea4( const Point<1U>& pt0,
                      const Point<1U>& pt1,
                      const Point<1U>& pt2,
                      const Point<1U>& pt3 )
{
    return 1.;
}



} // end namespace csmp
