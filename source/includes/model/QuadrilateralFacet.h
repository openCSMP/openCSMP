#ifndef QUADRILATERAL_FACET_H
#define QUADRILATERAL_FACET_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

/// SKM suggested workhorse: returns detJ into argument
Point<3U> normalAtFacetCenter( const Point<3U>&, 
                               const Point<3U>&, 
                               const Point<3U>&, 
                               const Point<3U>& );

Point<2U> normalAtFacetCenter( const Point<2U>&,
                               const Point<2U>&,
                               const Point<2U>&,
                               const Point<2U>& );

Point<1U> normalAtFacetCenter( const Point<1U>&,
                               const Point<1U>&,
                               const Point<1U>&,
                               const Point<1U>& );

/// 4-point quadrature
double64 facetArea4( const Point<3U>&, 
                     const Point<3U>&, 
                     const Point<3U>&, 
                     const Point<3U>& );

double64 facetArea4( const Point<2U>&,
                     const Point<2U>&,
                     const Point<2U>&,
                     const Point<2U>& );

double64 facetArea4( const Point<1U>&,
                     const Point<1U>&,
                     const Point<1U>&,
                     const Point<1U>& );

/// 1-point quadrature
double64 facetArea1( const Point<3U>&, 
                     const Point<3U>&, 
                     const Point<3U>&, 
                     const Point<3U>& );

/// analytical calculation of the area of a quadrilateral
double64 facetArea1( const Point<2U>&,
                     const Point<2U>&,
                     const Point<2U>&,
                     const Point<2U>& );


/**

@brief To compute area and unit-normal computations for linear quadrilaterals in 3D space.

@author Stephan Matthai
@date 2005

@note use triangular facet methods to do same for corresponding elements, but 
always prefer the native methods of the FiniteElement subclasses over this.

*/
class QuadrilateralFacet {
  public:
    QuadrilateralFacet( const Point<3U>&, 
                        const Point<3U>&, 
                        const Point<3U>&, 
                        const Point<3U>& );

    ~QuadrilateralFacet();

    /// returns detJ into its argument
    Point<3U>  NormalAtCenter() const;
    double64   Area() const;
  
  private:
    QuadrilateralFacet();
    
    Point<3U> pt0_, pt1_, pt2_, pt3_;
};

} // end namespace csmp

#endif
