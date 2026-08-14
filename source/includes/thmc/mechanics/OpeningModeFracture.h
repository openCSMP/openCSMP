#ifndef OPENING_MODE_FRACTURE_H
#define OPENING_MODE_FRACTURE_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

class OpeningModeFracture {
  public:
    OpeningModeFracture();
    OpeningModeFracture( const OpeningModeFracture& );
    OpeningModeFracture( double length, double Poissons_ratio, double Youngs_modulus );
    ~OpeningModeFracture() {};
  
    void      Length( double length );
    double  Length() const;

    /// convention syy, the fracture normal stress, is positive if compressive
    double  CenterAperture( double pf, double syy ) const;
    /// calculates fracture aperture at some distance cx (0..length/2) from the fracture center
    double  Aperture( double pf, double syy, double cx ) const;
    /// returns (permeability/viscosity) * aperture product 
    double  Transmissivity( double pf, double syy, double visc ) const;
  
    /// joint volume (assuming unit thickness)
    double  Mode_I_Volume( double pf, double syy ) const;
  
    /// volume of a joint in shear (assuming unit thickness)
    double  Mode_II_Volume( double sxy ) const;
  
    /// maximum volume in response to tensile and shear stresses
    double  MaximumVolume( double pf, double syy, double sxy ) const;
  
    /// dilatation of a square block of rock with fracture intensity (f/m), and dimensions = frac length^2
    double  Dilatation( double frac_intensity, double pf, double syy, double sxy ) const;
  
    /// 2D construction routines for fracture polygons
    void    SixteenPointConvexHull( double pf,
                                    double syy,
                                    const Point<2>& a, 
                                    const Point<2>& b, std::list<Point<2>>& chain );
                                    
    void    BluntTenPointHull( double pf,
                               double syy,
                               const Point<2>& b,
                               const Point<2>& c, 
                               std::list<Point<2>>& chain );
                                    
    void    RectangleHull( double pf,
                           double syy,
                           Point<2>& b, Point<2>& c,
                           std::list<Point<2>>& chain );                                

    void    RectangleHull( double fixed_aperture,
                           Point<2>& b, Point<2>& c,
                           std::list<Point<2>>& chain );
  private:
    double        nu,        // Poisson's ratio
                  E,         // Young's modulus (GPa)
                  a,         // half-length of joint
                  x,         // distance from fracture center
                  pos,
                  length;    // fracture length
    const double  MINIMUM_APERTURE;
 };


} // namespace csmp


#endif
