#ifndef OPENING_MODE_FRACTURE_H
#define OPENING_MODE_FRACTURE_H

#include "CSMP_definitions.h"
#include "MJL_Point.h"

namespace csmp {

class OpeningModeFracture {
  public:
    OpeningModeFracture();
    OpeningModeFracture( const OpeningModeFracture& );
    OpeningModeFracture( double64 length, double64 Poissons_ratio, double64 Youngs_modulus );
    ~OpeningModeFracture() {};
  
    void      Length( double64 length );
    double64  Length() const;

    /// convention syy, the fracture normal stress, is positive if compressive
    double64  CenterAperture( double64 pf, double64 syy ) const;
    /// calculates fracture aperture at some distance cx (0..length/2) from the fracture center
    double64  Aperture( double64 pf, double64 syy, double64 cx ) const;
    /// returns (permeability/viscosity) * aperture product 
    double64  Transmissivity( double64 pf, double64 syy, double64 visc ) const;
  
    /// joint volume (assuming unit thickness)
    double64  Mode_I_Volume( double64 pf, double64 syy ) const;
  
    /// volume of a joint in shear (assuming unit thickness)
    double64  Mode_II_Volume( double64 sxy ) const;
  
    /// maximum volume in response to tensile and shear stresses
    double64  MaximumVolume( double64 pf, double64 syy, double64 sxy ) const;
  
    /// dilatation of a square block of rock with fracture intensity (f/m), and dimensions = frac length^2
    double64  Dilatation( double64 frac_intensity, double64 pf, double64 syy, double64 sxy ) const;
  
    /// 2D construction routines for fracture polygons
    void    SixteenPointConvexHull( double64 pf,
                                    double64 syy,
                                    const mjl::Point& a, 
                                    const mjl::Point& b, std::list<mjl::Point>& chain );
                                    
    void    BluntTenPointHull( double64 pf,
                               double64 syy,
                               const mjl::Point& b,
                               const mjl::Point& c, 
                               std::list<mjl::Point>& chain );
                                    
    void    RectangleHull( double64 pf,
                           double64 syy,
                           mjl::Point& b, mjl::Point& c,
                           std::list<mjl::Point>& chain );                                

    void    RectangleHull( double64 fixed_aperture,
                           mjl::Point& b, mjl::Point& c,
                           std::list<mjl::Point>& chain );
  private:
    double64        nu,        // Poisson's ratio
                    E,         // Young's modulus (GPa)
                    a,         // half-length of joint
                    x,         // distance from fracture center
                    pos,
                    length;    // fracture length
    const double64  MINIMUM_APERTURE;
 };


} // namespace csmp


#endif