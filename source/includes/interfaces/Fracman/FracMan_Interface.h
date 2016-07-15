#ifndef CSMP_FRACMAN_INTERFACE_H
#define CSMP_FRACMAN_INTERFACE_H

#include "FracMan_Fracture.h"

namespace csmp {

class FRACMAN_Interface {
  public:
    FRACMAN_Interface();
    ~FRACMAN_Interface();
    
    void   InitializeFrom_FRACMAN_File( const char* ffb_file );
    
    /// returns volume of box
    double BoundingBox( mjl::Point3D& pmin, mjl::Point3D& pmax ) const;
    
    void   MoveGeometryToOrigin();
    
    void   ScaleGeometry( double xfac, double yfac, double zfac );
  
    void   OutputToDXF( const char* dxf_file ) const;
    
    void   OutputToTETIN( const char* tetin_file ) const;
    
    void   OutputSelectedRegions( const char* sregions ) const;

    void   OutputSelectedFractureDiameters( const char* sregions ) const;
    
    void   Out() const;
    
  private:
    /// returns number of points
    long ReadFORMAT( std::ifstream& ifs, int& scale, int& fracs, 
                     int& props, bool ascii );
    
    /// returns whether property specs were read correctly
    bool ReadPROPERTIES( std::ifstream& ifs, int n_props, std::list<std::string>& props );
    
    // stored data 
    std::map<std::string,FRACMAN_Fracture> fractures;
    std::list<std::string>                 properties;
    int     n_points,
            n_fractures,
            n_properties,
            scale_factor;
    std::string  dataset;
};

} // end namespace csmp

#endif
