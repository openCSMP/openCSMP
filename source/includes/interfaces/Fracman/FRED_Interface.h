#ifndef CSMP_FRED_INTERFACE_H
#define CSMP_FRED_INTERFACE_H

#include "FRED_Fracture.h"

namespace csmp {

/// FRED was a DFN-modelling consortium run by Golder Associates <2004
class FRED_Interface {
  public:
    FRED_Interface();
    ~FRED_Interface();
    
    void   InitializeFrom_FRED_File( const char* ffb_file );
    
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
    long ReadFORMAT( std::ifstream& ifs, double& scale, long& fracs,
                     int& props, bool ascii );
    
    /// returns whether property specs were read correctly
    bool ReadPROPERTIES( std::ifstream& ifs, int n_props, std::list<std::string>& props );
    
    // stored data 
    std::map<std::string,FRED_Fracture> fractures;
    std::list<std::string>              properties;
    long                                n_points,
                                        n_fractures;
    int                                 n_properties;
    double                              scale_factor;
    std::string  dataset;
};

} // end namespace csmp

#endif
