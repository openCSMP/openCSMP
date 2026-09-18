// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_FRACMAN_INTERFACE_H
#define CSMP_FRACMAN_INTERFACE_H

#include "FracMan_Fracture.h"
#include "Point.h"

namespace csmp {

class FRACMAN_Interface {
  public:
    FRACMAN_Interface();
    ~FRACMAN_Interface();
    
    void   InitializeFrom_FRACMAN_File( const char* ffb_file );
    
    /// returns volume of box
    double BoundingBox( Point<3>& pmin, Point<3>& pmax ) const;
    
    void   MoveGeometryToOrigin();
    
    void   ScaleGeometry( double xfac, double yfac, double zfac );
  
    void   OutputToDXF( const char* dxf_file ) const;
    
    void   OutputToTETIN( const char* tetin_file ) const;
    
    void   OutputSelectedRegions( const char* sregions ) const;

    void   OutputSelectedFractureDiameters( const char* sregions ) const;
    
    void   Out() const;
    
  private:
    /// returns number of points
    long ReadFORMAT( std::ifstream& ifs, long& scale, long& fracs,
                     long& props, bool ascii );
    
    /// returns whether property specs were read correctly
    bool ReadPROPERTIES( std::ifstream& ifs, long n_props, std::list<std::string>& props );
    
    // stored data 
    std::map<std::string,FRACMAN_Fracture> fractures;
    std::list<std::string>                 properties;
    long         n_points,
                 n_fractures,
                 n_properties,
                 scale_factor;
    std::string  dataset;
};

} // end namespace csmp

#endif
