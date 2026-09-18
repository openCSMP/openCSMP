// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FRACMAN_FRACTURE_H
#define FRACMAN_FRACTURE_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

class FRACMAN_Fracture {
  public:
    FRACMAN_Fracture();
    FRACMAN_Fracture( const FRACMAN_Fracture& ffr );
    
    FRACMAN_Fracture( double        ap,
                   double        compr,
                   double        perm,
                   const std::list<double>& properties,
                   const std::list<Point<3>>& boundary );
                    
    ~FRACMAN_Fracture();
    FRACMAN_Fracture& operator=( const FRACMAN_Fracture& ffr );
    
    bool    InitializeFrom( long nprops, std::ifstream& ifs );

    /// compare perimeter length
    bool operator<( const FRACMAN_Fracture& ffr ) const;
    
    uint32_t    ID() const;
    std::string  TextID() const;
    int     SetID() const;
    std::string  TextSetID() const;
    void    BaryCenter( Point<3>& ctr ) const;
    void    BaryCenter( double& x, double& y, double& z ) const;
    double  Perimeter() const;
    double  Diameter() const; // from perimeter assuming circle
    void    BoundingBox( Point<3>& cnr1, Point<3>& cnr8 ) const;
    void    Move( double dx, double dy, double dz );
    void    Scale( double xfac, double dfac, double zfac );
    
    std::list<Point<3>>::const_iterator  Begin() const;
    std::list<Point<3>>::const_iterator  End() const;
    size_t                               PolygonPoints() const;
    std::list<double>::const_iterator    PropertiesBegin() const;
    std::list<double>::const_iterator    PropertiesEnd() const;
    
    void    Erase();
    
    void    Out() const;
  
  private:
    double  aperture;        // always present default properties
    double  compressibility;
    double  permeability;
    uint32_t    id;
    int     fracture_set_id;
    static  uint32_t global_id;
    
    std::list<double>    props;
    std::list<Point<3>>  boundary;
    Point<3>             unit_normal;
};


inline  std::list<Point<3>>::const_iterator  FRACMAN_Fracture::Begin() const
 {
    return boundary.begin();
 }
 
 
inline  std::list<Point<3>>::const_iterator  FRACMAN_Fracture::End() const
 {
    return boundary.end();
 }


inline  size_t  FRACMAN_Fracture::PolygonPoints() const
 {
    return boundary.size();
 }


inline  std::list<double>::const_iterator  FRACMAN_Fracture::PropertiesBegin() const
 {
    return props.begin();
 }
 
 
inline  std::list<double>::const_iterator  FRACMAN_Fracture::PropertiesEnd() const
 {
    return props.end();
 }

} // namespace csmp

#endif

