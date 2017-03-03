#ifndef FRACMAN_FRACTURE_H
#define FRACMAN_FRACTURE_H

#include "MJL_Point3D.h"
#include "MJL_Edge3D.h"

#include "CSMP_definitions.h"

namespace csmp {

class FRACMAN_Fracture {
  public:
    FRACMAN_Fracture();
    FRACMAN_Fracture( const FRACMAN_Fracture& ffr );
    
    FRACMAN_Fracture( double64        ap,
                   double64        compr,
                   double64        perm,
                   const std::list<double64>& properties,
                   const std::list<mjl::Point3D>& boundary );
                    
    ~FRACMAN_Fracture();
    FRACMAN_Fracture& operator=( const FRACMAN_Fracture& ffr );
    
    bool    InitializeFrom( int nprops, std::ifstream& ifs );

    /// compare perimeter length
    bool operator<( const FRACMAN_Fracture& ffr ) const;
    
    uint32    ID() const;
    std::string  TextID() const;
    int     SetID() const;
    std::string  TextSetID() const;
    void    BaryCenter( mjl::Point3D& ctr ) const;
    void    BaryCenter( double64& x, double64& y, double64& z ) const;
    double64  Perimeter() const;
    double64  Diameter() const; // from perimeter assuming circle
    void    BoundingBox( mjl::Point3D& cnr1, mjl::Point3D& cnr8 ) const;
    void    Move( double64 dx, double64 dy, double64 dz );
    void    Scale( double64 xfac, double64 dfac, double64 zfac );
    
    std::list<mjl::Point3D>::const_iterator  Begin() const;
    std::list<mjl::Point3D>::const_iterator  End() const;
    size_t                                  PolygonPoints() const;
    std::list<double64>::const_iterator    PropertiesBegin() const;
    std::list<double64>::const_iterator    PropertiesEnd() const;
    
    void    Erase();
    
    void    Out() const;
  
  private:
    double64  aperture;        // always present default properties
    double64  compressibility;
    double64  permeability;
    uint32    id;
    int     fracture_set_id;
    static  uint32 global_id;
    
    std::list<double64>     props;
    std::list<mjl::Point3D>  boundary;
    mjl::Edge3D              unit_normal;
};


inline  std::list<mjl::Point3D>::const_iterator  FRACMAN_Fracture::Begin() const
 {
    return boundary.begin();
 }
 
 
inline  std::list<mjl::Point3D>::const_iterator  FRACMAN_Fracture::End() const
 {
    return boundary.end();
 }


inline  size_t  FRACMAN_Fracture::PolygonPoints() const
 {
    return boundary.size();
 }


inline  std::list<double64>::const_iterator  FRACMAN_Fracture::PropertiesBegin() const
 {
    return props.begin();
 }
 
 
inline  std::list<double64>::const_iterator  FRACMAN_Fracture::PropertiesEnd() const
 {
    return props.end();
 }

} // namespace csmp

#endif

