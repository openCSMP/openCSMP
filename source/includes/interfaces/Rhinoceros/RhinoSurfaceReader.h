// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_RHINO_SURFACE_READER_H
#define CSMP_RHINO_SURFACE_READER_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

class LightWeightTriangle {
public:
    // Member initializer list with perfect forwarding or moves 
    // for maximum efficiency
    LightWeightTriangle(const Point<3>& p1, const Point<3>& p2, const Point<3>& p3, int64_t id) noexcept
        : id_(id), nodes_{p1, p2, p3} {}

    // Default constructor for container pre-allocation
    LightWeightTriangle() noexcept : id_(-1) {}

    // Coordinate access - noexcept for high-performance loops
    // Using const and non-const overloads
    [[nodiscard]] const Point<3>& operator[](size_t i) const noexcept {  return nodes_[i];  }
    
    [[nodiscard]] Point<3>& operator[](size_t i) noexcept { return nodes_[i]; }

    // Modern SetVertex using move semantics
    void SetVertex(size_t i, Point<3> p) noexcept { if (i < 3) nodes_[i] = std::move(p); }

    // Efficient swap for winding reversal
    void ReverseWinding() noexcept { std::swap(nodes_[1], nodes_[2]); }

    // Bulk update method
    void Set(Point<3> p1, Point<3> p2, Point<3> p3, int64_t id) noexcept {
          nodes_[0] = std::move(p1);
          nodes_[1] = std::move(p2);
          nodes_[2] = std::move(p3);
          id_ = id;
      }

    int64_t id_;

private:
    Point<3> nodes_[3];
};


template<uint32_t> class VSet;

/**
    Reads triangulated render meshes of surfaces output from Rhinoceros (McNeel&Associates) and
    converts them into a 3D VSet that can be read by CSMP.
    
    Multiple surfaces are distinguished using the object names from Rhino.
    
     @author S. K. Matthai
     @date 2002
*/
class RhinoSurfaceReader {
   public:
     RhinoSurfaceReader();
     RhinoSurfaceReader( const char* raw_file );
     RhinoSurfaceReader( const RhinoSurfaceReader& );
     ~RhinoSurfaceReader();
     RhinoSurfaceReader& operator=( const RhinoSurfaceReader& );
     
     /// conversion of triangulated surface object 'obj' into polygonal data set stored in NCSA VSet format
     void OutputObjectTo( const char* obj, VSet<3U>& ) const;
     
     void EraseObjects();   
     
     /// GoCAD output: TSurf format: a single selected surface
     void WriteObjectToTSurf( const char* obj, std::ofstream& ) const;
     
     /// GoCAD output: TSurf format: all surfaces in the file
     void WriteObjectsToTSurf( const char* tsurf_file, const char* GEOLOGICAL_TYPE="fault" ) const;

     /// GoCAD output: TSurf format: surfaces selected by the user through stdin
     void WriteSelectedObjectsToTSurf( const char* tsurf_file ) const;

     void ExchangeCoordinateAxes( int axis_a, int axis_b );
     
     void ScaleCoordinates( double x_fac, double y_fac, double z_fac );
     
     void MoveCoordinates( double x_move, double y_move, double z_move );
     
     size_t Objects() const;
     
   private:
     bool InitializeFrom( const char* raw_file, bool erase_old=true );

     bool PopObject( const char *obj_name,
                     std::map<size_t,Point<3>>& points,
                     std::map<size_t,std::vector<size_t> >& plist,
                     size_t poffset=0 ) const; ///< node/element numbering (0..n-1)
                     
     void ObjectToPData( const std::string& obj_name,
                         std::map<size_t,Point<3>>& pxyz,
                         std::map<size_t,std::vector<size_t> >& plist,
                         size_t poffset ) const; ///< node/element numbering (0..n-1)
                         
     void CreateNeighborPData( const std::map<size_t,std::vector<size_t> >& plist,
                               std::map<size_t,std::vector<int64_t> >& pfverts,
                               std::vector<std::int8_t>& pbflags ) const;
                         
     void WriteGocadHeader( const char* surf_name, const char* GEOLOGICAL_TYPE, std::ofstream& ofs ) const;
     
     ///                                                          name      geological type
     size_t  SelectAndDescribeObjects( std::list<std::pair<std::string,std::string> >&  selections, 
                                       bool erase_list_before ) const;
   private:
     ///         obj.name   contained triangles
     std::map<std::string, std::vector<LightWeightTriangle>> objects;
};

} // csmp

#endif
