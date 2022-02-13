#ifndef CSMP_RHINO_SURFACE_READER_H
#define CSMP_RHINO_SURFACE_READER_H

#include "CSMP_definitions.h"
#include "MJL_geometry.h"
#include "MJL_Point3D.h"
#include "MJL_Triangle3D.h"

namespace csmp {

template<uint32_t> class VSet;

/**
    Reads triangulated render meshes of surfaces output from Rhinoceros (McNeel&Associates) and
    converts them into a 3D VSet that can be read by CSMP.
    
    Multiple surfaces are distinguished using the object names from Rhino.
    
     @author S. K. Matthai
     @date 2002
*/
class SKM_RhinoSurfaceReader {
   public:
     SKM_RhinoSurfaceReader();
     SKM_RhinoSurfaceReader( const char* raw_file );
     SKM_RhinoSurfaceReader( const SKM_RhinoSurfaceReader& );
     ~SKM_RhinoSurfaceReader();
     SKM_RhinoSurfaceReader& operator=( const SKM_RhinoSurfaceReader& );
     
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
                     std::map<size_t,mjl::Point3D>& points,
                     std::map<size_t,std::vector<int64_t> >& plist,
                     size_t poffset=0 ) const; ///< node/element numbering (0..n-1)
                     
     void ObjectToPData( const std::string& obj_name,
                         std::map<size_t,mjl::Point3D>& pxyz,
                         std::map<size_t,std::vector<int64_t> >& plist,
                         size_t poffset ) const; ///< node/element numbering (0..n-1)
                         
     void CreateNeighborPData( const std::map<size_t,std::vector<int64_t> >& plist,
                               std::map<size_t,std::vector<int64_t> >& pfverts,
                               std::vector<std::int8_t>& pbflags ) const;
                         
     void WriteGocadHeader( const char* surf_name, const char* GEOLOGICAL_TYPE, std::ofstream& ofs ) const;
     
     ///                                                          name      geological type
     size_t  SelectAndDescribeObjects( std::list<std::pair<std::string,std::string> >&  selections, 
                                       bool erase_list_before ) const;
   private:
     ///         obj.name   contained triangles
     std::map<std::string,std::list<mjl::Triangle3D> >  objects;
 };

} // csmp

#endif
