#ifndef CSMP_VSET_CONVERTER_H
#define CSMP_VSET_CONVERTER_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class VSet;

/**
    To convert linear to quadratic elements and to remove degenerate elements 
    from the input mesh. Methods are only for triangles and tetrahedra.
    
    @author Stephan Matthai
    @date 2001
 
*/
template<uint32_t dim>
class VSetConverter {
  public:
    /// the coordinates of the bounding box are initialized to Not A Number
    VSetConverter() : xmin(std::numeric_limits<double>::quiet_NaN()), xmax(std::numeric_limits<double>::quiet_NaN()),
                      ymin(std::numeric_limits<double>::quiet_NaN()), ymax(std::numeric_limits<double>::quiet_NaN()),
                      zmin(std::numeric_limits<double>::quiet_NaN()), zmax(std::numeric_limits<double>::quiet_NaN()) {};

    ~VSetConverter() {};
  
    /// applies the naming conventions TOP, BOTTOM etc. to sides, edges and cornier points of the box
    void EstablishBoundaryFlagsForBoxModel( VSet<dim>&, double tolerance=0.2 );
    
    /// turn triangular element mesh data to target element data
    void ConvertLinearToQuadraticTriangles( VSet<dim>& );
    
    /// turn triangular element mesh data to target 3D surface element data
    void ConvertLinearToQuadraticTriangles3D( VSet<dim>& );

    /// create a mesh of 7-noded barycentric triangular finite elements
    void ConvertLinearToBarycentricTriangles( VSet<dim>& );
    
    /// create 10-noded tetrahedra from 4-noded ones
    void ConvertLinearToQuadraticTetrahedra( VSet<dim>& );
    
    /// turn linear tetrahedral mesh into 11-noded barycentric tetrahedral mesh
    void ConvertLinearToBarycentricTetrahedra( VSet<dim>& );
    
  private:
    double xmin, xmax, ymin, ymax, zmin, zmax; ///< coordinates values of bounding box of model
  
    // 2D triangular meshes
    // --------------------
  
    /// rectangular model in 2 and 3D, flagging of corner node
    void FlagCornerNodes( VSet<dim>&, bool three_dimensional=false ) const;
    
    /// returns 0 if face is not at the model boundary
    int8_t  TestForBoundaryFlags( const std::vector<std::int8_t>& bflags,
                                  size_t nID1, size_t nID2 ) const;
  
    /// interpolation between value pairs
    double BoundaryValue( const std::map<size_t,double>& bvals, 
                            size_t nID1, size_t nID2 ) const;
  
    // from linear to quadratic mesh
    //void InterpolateNodeProperties( size_t nodes, int32_t enodes, VSet<dim>& ) const;
    
    void OrderQuadraticTriangleCoordinateOrigins( VSet<dim>& ) const;

    void OrderBarycentricQuadraticTriangleCoordinateOrigins( VSet<dim>& ) const;
    
    // 3D tetrahedral meshes
    // ---------------------                    default as in most Rhino models
    void   FlagEdges( VSet<dim>& vset, double tolerance=1.0e-4 ) const;

    /// flag-based check
    int8_t  BoundaryFlags3D( const std::vector<std::int8_t>& bflags,
                             size_t nID1, size_t nID2 ) const;
  
    /// coordinate-based check of whether a node is in on the boundary of a box shaped model
    int8_t  TestForBoundaryFlags3D( double x, double y, double z ) const;
 };


// inline function definitions

/**
    Interpolation of variable values on the boundary, assuming that it lies in one 
    of the coordinate planes.
*/
template<uint32_t dim>
inline double VSetConverter<dim>::BoundaryValue( const std::map<size_t,double>& bvals, 
                                                   size_t nID1, size_t nID2 ) 
 const
  {
      typename std::map<size_t,double>::const_iterator  bvit1(bvals.find(nID1)), 
                                                          bvit2(bvals.find(nID2));
      assert( bvit1 != bvals.end() );
      assert( bvit2 != bvals.end() );

      return ((*bvit1).second + (*bvit2).second) / 2.;

  } // end BoundaryValue
                                         



} // csmp

#endif
