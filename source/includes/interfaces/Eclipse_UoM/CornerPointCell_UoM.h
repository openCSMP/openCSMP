#ifndef CORNER_POINT_CELL_UOM_H
#define CORNER_POINT_CELL_UOM_H

#include "PolygonGrid_UoM.h"
#include "PolygonCell_UoM.h"
#include "Pillar_UoM.h"
#include "MeshingTools_UoM.h"

namespace csmp {

namespace eclipse {

enum CORNER_POINT_CELL_TYPE
{
    CORNER_POINT_CELL_UNDEFINED          = 0,
    CORNER_POINT_CELL_8_HEXAHEDRON       = 1,    /// observed and tested
    CORNER_POINT_CELL_7_AUGMENTED_PRISM  = 2,    /// observed and tested
    CORNER_POINT_CELL_6_WEDGE            = 3,    /// observed and tested
    CORNER_POINT_CELL_6_OCTAHEDRON       = 4,    /// observed and tested
    CORNER_POINT_CELL_5_AUGMENTED_TETRA  = 5,    /// observed and tested
    CORNER_POINT_CELL_5_PYRAMID          = 6,    /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
    CORNER_POINT_CELL_4_TETRAHEDRON      = 7,    /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
    CORNER_POINT_CELL_4_QUADRILATERAL    = 8,    /// observed and tested
    CORNER_POINT_CELL_3_TRIANGLE         = 9,    /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
    CORNER_POINT_CELL_2_BAR              = 10,   /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
    CORNER_POINT_CELL_1_POINT            = 11    /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
};

std::string toString( CORNER_POINT_CELL_TYPE );



enum CORNER_POINT_CELL_CATEGORY
{
    CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING = 0,
    CORNER_POINT_CELL_DIM_3_OVERLAPPING     = 1,
    CORNER_POINT_CELL_DIM_2                 = 2,
    CORNER_POINT_CELL_DIM_1                 = 3,
    CORNER_POINT_CELL_DIM_0                 = 4
};

std::string toString( CORNER_POINT_CELL_CATEGORY );



/// the order is important since it affects the remeshing functionality and inclusion of wells
///
///   Ymin  Zmin
///       \ |
/// Xmin __\|______ Xplus
///         |\
///         | \ Yplus
///         |
///         Zplus
///
///   0 ______________ 1
///    |\             |\
///    | \            | \
///    |  \           |  \
///    |   \2_________|___\3
///    |    |         |    |
///    |4___|________5|    |
///     \   |          \   |
///      \  |           \  |
///       \ |            \ |
///        \|_____________\|
///         6              7
///

enum CORNER_POINT_CELL_FACE_INDEX {
    CORNER_POINT_CELL_FACE_Xplus  = 4,
    CORNER_POINT_CELL_FACE_Xminus = 2,
    CORNER_POINT_CELL_FACE_Yplus  = 3,
    CORNER_POINT_CELL_FACE_Yminus = 1,
    CORNER_POINT_CELL_FACE_Zplus  = 0,
    CORNER_POINT_CELL_FACE_Zminus = 5
};

std::string toString( CORNER_POINT_CELL_FACE_INDEX );


/**

@brief CornerPointCell  provides the necessary functionality
to construct valid Finite Elements out of distorted hexahedrons.

@attention the CornerPoint cell is inherited from the PolygonCell !

@author R. Manasipov
@date 2015

*/
class CornerPointCell: public PolygonCell
{
public:
    explicit CornerPointCell( PolygonGridManager& );
    CornerPointCell( const CornerPointCell& poly );
    CornerPointCell& operator=( const CornerPointCell& poly );
    ~CornerPointCell();

    /// initialization of cell ( assign nodes, define topology )
    void   AssignNode( size_t local_nid, const csmp::Point<3U>& pt );
    void   InitializeCornerPointCell( );

    /// define cell which contains well
    void   AddCornerPointWellCell( size_t face_id_org, size_t face_id_dst );

    /// resolve mesh conflicts for cells
    void   ResolveCornerPointCellConflicts( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges );

    /// process cells
    void   ProcessCornerPointCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
  
    void   ProcessCornerPointRegularCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
  
    void   ProcessCornerPointDegenerateCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
  
    void   ProcessCornerPointDegenerateVolumetricCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
  
    void   ProcessCornerPointDegenerateVolumetricNonOverlappingCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
  
    void   ProcessCornerPointDegenerateVolumetricOverlappingCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );
    
    void   ProcessCornerPointDegenerateLowDimensionalCell( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element = false, bool tetra_mesh = false );

    /// returns the number of faces that make up the corner-point cell (only in degenerate cells can this be less than 6!)
    size_t GetNumFaces( size_t fid ) const;

    /// type of element ( dimension )
    bool   IsVolumeElement() const;
    bool   IsSurfaceElement() const;
    bool   IsLineElement() const;
    bool   IsPointElement() const;

    /// type of element ( topology )
    bool   IsVolumetricNonOverlapingElement() const;
    bool   IsVolumetricOverlapingElement() const;
    bool   IsLowDimensionalElement() const;

    /// wells info
    size_t GetNumWellSections() const;
    size_t GetNumWellSectionNodes( size_t eid ) const;
    const csmp::Point<3U>& GetWellSectionNode( size_t eid, size_t nid ) const;
    size_t GetWellSectionNodeGlobalId( size_t eid, size_t nid ) const;

    /// pillar nodes
    size_t GetPillarNodeGlobalIdOriginalOrder( size_t nid ) const;
    size_t GetPillarNodeGlobalIdCustomOrder( size_t nid ) const;
    const csmp::Point<3U>& GetPillarPointOriginalOrder( size_t nid ) const;
    const csmp::Point<3U>& GetPillarPointCustomOrder( size_t nid ) const;
    GridNode* GetPillarNodeOriginalOrder( size_t nid );
    GridNode* GetPillarNodeCustomOrder( size_t nid );
  
    /// prints object state to screen
    void Out() const;

protected:

    /// add new well
    void AddWell( std::map<std::set<size_t>,GridNode*>& additional_points );
    void AddWell( GridNode* n0, GridNode* n1 );

    /// initialize elements ( renumber node and eliminate duplicates )

    void AssignPillarNodesOrder( );

    void InitializeNodes( std::map<csmp::Point<3U>,size_t>& cell_nodes );
    void CheckNodeOrder();

    void InitializeFaces();
    void CheckFaceOrder();

    void InitializePoly8Element( );

    void InitializePoly7Element( const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                                 const std::vector<size_t>&                             tri_faces,
                                 const std::vector<size_t>&                             quad_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_tri_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_quad_faces,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_edges );

    void InitializePoly6Element( const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                                 const std::vector<size_t>&                             tri_faces,
                                 const std::vector<size_t>&                             quad_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_tri_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_quad_faces,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_edges );

    void InitializePoly5Element( const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                                 const std::vector<size_t>&                             tri_faces,
                                 const std::vector<size_t>&                             quad_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_tri_faces,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_quad_faces,
                                 const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                                 const std::map<csmp::Point<3U>,std::vector<size_t> >& point_edges );

    void InitializePoly4Element( const std::vector<size_t>& tri_faces,
                                 const std::vector<size_t>& quad_faces );

    void InitializePoly3Element( );

    void InitializePoly2Element( );

    void InitializePoly1Element( );

    void GetTopologycalInfo(    std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                                std::set<std::set<csmp::Point<3U> > >&          unique_faces,
                                std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                                std::vector<std::set<csmp::Point<3U> > >&       face_unique_points,
                                std::vector<size_t>&                             tri_faces,
                                std::vector<size_t>&                             quad_faces,
                                std::map<csmp::Point<3U>,std::vector<size_t> >& point_tri_faces,
                                std::map<csmp::Point<3U>,std::vector<size_t> >& point_quad_faces,
                                std::set<std::set<csmp::Point<3U> > >&          unique_edges,
                                std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                                std::vector<std::set<csmp::Point<3U> > >&       edge_unique_points,
                                std::map<csmp::Point<3U>,std::vector<size_t> >& point_edges );

    /// process cells
    void AddQuadFaceCentroids( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges );
    void ProcessPolyElement(  std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element );
    void ProcessPolyNonWellElement( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool add_edge_centroids, bool add_face_centroids, bool add_cell_centroids );
    void ProcessPolyWellElement( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges );
    void ProcessPoly8Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly7Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly6Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly5Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly4Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly3Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly2Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void ProcessPoly1Element( std::map<std::set<size_t>,GridNode*>& additional_points, std::set<std::set<size_t> >& additional_edges, bool well_element, bool tetra_mesh );
    void FinishMeshing();

private:

    /// cell specs
    CORNER_POINT_CELL_TYPE           cell_type_;
    CORNER_POINT_CELL_CATEGORY       cell_category_;
    size_t                           num_quad_faces_;
    size_t                           num_tri_faces_;
    int                              meshing_cycle_;
    /// pillar nodes data
    std::vector<GridNode*>     pillar_nodes_;
    std::vector<size_t>              pillar_nodes_order_;
    /// wells
    long                             well_face_org_;
    long                             well_face_dst_;
    std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >    wells_;
};

} // eclipse

}// end namespace csmp

#endif

