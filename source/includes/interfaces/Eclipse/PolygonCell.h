// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POLYGON_CELL_UOM_H
#define POLYGON_CELL_UOM_H

#include "Point.h"
#include "PolygonGrid.h"

namespace csmp {

/**

@class PolygonCell  PolygonCell "PolygonCell.h"
@author R. Manasipov
@date 2015

PolygonCell stores the cell topology.

*/
class PolygonCell
{
public:
    // TODO: dependency is the wrong way around!
    PolygonCell( PolygonGridManager& pgm );
    PolygonCell( const PolygonCell& cell );
    PolygonCell& operator=( const PolygonCell& cell );
    ~PolygonCell();

    /// cell dim
    size_t GetElementDim( ) const;

    /// elements info
    size_t GetNumElements() const;
    size_t GetNumElementNodes( size_t eid ) const;
    const csmp::Point<3U>& GetElementPoint( size_t eid, size_t nid ) const;
    GridNode* GetElementNode( size_t eid, size_t nid );
    size_t GetElementNodeGlobalId( size_t eid, size_t nid ) const;
    const csmp::CSMP_FEM_TYPE& GetElementType( size_t eid ) const;
    size_t GetElementDim( size_t eid ) const;

    /// faces info
    size_t GetNumFaces() const;
    /// returns into how many faces a particular face is split?
    size_t GetNumSubFaces( size_t fid ) const;
    /// should this be GetNumSubFaceNodes() where subfaces only exist when a quadrilateral faces is split into triangles?
    size_t GetNumFaceNodes( size_t fid, size_t sfid ) const;
    const csmp::Point<3U>& GetFacePoint( size_t fid, size_t sfid, size_t nid ) const;
    GridNode* GetFaceNode( size_t fid, size_t sfid, size_t nid );
    size_t GetFaceNodeGlobalId( size_t fid, size_t sfid, size_t nid ) const;
    const csmp::CSMP_FEM_TYPE& GetFaceType( size_t fid, size_t sfid ) const;
    size_t GetFaceDim( size_t fid, size_t sfid ) const;

    /// nodes info
    size_t GetNumNodes() const;
    size_t GetNodeGlobalIdOriginalOrder( size_t nid ) const;
    size_t GetNodeGlobalIdCustomOrder( size_t nid ) const;
    const csmp::Point<3U>& GetPointOriginalOrder( size_t nid ) const;
    const csmp::Point<3U>& GetPointCustomOrder( size_t nid ) const;
    GridNode* GetNodeOriginalOrder( size_t nid );
    GridNode* GetNodeCustomOrder( size_t nid );

    /// extra nodes info
    size_t GetNumExtraNodes() const;
    const csmp::Point<3U>& GetExtraPoint( size_t nid ) const;
    GridNode* GetExtraNode( size_t nid );
    void AddExtraNode( const csmp::Point<3U>& pt );

    /// cell centroid ( if exist: by convention it's a first node of extra nodes arrays )
    size_t GetCellCentroidGlobalId() const;
    const csmp::Point<3U>& GetCellCentroidPoint() const;
    GridNode* GetCellCentroidNode();

    /// face and node id's
    size_t GetOriginalNodeId( size_t custom_nid ) const;
    size_t GetCustomNodeId( size_t custom_nid ) const;
    size_t GetOriginalFaceId( size_t custom_fid ) const;
    size_t GetCustomFaceId( size_t custom_fid ) const;

    /// polygon faces info
    size_t GetNumPolygonFaceNodes( size_t fid ) const;
    size_t GetPolygonFaceNodeLocalId( size_t fid, size_t nid ) const;
    size_t GetPolygonFaceNodeGlobalId( size_t fid, size_t nid ) const;
    const csmp::Point<3U>& GetPolygonFacePoint( size_t fid, size_t nid ) const;
    GridNode* const GetPolygonFaceNode( size_t fid, size_t nid ) const;
  
    /// print object state to screen
    void Out() const;

protected:

    PolygonCell();

    /// polygon faces info
    void AssignPolygonFaceNode( size_t fid, size_t lnid, size_t cell_nid );
    void AddPolygonFaceNode( size_t fid, size_t cell_nid );
    void ClearPolygonFaceNodes( size_t fid );
    void ResizePolygonFaceNodes( size_t fid, size_t size );

    /// nodes
    void InitializeNodeOrder( size_t num_nodes );
    void AssignNodeOrder( size_t cnid, size_t onid );

    /// elements
    void AddCell( const std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> >& cell );
    void EraseCell( size_t position );
    void ClearCells();

    /// faces
    void InitializeFaceOrder( size_t num_faces );
    void AssignFaceOrder( size_t cfid, size_t ofid );
    void AddFace( size_t fid, const std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> >& face );
    void EraseFace( size_t fid );
    void EraseFace( size_t fid, size_t position );
    void ClearFaces( size_t fid );
    void ClearFaces();

protected:

    /// grid manager TODO: this class is only used by AddExtraNode(), not clear why a pointer to it is needed
    PolygonGridManager* grid_;

    /// cell nodes data
    size_t  num_nodes_;
    std::vector<GridNode*>   extra_nodes_;
    std::vector<GridNode*>   nodes_;
    std::vector<GridNode*>   nodes_in_custom_order_;
    std::vector<size_t>      custom_node_order_;

    /// cells TODO: what is the physical meaning of this data structure?
    std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > elements_;

    /// faces
    size_t                             num_faces_;
    std::vector<std::vector<size_t> >  face_nodes_;
    std::vector<std::vector<size_t> >  face_nodes_in_custom_order_;
    std::vector<size_t>                custom_face_order_;

    /// faces[i][j] is a 6 x 5 matrix of faces defined by face-type and corresponding node pointers as entries (type-value pairs)
    std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > > faces_;
};

}// end namespace csmp

#endif

