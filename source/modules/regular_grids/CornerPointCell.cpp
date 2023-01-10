#include "CornerPointCell.h"

#include "ErrorHandler.h"

namespace csmp {

// CORNER POINT CELL

CornerPointCell::CornerPointCell( PolygonGridManager& pgm )
: PolygonCell( pgm ),
 meshing_cycle_(0),
 pillar_nodes_(8U),
 pillar_nodes_order_(8U,0),
 well_face_org_(-1),
 well_face_dst_(-1)
{
}




CornerPointCell::CornerPointCell( const CornerPointCell&  poly )
: PolygonCell( poly ),
  cell_type_         ( poly.cell_type_ ),
  cell_category_     ( poly.cell_category_ ),
  num_quad_faces_    ( poly.num_quad_faces_),
  num_tri_faces_     ( poly.num_tri_faces_ ),
  meshing_cycle_     ( poly.meshing_cycle_ ),
  wells_             ( poly.wells_ ),
  well_face_org_     ( poly.well_face_org_ ),
  well_face_dst_     ( poly.well_face_dst_ ),
  pillar_nodes_      ( poly.pillar_nodes_ ),
  pillar_nodes_order_( poly.pillar_nodes_order_ )
{
}




CornerPointCell& CornerPointCell::operator=( const CornerPointCell&  poly )
{
    if ( &poly != this )
    {
        PolygonCell::operator =( poly );
        cell_type_          = poly.cell_type_;
        cell_category_      = poly.cell_category_;
        num_quad_faces_     = poly.num_quad_faces_;
        num_tri_faces_      = poly.num_tri_faces_;

        meshing_cycle_      = poly.meshing_cycle_;

        wells_              = poly.wells_;
        well_face_org_      = poly.well_face_org_;
        well_face_dst_      = poly.well_face_dst_;

        pillar_nodes_       = poly.pillar_nodes_;
        pillar_nodes_order_ = poly.pillar_nodes_order_;
      }
    return *this;
}

CornerPointCell
::~CornerPointCell()
{
}

/**
   function replaces method in base class
     
   @todo remove argument because it has no impact on number of faces returned.
*/
size_t CornerPointCell::GetNumFaces( size_t ) const
 {
    return num_quad_faces_ + num_tri_faces_;
 }


// TYPE OF ELEMENT ( DIMENSION )

bool CornerPointCell::IsVolume() const
{
    return ( ( cell_category_ == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) || ( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING ) );
}

bool CornerPointCell::IsSurface() const
{
    return ( cell_category_ == CORNER_POINT_CELL_DIM_2 );
}

bool CornerPointCell::IsLine() const
{
    return ( cell_category_ == CORNER_POINT_CELL_DIM_1 );
}

bool CornerPointCell::IsPointElement() const
{
    return ( cell_category_ == CORNER_POINT_CELL_DIM_0 );
}

// TYPE OF ELEMENT ( TOPOLOGY )

bool CornerPointCell::IsVolumetricNonOverlapingElement() const
{
    return ( cell_category_ == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING );
}

bool CornerPointCell::IsVolumetricOverlapingElement() const
{
    return ( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING );
}

bool CornerPointCell::IsLowDimensionalElement() const
{
    return ( ( cell_category_ != CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) && ( cell_category_ != CORNER_POINT_CELL_DIM_3_OVERLAPPING ) );
}

// WELLS

size_t CornerPointCell::GetNumWellSections() const
{
    return wells_.size();
}

size_t CornerPointCell::GetNumWellSectionNodes( size_t eid ) const
{
    return wells_[eid].second.size();
}

const csmp::Point<3U>& CornerPointCell::GetWellSectionNode( size_t eid, size_t nid ) const
{
    return wells_[eid].second[nid]->GetPoint();
}

size_t CornerPointCell::GetWellSectionNodeGlobalId( size_t eid, size_t nid ) const
{
    return wells_[eid].second[nid]->GetIdx();
}

// PILLARS

size_t CornerPointCell::GetPillarNodeGlobalIdOriginalOrder( size_t nid ) const
{
    return pillar_nodes_[nid]->GetIdx();
}

size_t CornerPointCell::GetPillarNodeGlobalIdCustomOrder( size_t nid ) const
{
    return pillar_nodes_[nid]->GetIdx();
}

const csmp::Point<3U>& CornerPointCell::GetPillarPointOriginalOrder( size_t nid ) const
{
    return pillar_nodes_[nid]->GetPoint();
}

const csmp::Point<3U>& CornerPointCell::GetPillarPointCustomOrder( size_t nid ) const
{
    return pillar_nodes_[ pillar_nodes_order_[ nid ] ]->GetPoint();
}

GridNode* CornerPointCell::GetPillarNodeOriginalOrder( size_t nid )
{
    return pillar_nodes_[nid];
}

GridNode* CornerPointCell::GetPillarNodeCustomOrder( size_t nid )
{
    return pillar_nodes_[ pillar_nodes_order_[ nid ] ];
}


void CornerPointCell::AssignPillarNodesOrder()
{
    /// reorder nodes
    pillar_nodes_order_.clear();
    pillar_nodes_order_.resize(8);
    pillar_nodes_order_[0] = 4;
    pillar_nodes_order_[1] = 6;
    pillar_nodes_order_[2] = 7;
    pillar_nodes_order_[3] = 5;
    pillar_nodes_order_[4] = 0;
    pillar_nodes_order_[5] = 2;
    pillar_nodes_order_[6] = 3;
    pillar_nodes_order_[7] = 1;
}

void CornerPointCell::AssignNode( size_t local_nid, const csmp::Point<3U>& pt )
{
    assert( local_nid < grid_->GetNumNodes() );
    assert( local_nid < pillar_nodes_.size() );
    pillar_nodes_[local_nid] = this->grid_->AddNode(pt);
}

void CornerPointCell::CheckFaceOrder()
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    std::set<uint32_t> order;
    const size_t num_faces( 6 );
    for( size_t fid = 0; fid<num_faces; ++fid )
        order.insert( this->GetOriginalFaceId(fid) );
    if( order.size() != num_faces )
        csmp_error.Note( csmp::ERROR, "CornerPointCell::CheckFaceOrder()",
                           "Face order is not unique!!!" );
}

void CornerPointCell::CheckNodeOrder()
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );
    std::set<uint32_t> order;
    for( size_t nid = 0; nid<this->GetNumNodes(); ++nid )
        order.insert( this->GetOriginalNodeId(nid) );
    if( order.size() != this->GetNumNodes() )
        csmp_error.Note( csmp::ERROR, "CornerPointCell::CheckNodeOrder()",
                           "Node order is not unique!!!" );
}

void CornerPointCell::AddCornerPointWellCell( size_t face_id_org, size_t face_id_dst )
{
    well_face_org_ = face_id_org;
    well_face_dst_ = face_id_dst;
}



void CornerPointCell::InitializeCornerPointCell()
{
    std::map<csmp::Point<3U>,size_t> cell_nodes;
    InitializeNodes( cell_nodes );
    InitializeFaces();

    if ( this->GetNumNodes() == 8U )
      {
          InitializePoly8Element();
          return;
      }

    /// treat degenerate cells

    /// faces
    std::vector<std::vector<csmp::Point<3U> > >    face_ordered_points;
    std::vector<std::set<csmp::Point<3U> > >       face_unique_points;
    std::vector<uint32_t>                            tri_faces;
    std::vector<uint32_t>                            quad_faces;
    std::map<csmp::Point<3U>,std::vector<uint32_t> > point_tri_faces;
    std::map<csmp::Point<3U>,std::vector<uint32_t> > point_quad_faces;

    /// edges
    std::vector<std::vector<csmp::Point<3U> > >    edge_ordered_points;
    std::vector<std::set<csmp::Point<3U> > >       edge_unique_points;
    std::map<csmp::Point<3U>,std::vector<uint32_t> > point_edges;

    /// unique edges and faces
    std::set<std::set<csmp::Point<3U> > >          unique_faces;
    std::set<std::set<csmp::Point<3U> > >          unique_edges;

    GetTopologycalInfo( cell_nodes,
                        unique_faces, face_ordered_points, face_unique_points, tri_faces, quad_faces,
                        point_tri_faces, point_quad_faces,
                        unique_edges, edge_ordered_points, edge_unique_points,
                        point_edges );

    switch ( this->GetNumNodes() ) {
    case 7U:
        InitializePoly7Element( cell_nodes,
                                face_ordered_points, tri_faces, quad_faces,
                                point_tri_faces, point_quad_faces,
                                edge_ordered_points,
                                point_edges );
        break;
    case 6U:
        InitializePoly6Element( cell_nodes,
                                face_ordered_points, tri_faces, quad_faces,
                                point_tri_faces, point_quad_faces,
                                edge_ordered_points,
                                point_edges );
        break;
    case 5U:
        InitializePoly5Element( cell_nodes,
                                face_ordered_points, tri_faces, quad_faces,
                                point_tri_faces, point_quad_faces,
                                edge_ordered_points,
                                point_edges );
        break;
    case 4U:
        InitializePoly4Element( tri_faces, quad_faces );
        break;
    case 3U:
        InitializePoly3Element();
        break;
    case 2U:
        InitializePoly2Element();
        break;
    case 1U:
        InitializePoly1Element();
        break;
    default:
        break;
    }

} // end InitializeCornerPointCell




// RESOLVE MESH CONFLICTS

void CornerPointCell
::ResolveCornerPointCellConflicts( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges )
{
    if( ( this->GetNumNodes() < 8 ) && ( ( cell_category_ == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) || ( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING ) ) )
    //if( ( this->GetNumNodes() < 8 ) && ( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING ) )
        AddQuadFaceCentroids(additional_points,additional_edges);
}


// PROCESS CELLS

void CornerPointCell
::ProcessCornerPointCell( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    switch( this->GetNumNodes() )
    {
    case 8U:
        ProcessPoly8Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 7U:
        ProcessPoly7Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 6U:
        ProcessPoly6Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 5U:
        ProcessPoly5Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 4U:
        ProcessPoly4Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 3U:
        ProcessPoly3Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 2U:
        ProcessPoly2Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 1U:
        ProcessPoly1Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    default:
        break;
    }
}


void CornerPointCell::ProcessCornerPointRegularCell( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                          std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    if( this->GetNumNodes() == 8U )
        ProcessPoly8Element( additional_points, additional_edges, well_element, tetra_mesh );
}


void CornerPointCell::ProcessCornerPointDegenerateCell( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                             std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    switch( this->GetNumNodes() )
    {
    case 7U:
        ProcessPoly7Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 6U:
        ProcessPoly6Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 5U:
        ProcessPoly5Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 4U:
        ProcessPoly4Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 3U:
        ProcessPoly3Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 2U:
        ProcessPoly2Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    case 1U:
        ProcessPoly1Element( additional_points, additional_edges, well_element, tetra_mesh );
        break;
    default:
        break;
    }
}


void CornerPointCell
::ProcessCornerPointDegenerateVolumetricCell( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    if( ( cell_category_ == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) || ( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING ) )
    {
        switch( this->GetNumNodes() )
        {
        case 7U:
            ProcessPoly7Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 6U:
            ProcessPoly6Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 5U:
            ProcessPoly5Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 4U:
            ProcessPoly4Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        default:
            break;
        }
    }
}


void CornerPointCell
::ProcessCornerPointDegenerateVolumetricNonOverlappingCell( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges,
                                                            bool well_element, bool tetra_mesh )
{
    if( cell_category_ == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING )
    {
        switch( this->GetNumNodes() )
        {
        case 7U:
            ProcessPoly7Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 6U:
            ProcessPoly6Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 5U:
            ProcessPoly5Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 4U:
            ProcessPoly4Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        default:
            break;
        }
    }
}

void CornerPointCell::ProcessCornerPointDegenerateVolumetricOverlappingCell( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges,
                                                                                  bool well_element, bool tetra_mesh )
{
    if( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING )
    {
        switch( this->GetNumNodes() )
        {
        case 6U:
            ProcessPoly6Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 5U:
            ProcessPoly5Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        default:
            break;
        }
    }
}

void CornerPointCell
::ProcessCornerPointDegenerateLowDimensionalCell( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    if( ( cell_category_ != CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) && ( cell_category_ != CORNER_POINT_CELL_DIM_3_OVERLAPPING ) )
    {
        switch( this->GetNumNodes() )
        {
        case 4U:
            ProcessPoly4Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 3U:
            ProcessPoly3Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 2U:
            ProcessPoly2Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        case 1U:
            ProcessPoly1Element( additional_points, additional_edges, well_element, tetra_mesh );
            break;
        default:
            break;
        }
    }
}


// INITIALIZATION

void CornerPointCell
::InitializeNodes( std::map<csmp::Point<3U>,size_t>& cell_nodes )
{
    /// general hexa fem info
    const size_t hexa_fem_nodes( 8U );

    /// assign pillar nodes order
    this->AssignPillarNodesOrder();

    /// initial numbering of nodes
    cell_nodes.clear();
    size_t current_id = 0;
    std::pair<typename std::map<csmp::Point<3U>,size_t>::iterator,bool> unit;
    for( size_t nid = 0; nid < hexa_fem_nodes; ++nid )
    {
        unit= cell_nodes.insert( std::make_pair( GetPillarPointCustomOrder(nid), current_id ) );
        if( unit.second == true )
        {
            this->nodes_.push_back( pillar_nodes_[ pillar_nodes_order_[ nid ] ] );
            current_id++;
        }
    }

    /// reorder nodes
    this->InitializeNodeOrder( this->nodes_.size() );
}



void CornerPointCell::InitializeFaces( )
{
    const size_t num_faces( 6U );

    /// reorder nodes
    this->InitializeFaceOrder( num_faces );

    /// check correspondance of corner point cell enums to original face order
    assert( CORNER_POINT_CELL_FACE_Zplus   == 0 );
    assert( CORNER_POINT_CELL_FACE_Yminus  == 1 );
    assert( CORNER_POINT_CELL_FACE_Xminus  == 2 );
    assert( CORNER_POINT_CELL_FACE_Yplus   == 3 );
    assert( CORNER_POINT_CELL_FACE_Xplus   == 4 );
    assert( CORNER_POINT_CELL_FACE_Zminus  == 5 );
}


// PROCESSING CERTAIN CELLS

void CornerPointCell
::GetTopologycalInfo(   std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                        std::set<std::set<csmp::Point<3U> > >&          unique_faces,
                        std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                        std::vector<std::set<csmp::Point<3U> > >&       face_unique_points,
                        std::vector<uint32_t>&                             tri_faces,
                        std::vector<uint32_t>&                             quad_faces,
                        std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_tri_faces,
                        std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_quad_faces,
                        std::set<std::set<csmp::Point<3U> > >&          unique_edges,
                        std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                        std::vector<std::set<csmp::Point<3U> > >&       edge_unique_points,
                        std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_edges)
{
    /// hexa cell info
    const size_t quad_fem_nodes( 4U  );
    const size_t hexa_fem_nodes( 8U  );
    const size_t hexa_fem_faces( 6U  );
    const size_t hexa_fem_edges( 12U );

    /// topological info
    /// faces
    face_ordered_points.clear();
    face_ordered_points.resize( hexa_fem_nodes );
    face_unique_points.clear();
    face_unique_points.resize( hexa_fem_faces );
    tri_faces.clear();
    quad_faces.clear();
    point_tri_faces.clear();
    point_quad_faces.clear();

    /// edges
    edge_ordered_points.clear();
    edge_ordered_points.resize( hexa_fem_edges );
    edge_unique_points.clear();
    edge_unique_points.resize( hexa_fem_edges );
    point_edges.clear();

    /// unique edges and faces
    std::pair<typename std::set<csmp::Point<3U> >::iterator,bool >           pit;
    std::pair<typename std::set<std::set<csmp::Point<3U> > >::iterator,bool> upit;

    /// point sets
    std::set<csmp::Point<3U> >  points_set;

    /// initialize faces

    /// bottom face
    this->ClearPolygonFaceNodes(0);
    pit = face_unique_points[ 0 ].insert( GetPillarPointCustomOrder( 0 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 0 ].push_back( GetPillarPointCustomOrder( 0 ) );
        this->AddPolygonFaceNode( 0, cell_nodes[ GetPillarPointCustomOrder( 0 ) ] );
    }
    pit = face_unique_points[ 0 ].insert( GetPillarPointCustomOrder( 3 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 0 ].push_back( GetPillarPointCustomOrder( 3 ) );
        this->AddPolygonFaceNode( 0, cell_nodes[ GetPillarPointCustomOrder( 3 ) ] );
    }
    pit = face_unique_points[ 0 ].insert( GetPillarPointCustomOrder( 2 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 0 ].push_back( GetPillarPointCustomOrder( 2 ) );
        this->AddPolygonFaceNode( 0, cell_nodes[ GetPillarPointCustomOrder( 2 ) ] );
    }
    pit = face_unique_points[ 0 ].insert( GetPillarPointCustomOrder( 1 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 0 ].push_back( GetPillarPointCustomOrder( 1 ) );
        this->AddPolygonFaceNode( 0, cell_nodes[ GetPillarPointCustomOrder( 1 ) ] );
    }
    if( face_unique_points[ 0 ].size() > 2U )
    {
        upit = unique_faces.insert( face_unique_points[ 0 ] );
        if( (*upit.first).size() == 4 )
        {
            if( upit.second == true )
                quad_faces.push_back( 0 );
            for( typename std::set<csmp::Point<3U> >::const_iterator
                 it = face_unique_points[ 0 ].begin(); it!=face_unique_points[ 0 ].end(); ++it )
                point_quad_faces[ *it ].push_back( 0 );
            //addQuadFace(this->faces_,0,
            //            this->GetPolygonFaceNode(0,0),
            //            this->GetPolygonFaceNode(0,1),
            //            this->GetPolygonFaceNode(0,2),
            //            this->GetPolygonFaceNode(0,3) );
        }
        else if( (*upit.first).size() == 3 )
        {
            if( upit.second == true )
                tri_faces.push_back( 0 );
            for( typename std::set<csmp::Point<3U> >::const_iterator
                 it = face_unique_points[ 0 ].begin(); it!=face_unique_points[ 0 ].end(); ++it )
                point_tri_faces[ *it ].push_back( 0 );
            //addTriFace(this->faces_,0,
            //           this->GetPolygonFaceNode(0,0),
            //           this->GetPolygonFaceNode(0,1),
            //           this->GetPolygonFaceNode(0,2) );
        }
    }
    else if( face_unique_points[ 0 ].size() == 2U )
    {
        //addBarFace(this->faces_,0,
        //           this->GetPolygonFaceNode(0,0),
        //           this->GetPolygonFaceNode(0,1) );
    }
    else
    {
        //addPointFace(this->faces_,0,
        //             this->GetPolygonFaceNode(0,0) );
    }

    /// lateral faces
    for( uint32_t fid=0U; fid<4; ++fid )
    {
        const size_t nid = (fid + 3)%quad_fem_nodes;
        this->ClearPolygonFaceNodes(fid+1);
        pit = face_unique_points[ fid + 1 ].insert( GetPillarPointCustomOrder( nid ) );
        if( pit.second == true )
        {
            face_ordered_points[ fid + 1 ].push_back( GetPillarPointCustomOrder( nid ) );
            this->AddPolygonFaceNode( fid+1, cell_nodes[ GetPillarPointCustomOrder( nid ) ] );
        }
        pit = face_unique_points[ fid + 1 ].insert( GetPillarPointCustomOrder( (nid + 1)%quad_fem_nodes ) );
        if( pit.second == true )
        {
            face_ordered_points[ fid + 1 ].push_back( GetPillarPointCustomOrder( (nid + 1)%quad_fem_nodes ) );
            this->AddPolygonFaceNode( fid+1, cell_nodes[ GetPillarPointCustomOrder( (nid + 1)%quad_fem_nodes ) ] );
        }
        pit = face_unique_points[ fid + 1 ].insert( GetPillarPointCustomOrder( ( (nid + 1)%quad_fem_nodes ) + quad_fem_nodes ) );
        if( pit.second == true )
        {
            face_ordered_points[ fid + 1 ].push_back( GetPillarPointCustomOrder( ( (nid + 1)%quad_fem_nodes ) + quad_fem_nodes ) );
            this->AddPolygonFaceNode( fid+1, cell_nodes[ GetPillarPointCustomOrder( ( (nid + 1)%quad_fem_nodes ) + quad_fem_nodes ) ] );
        }
        pit = face_unique_points[ fid + 1 ].insert( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
        if( pit.second == true )
        {
            face_ordered_points[ fid + 1 ].push_back( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
            this->AddPolygonFaceNode( fid+1, cell_nodes[ GetPillarPointCustomOrder( nid + quad_fem_nodes ) ] );
        }
        if( face_unique_points[ fid + 1 ].size() > 2U )
        {
            upit = unique_faces.insert( face_unique_points[ fid + 1 ] );
            if( (*upit.first).size() == 4 )
            {
                if( upit.second == true )
                    quad_faces.push_back( fid + 1 );
                for( typename std::set<csmp::Point<3U> >::const_iterator
                     it = face_unique_points[ fid + 1 ].begin(); it!=face_unique_points[ fid + 1 ].end(); ++it )
                    point_quad_faces[ *it ].push_back( fid + 1 );
                //addQuadFace(this->faces_,fid+1,
                //            this->GetPolygonFaceNode(fid+1,0),
                //            this->GetPolygonFaceNode(fid+1,1),
                //            this->GetPolygonFaceNode(fid+1,2),
                //            this->GetPolygonFaceNode(fid+1,3) );
            }
            else if( (*upit.first).size() == 3 )
            {
                if( upit.second == true )
                    tri_faces.push_back( fid + 1 );
                for( typename std::set<csmp::Point<3U> >::const_iterator
                     it = face_unique_points[ fid + 1 ].begin(); it!=face_unique_points[ fid + 1 ].end(); ++it )
                    point_tri_faces[ *it ].push_back( fid + 1 );
                //addTriFace(this->faces_,fid+1,
                //           this->GetPolygonFaceNode(fid+1,0),
                //           this->GetPolygonFaceNode(fid+1,1),
                //           this->GetPolygonFaceNode(fid+1,2) );
            }
        }
        else if( face_unique_points[ fid + 1 ].size() == 2U )
        {
            //addBarFace(this->faces_,fid+1,
            //           this->GetPolygonFaceNode(fid+1,0),
            //           this->GetPolygonFaceNode(fid+1,1) );
        }
        else
        {
            //addPointFace(this->faces_,fid+1,
            //             this->GetPolygonFaceNode(fid+1,0) );
        }
    }
    /// top face (4,5,6,7)
    this->ClearPolygonFaceNodes(5);
    pit = face_unique_points[ 5 ].insert( GetPillarPointCustomOrder( 4 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 5 ].push_back( GetPillarPointCustomOrder( 4 ) );
        this->AddPolygonFaceNode( 5, cell_nodes[ GetPillarPointCustomOrder( 4 ) ] );
    }
    pit = face_unique_points[ 5 ].insert( GetPillarPointCustomOrder( 5 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 5 ].push_back( GetPillarPointCustomOrder( 5 ) );
        this->AddPolygonFaceNode( 5, cell_nodes[ GetPillarPointCustomOrder( 5 ) ] );
    }
    pit = face_unique_points[ 5 ].insert( GetPillarPointCustomOrder( 6 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 5 ].push_back( GetPillarPointCustomOrder( 6 ) );
        this->AddPolygonFaceNode( 5, cell_nodes[ GetPillarPointCustomOrder( 6 ) ] );
    }
    pit = face_unique_points[ 5 ].insert( GetPillarPointCustomOrder( 7 ) );
    if( pit.second == true )
    {
        face_ordered_points[ 5 ].push_back( GetPillarPointCustomOrder( 7 ) );
        this->AddPolygonFaceNode( 5, cell_nodes[ GetPillarPointCustomOrder( 7 ) ] );
    }
    if( face_unique_points[ 5 ].size() > 2U )
    {
        upit = unique_faces.insert( face_unique_points[ 5 ] );
        if( (*upit.first).size() == 4 )
        {
            if( upit.second == true )
                quad_faces.push_back( 5 );
            for( typename std::set<csmp::Point<3U> >::const_iterator
                 it = face_unique_points[ 5 ].begin(); it!=face_unique_points[ 5 ].end(); ++it )
                point_quad_faces[ *it ].push_back( 5 );
            //addQuadFace(this->faces_,5,
            //            this->GetPolygonFaceNode(5,0),
            //            this->GetPolygonFaceNode(5,1),
            //            this->GetPolygonFaceNode(5,2),
            //            this->GetPolygonFaceNode(5,3) );
        }
        else if( (*upit.first).size() == 3 )
        {
            if( upit.second == true )
                tri_faces.push_back( 5 );
            for( typename std::set<csmp::Point<3U> >::const_iterator
                 it = face_unique_points[ 5 ].begin(); it!=face_unique_points[ 5 ].end(); ++it )
                point_tri_faces[ *it ].push_back( 5 );
            //addTriFace(this->faces_,5,
            //           this->GetPolygonFaceNode(5,0),
            //           this->GetPolygonFaceNode(5,1),
            //           this->GetPolygonFaceNode(5,2) );
        }
    }
    else if( face_unique_points[ 5 ].size() == 2U )
    {
        //addBarFace(this->faces_,5,
        //           this->GetPolygonFaceNode(5,0),
        //           this->GetPolygonFaceNode(5,1) );
    }
    else
    {
        //addPointFace(this->faces_,5,
        //             this->GetPolygonFaceNode(5,0) );
    }

    /// edges
    for( uint32_t nid = 0U; nid <quad_fem_nodes; ++nid )
    {
        /// bottom edge
        points_set.clear();
        points_set.insert( GetPillarPointCustomOrder( nid ) );
        points_set.insert( GetPillarPointCustomOrder( (nid + 1)%quad_fem_nodes ) );
        edge_unique_points[ nid ].insert( points_set.begin(), points_set.end() );
        if( points_set.size() == 2 )
        {
            edge_ordered_points[ nid ].push_back( GetPillarPointCustomOrder( nid ) );
            edge_ordered_points[ nid ].push_back( GetPillarPointCustomOrder( (nid + 1)%quad_fem_nodes ) );
            upit = unique_edges.insert( points_set );
            if( upit.second == true )
            {
                for( typename std::set<csmp::Point<3U> >::const_iterator
                     it = points_set.begin(); it != points_set.end(); ++it )
                    point_edges[ *it ].push_back( nid );
            }
        }

        /// top edge
        points_set.clear();
        points_set.insert( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
        points_set.insert( GetPillarPointCustomOrder( ((nid + 1)%quad_fem_nodes) + quad_fem_nodes ) );
        edge_unique_points[  nid + quad_fem_nodes ].insert( points_set.begin(), points_set.end() );
        if( points_set.size() == 2 )
        {
            edge_ordered_points[ nid + quad_fem_nodes ].push_back( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
            edge_ordered_points[ nid + quad_fem_nodes ].push_back( GetPillarPointCustomOrder( ((nid + 1)%quad_fem_nodes) + quad_fem_nodes ) );
            upit = unique_edges.insert( points_set );
            if( upit.second == true )
            {
                for( typename std::set<csmp::Point<3U> >::const_iterator
                     it = points_set.begin(); it != points_set.end(); ++it )
                    point_edges[ *it ].push_back( nid + quad_fem_nodes );
            }
        }

        /// lateral edge
        points_set.clear();
        points_set.insert( GetPillarPointCustomOrder( nid ) );
        points_set.insert( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
        edge_unique_points[ nid + hexa_fem_nodes ].insert( points_set.begin(), points_set.end() );
        if( points_set.size() == 2 )
        {
            edge_ordered_points[ nid + hexa_fem_nodes ].push_back( GetPillarPointCustomOrder( nid ) );
            edge_ordered_points[ nid + hexa_fem_nodes ].push_back( GetPillarPointCustomOrder( nid + quad_fem_nodes ) );
            upit = unique_edges.insert( points_set );
            if( upit.second == true )
            {
                for( typename std::set<csmp::Point<3U> >::const_iterator
                     it = points_set.begin(); it != points_set.end(); ++it )
                    point_edges[ *it ].push_back( nid + hexa_fem_nodes );
            }
        }
    }

    /// assign number of unique faces
    num_quad_faces_= quad_faces.size();
    num_tri_faces_ = tri_faces.size();
}



void CornerPointCell::AddQuadFaceCentroids( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges )
{
    std::vector<std::vector<uint32_t> > quad_face_nodes;
    for( uint32_t fid = 0U; fid<PolygonCell::GetNumFaces(); ++fid )
    {
        const size_t num_face_nodes( this->GetNumPolygonFaceNodes( fid ) );
        if( num_face_nodes == 4 )
        {
            std::vector<uint32_t> fnodes;
            for( size_t nid = 0; nid < num_face_nodes; ++nid )
                fnodes.push_back( this->GetPolygonFaceNodeLocalId(fid,nid) );
            quad_face_nodes.push_back(fnodes);
        }
    }
    if( cell_category_ == CORNER_POINT_CELL_DIM_3_OVERLAPPING )
    {
        GridNode* gn;
        addOverlappingCellFaceCentroids( this->grid_, additional_points, additional_edges, gn, this->nodes_, quad_face_nodes );
    }
    else
    {
        addCellFaceCentroids( this->grid_, additional_points, additional_edges, this->nodes_, quad_face_nodes );
    }
}

void CornerPointCell::FinishMeshing()
{
    meshing_cycle_ = -1U;

    const size_t num_faces( this->faces_.size() );
    std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > empty_faces;
    std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > > vec_temp_faces( num_faces, empty_faces );
    for( size_t fid = 0; fid < num_faces; ++fid )
        vec_temp_faces[ this->GetOriginalFaceId(fid) ] = this->faces_[ fid ];
    this->faces_.clear();
    this->faces_ = vec_temp_faces;
    return;
}


void CornerPointCell::AddWell( std::map<std::set<uint32_t>,GridNode*>& additional_points )
{
    GridNode* cgn( this->GetCellCentroidNode() );
    const size_t cgn_idx( this->GetCellCentroidGlobalId() );

    GridNode* org_gn;
    const size_t num_org_face_nodes( this->GetNumPolygonFaceNodes(well_face_org_) );
    if( num_org_face_nodes == 4 )
    {
        doesPointExist( additional_points,org_gn,
                      this->GetPolygonFaceNode(well_face_org_,0),
                      this->GetPolygonFaceNode(well_face_org_,1),
                      this->GetPolygonFaceNode(well_face_org_,2),
                      this->GetPolygonFaceNode(well_face_org_,3) );
    }
    else if( num_org_face_nodes == 3 )
    {
        doesPointExist( additional_points,org_gn,
                      this->GetPolygonFaceNode(well_face_org_,0),
                      this->GetPolygonFaceNode(well_face_org_,1),
                      this->GetPolygonFaceNode(well_face_org_,2) );
    }
    else if( num_org_face_nodes == 2 )
    {
        doesPointExist( additional_points,org_gn,
                      this->GetPolygonFaceNode(well_face_org_,0),
                      this->GetPolygonFaceNode(well_face_org_,1) );
    }
    else
        org_gn = this->GetPolygonFaceNode(well_face_org_,0);

// SKM FIX:    if( ( cgn_idx != NULL_IDX ) && ( cgn_idx != org_gn->GetIdx() ) )
    if( ( cgn_idx != UINT_MAX ) && ( cgn_idx != org_gn->GetIdx() ) )
        addBarCell( wells_, org_gn, cgn );

    GridNode* dst_gn;
    const size_t num_dst_face_nodes( this->GetNumPolygonFaceNodes(well_face_dst_) );
    if( num_dst_face_nodes ==4 )
    {
        doesPointExist( additional_points,dst_gn,
                      this->GetPolygonFaceNode(well_face_dst_,0),
                      this->GetPolygonFaceNode(well_face_dst_,1),
                      this->GetPolygonFaceNode(well_face_dst_,2),
                      this->GetPolygonFaceNode(well_face_dst_,3) );
    }
    else if( num_dst_face_nodes == 3 )
    {
        doesPointExist( additional_points,dst_gn,
                      this->GetPolygonFaceNode(well_face_dst_,0),
                      this->GetPolygonFaceNode(well_face_dst_,1),
                      this->GetPolygonFaceNode(well_face_dst_,2) );
    }
    else if( num_dst_face_nodes == 2 )
    {
        doesPointExist( additional_points,dst_gn,
                      this->GetPolygonFaceNode(well_face_dst_,0),
                      this->GetPolygonFaceNode(well_face_dst_,1) );
    }
    else
        dst_gn = this->GetPolygonFaceNode(well_face_dst_,0);

// SKM FIX    if( ( cgn_idx != NULL_IDX ) && ( cgn_idx != dst_gn->GetIdx() ) )
    if( ( cgn_idx != UINT_MAX ) && ( cgn_idx != dst_gn->GetIdx() ) )
        addBarCell( wells_, cgn, dst_gn );
}

void CornerPointCell::ProcessPolyWellElement( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges )
{
    const bool tetra_mesh( true );
    const bool add_edge_centroids( true );
    const bool add_face_centroids( true );
    const bool add_cell_centroids( true );
    const bool is_volumetric_element( IsVolume() );
    processCellWithCentroids( add_edge_centroids, add_face_centroids, add_cell_centroids, is_volumetric_element, tetra_mesh,
                              this->grid_, additional_points, additional_edges,
                              this->elements_,this->faces_,this->extra_nodes_,
                              this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
    /// add well elements
    AddWell(additional_points);
}

void CornerPointCell::ProcessPolyNonWellElement( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges,
                                                      bool add_edge_centroids, bool add_face_centroids, bool add_cell_centroids )
{
    const bool tetra_mesh( true );
    const bool do_not_add_cell_centroids( false );
    const bool is_volumetric_element( IsVolume() );
    if( this->GetNumNodes() == 8U )
        processCellWithCentroids( add_edge_centroids, add_face_centroids, add_cell_centroids, is_volumetric_element, tetra_mesh,
                                  this->grid_, additional_points, additional_edges,
                                  this->elements_,this->faces_,this->extra_nodes_,
                                  this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
    else
        processCellWithCentroids( add_edge_centroids, add_face_centroids, do_not_add_cell_centroids, is_volumetric_element, tetra_mesh,
                                  this->grid_, additional_points, additional_edges,
                                  this->elements_,this->faces_,this->extra_nodes_,
                                  this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
}

void CornerPointCell::ProcessPolyElement( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element )
{
    if( ( meshing_cycle_ == 0 ) && well_element )
    {
        /// cell with well
        ProcessPolyWellElement(additional_points,additional_edges);

        FinishMeshing();
        return;
    }
    else if( !well_element )
    {
        if( meshing_cycle_ == 0 )
        {
            /// cell with meshed edges
            if( containRemeshedEdges( additional_points, this->nodes_in_custom_order_, this->face_nodes_in_custom_order_ ) )
            {
                const bool add_edge_centroids( false );
                const bool add_face_centroids( true );
                const bool add_cell_centroids( true );
                ProcessPolyNonWellElement(additional_points,additional_edges,add_edge_centroids,add_face_centroids,add_cell_centroids);
                FinishMeshing();
                return;
            }
            else
                meshing_cycle_ += 1U;
        }
        else if( meshing_cycle_ == 1 )
        {
            /// cell with meshed faces
            if( containRemeshedFaces( additional_points, this->nodes_in_custom_order_, this->face_nodes_in_custom_order_ ) )
            {
                const bool add_edge_centroids( false );
                const bool add_face_centroids( false );
                const bool add_cell_centroids( true );
                ProcessPolyNonWellElement(additional_points,additional_edges,add_edge_centroids,add_face_centroids,add_cell_centroids);
                FinishMeshing();
                return;
            }
            else
                meshing_cycle_ += 1U;
        }
        else if( meshing_cycle_ == 2 )
        {
            /// edges and faces are not meshed
            const bool add_edge_centroids( false );
            const bool add_face_centroids( false );
            const bool add_cell_centroids( true );
            ProcessPolyNonWellElement(additional_points,additional_edges,add_edge_centroids,add_face_centroids,add_cell_centroids);
            FinishMeshing();
            return;
        }
    }
}



/**
    InitializePoly8Element():
    Poly8: sketch

      ______ X
     |\
     | \ Y
     |
     Z

     hexahedron:

       |              |
       |              |
       |              |
       |              |
     7 x____|_________x 6  |
       |\   |         |\   |
       | \  |         | \  |
       |  \ |         |  \ |
       |   4x_________|___\x 5
       |    |         |    |
       x____|_________x    |
     3 |\   |        2|\   |
       | \  |         | \  |
       |  \ |         |  \ |
       |   \x_________|___\x
       |  0 |         |    | 1
            |              |
            |              |
            |              |
            |              |
*/

void CornerPointCell::InitializePoly8Element()
{
    /// hexa cell info
    const size_t hexa_fem_nodes( 8U  );
    const size_t quad_fem_nodes( 4U  );
    const size_t hexa_fem_faces( 6U  );

    /// (1) 5 or 6 tetrahedrons
    /// observed and tested
    cell_type_     = CORNER_POINT_CELL_8_HEXAHEDRON;
    cell_category_ = CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING;
    num_quad_faces_= 6;
    num_tri_faces_ = 0;

    /// initialize faces
    for( size_t fid=0; fid<hexa_fem_faces; ++fid )
        this->ResizePolygonFaceNodes( fid, quad_fem_nodes );

    /// bottom face0 ( 0321 )
    this->AssignPolygonFaceNode( 0, 0, 0 );
    this->AssignPolygonFaceNode( 0, 1, 3 );
    this->AssignPolygonFaceNode( 0, 2, 2 );
    this->AssignPolygonFaceNode( 0, 3, 1 );
    /// left face1   ( 3047 )
    this->AssignPolygonFaceNode( 1, 0, 3 );
    this->AssignPolygonFaceNode( 1, 1, 0 );
    this->AssignPolygonFaceNode( 1, 2, 4 );
    this->AssignPolygonFaceNode( 1, 3, 7 );
    /// front face2  ( 0154 )
    this->AssignPolygonFaceNode( 2, 0, 0 );
    this->AssignPolygonFaceNode( 2, 1, 1 );
    this->AssignPolygonFaceNode( 2, 2, 5 );
    this->AssignPolygonFaceNode( 2, 3, 4 );
    /// right face3  ( 1265 )
    this->AssignPolygonFaceNode( 3, 0, 1 );
    this->AssignPolygonFaceNode( 3, 1, 2 );
    this->AssignPolygonFaceNode( 3, 2, 6 );
    this->AssignPolygonFaceNode( 3, 3, 5 );
    /// back face4  ( 2376 )
    this->AssignPolygonFaceNode( 4, 0, 2 );
    this->AssignPolygonFaceNode( 4, 1, 3 );
    this->AssignPolygonFaceNode( 4, 2, 7 );
    this->AssignPolygonFaceNode( 4, 3, 6 );
    /// top face5   ( 4567 )
    this->AssignPolygonFaceNode( 5, 0, 4 );
    this->AssignPolygonFaceNode( 5, 1, 5 );
    this->AssignPolygonFaceNode( 5, 2, 6 );
    this->AssignPolygonFaceNode( 5, 3, 7 );

    /// add faces
    //for( size_t fid=0; fid<hexa_fem_faces; ++fid )
    //    addQuadFace( this->faces_, fid, this->GetPolygonFaceNode( fid, 0 ), this->GetPolygonFaceNode( fid, 1 ), this->GetPolygonFaceNode( fid, 2 ), this->GetPolygonFaceNode( fid, 3 ) );

    /// node order
    for( size_t nid = 0; nid < hexa_fem_nodes; ++nid )
        this->AssignNodeOrder( nid, nid );

    /// face order
    for( size_t fid = 0; fid < hexa_fem_faces; ++fid )
        this->AssignFaceOrder( fid, fid );

    CheckFaceOrder();
    CheckNodeOrder();
}



void CornerPointCell::ProcessPoly8Element( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                std::set<std::set<uint32_t> >& additional_edges,
                                                bool well_element, bool tetra_mesh )
{
    /// (1) 5 or 6 tetrahedrons
    if( meshing_cycle_ < 2 )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( meshing_cycle_ == 2U )
    {
        /// edges and faces are not meshed
        processHexahedronElement(tetra_mesh,
                                 this->grid_, additional_points, additional_edges,
                                 this->elements_,this->faces_,this->extra_nodes_,
                                 this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}

/**
    InitializePoly7Element():
    Poly7: sketch

      ______ X
     |\
     | \ Y
     |
     Z

     augmented prism:

       |              |
       |              |
       |              |
       |              |
     5 x____|_________x 4  |
       |\   |         |    |
       | \  |         |    |
       |  \ |         |    |
       |  6 x         |    |
       |    |         |    |
       x____|_________x    |
     2 |\   |        1|\   |
       | \  |         | \  |
       |  \ |         |  \ |
       |   \x_________|___\x
       |  3 |         |    | 0
            |              |
            |              |
            |              |
            |              |

*/

void CornerPointCell
::InitializePoly7Element(  const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                           const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                           const std::vector<uint32_t>&                            tri_faces,
                           const std::vector<uint32_t>&                            quad_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_tri_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_quad_faces,
                           const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_edges )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// point sets
    std::set<csmp::Point<3U> > points_set;

    /// node id
    size_t current_id;

    /// hexa cell info
    const size_t quad_fem_nodes( 4U  );

    /// (1) 4 or 5 tetrahedrons

    if( ( num_quad_faces_ == 4 ) && ( num_tri_faces_  == 2 ) )
    {
        /// observed and tested
        cell_type_     = CORNER_POINT_CELL_7_AUGMENTED_PRISM;
        cell_category_ = CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING;

        /// renumbering nodes
        /// 0 - point where 2 triangles meet
        current_id = 0;
        points_set.clear();
        for( size_t nid =0; nid < this->GetNumNodes(); ++nid )
            if( point_edges.at( this->GetPointOriginalOrder( nid ) ).size() == 4 )
            {
                points_set.insert( this->GetPointOriginalOrder( nid ) );
                this->AssignNodeOrder( current_id++, nid );
                break;
            }
        assert( current_id == 1 );
        /// 1,2,3 - bottom quad
        for( size_t nid = 0; nid < 4U; ++nid )
            if( points_set.find( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ nid ] ) != points_set.end() )
            {
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] ) );
                break;
            }
        assert( current_id == 4 );
        /// 4,5,6 - top quad
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(1) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(1) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 5 );
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(2) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(2) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 6 );
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(3) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(3) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 7 );
        assert( points_set.size() == 7 );

        /// reorder faces
        /// bottom quad face ( 0123 )
        this->AssignFaceOrder( 0, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] );
        /// left quad face   ( 1452 )
        this->AssignFaceOrder( 1, ( point_quad_faces.at( this->GetPointCustomOrder(1) )[ 0 ] != this->GetOriginalFaceId(0) ? point_quad_faces.at( this->GetPointCustomOrder(1) )[ 0 ] : point_quad_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ) );
        /// right quad face  ( 2563 )
        this->AssignFaceOrder( 2, ( point_quad_faces.at( this->GetPointCustomOrder(3) )[ 0 ] != this->GetOriginalFaceId(0) ? point_quad_faces.at( this->GetPointCustomOrder(3) )[ 0 ] : point_quad_faces.at( this->GetPointCustomOrder(3) )[ 1 ] ) );
        /// top quud face    ( 0654 )
        this->AssignFaceOrder( 3, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 1 ] );
        /// left tri face    ( 014 )
        this->AssignFaceOrder( 4, point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] );
        /// right tri face   ( 036 )
        this->AssignFaceOrder( 5, point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] );
    }
    else
        csmp_error.Note( csmp::FATAL_ERROR,
                           "CornerPointCell::InitializePoly7Element():", "Unpredictable case for 7 point cell!!");
    CheckFaceOrder();
    CheckNodeOrder();
}

void CornerPointCell::ProcessPoly7Element( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges,
                                                bool well_element, bool tetra_mesh )
{
    /// (1) 4 or 5 tetrahedrons
    if( meshing_cycle_ < 2 )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( meshing_cycle_ == 2U )
    {
        /// edges and faces are not meshed
        processDegenerateHexahedronElement(tetra_mesh,
                                           this->grid_, additional_points, additional_edges,
                                           this->elements_,this->faces_,this->extra_nodes_,
                                           this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}

/**
    InitializePoly6Element():
    Poly6: sketch

      ______ X
     |\
     | \ Y
     |
     Z

     wedge:

     |          |
     |          |
     |          | 5
     |          x
     |   |     /|\   |
     |   |    / | \  |
     |   |   /  |  \ |
     |   |  /   |   \x 2
     |   | /    |  / |
     |   |/     | /  |
     |   |      |/   |
     |  /|      /    |
     | / |     /|    |
     |/  |    / |    |
   3 x___|___/__x 4  |
     |\  |  /    \   |
     | \ | /      \  |
     |  \|/        \ |
     |   x__________x|
       0 |           | 1
         |           |
         |           |

    octahedron:

       |              |
     5 x              |
       |              |
       |    |         |    |
       |    |         |    |
       |    |         |    x 4
       x____|_________x    |
     3 |\   |        2|\   |
       | \  |         | \  |
       |  \ |         |  \ |
       |   \x_________|___\x
       |  0 |         |    | 1
            |              |
            |              |
            |              |

*/

void CornerPointCell
::InitializePoly6Element(  const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                           const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                           const std::vector<uint32_t>&                            tri_faces,
                           const std::vector<uint32_t>&                            quad_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_tri_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_quad_faces,
                           const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_edges )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// point sets
    std::set<csmp::Point<3U> >      points_set;

    /// node id
    size_t current_id;

    /// hexa cell info
    const size_t quad_fem_nodes( 4U  );

    /// (1) wedge ( 3 tetrahedrons )
    /// (2) octahedron ( 2 tetrahedrons )

    /// 1st case: wedge
    if( ( num_quad_faces_ == 3 ) && ( num_tri_faces_ == 2 ) )
    {
        /// observed and tested
        cell_type_     = CORNER_POINT_CELL_6_WEDGE;
        cell_category_ = CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING;

        /// renumbering nodes
        /// 0,1,2 - front triangle
        current_id = 0;
        points_set.clear();
        for( typename std::vector<csmp::Point<3U> >::const_iterator
             it = face_ordered_points.at( tri_faces[ 0 ] ).begin(); it != face_ordered_points.at( tri_faces[ 0 ] ).end(); ++it )
            if( points_set.find( *it ) == points_set.end() )
            {
                points_set.insert( *it );
                this->AssignNodeOrder( current_id++, cell_nodes.at( *it ) );
            }
        assert( current_id == 3 );
        /// 3,4,5 - back triangle
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(0) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(0) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 4 );
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(1) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(1) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 5 );
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(2) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(2) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
         assert( current_id == 6 );
         assert( points_set.size() == 6 );

         /// reorder faces
         std::set<uint32_t> faces;
         for( uint32_t fid{0U}; fid<6; ++fid )
             faces.insert(fid);
         /// bottom quad face ( 0253 )
         this->AssignFaceOrder( 0, ( ( point_quad_faces.at( this->GetPointCustomOrder(2) )[ 0 ] == point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] ) ? point_quad_faces.at( this->GetPointCustomOrder(2) )[ 0 ] : ( point_quad_faces.at( this->GetPointCustomOrder(2) )[ 0 ] == point_quad_faces.at( this->GetPointCustomOrder(0) )[ 1 ] ) ? point_quad_faces.at( this->GetPointCustomOrder(2) )[ 0 ] : point_quad_faces.at( this->GetPointCustomOrder(2) )[ 1 ] ) );
         faces.erase(this->GetOriginalFaceId(0));
         /// front quad face  ( 0143 )
         this->AssignFaceOrder( 1, ( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] != this->GetOriginalFaceId(0) ? point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] : point_quad_faces.at( this->GetPointCustomOrder(0) )[ 1 ] ) );
         faces.erase(this->GetOriginalFaceId(1));
         /// back quad face   ( 1452 )
         this->AssignFaceOrder( 2, ( point_quad_faces.at( this->GetPointCustomOrder(1) )[ 0 ] != this->GetOriginalFaceId(1) ? point_quad_faces.at( this->GetPointCustomOrder(1) )[ 0 ] : point_quad_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ) );
         faces.erase(this->GetOriginalFaceId(2));
         /// left tri face   ( 012 )
         this->AssignFaceOrder( 3, tri_faces[ 0 ] );
         faces.erase(this->GetOriginalFaceId(3));
         /// right tri face  ( 354 )
         this->AssignFaceOrder( 4, tri_faces[ 1 ] );
         faces.erase(this->GetOriginalFaceId(4));
         /// edge face ( 03 )
         size_t fid = 4;
         for(std::set<uint32_t>::const_iterator fit = faces.begin(); fit != faces.end(); ++fit )
             this->AssignFaceOrder( ++fid, *fit );
    }
    /// 2nd case: octahedron( 4 or 2 tetrahedrons )
    else if( ( num_quad_faces_ == 2 ) && ( num_tri_faces_ == 4 ) )
    {
        /// observed and tested
        cell_type_     = CORNER_POINT_CELL_6_OCTAHEDRON;
        cell_category_ = CORNER_POINT_CELL_DIM_3_OVERLAPPING;

        /// renumbering nodes
        /// 0 - node where quad's collapse
        current_id = 0;
        points_set.clear();
        for( size_t nid =0; nid < this->GetNumNodes(); ++nid )
            if( point_edges.at( this->GetPointOriginalOrder( nid ) ).size() == 4 )
            {
                points_set.insert( this->GetPointOriginalOrder( nid ) );
                this->AssignNodeOrder( current_id++, nid );
                break;
            }
        assert( current_id == 1 );
        /// 1,2,3 - bottom quad
        for( size_t nid = 0; nid < 4U; ++nid )
            if( points_set.find( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ nid ] ) != points_set.end() )
            {
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] ) );
                break;
            }
        assert( current_id == 4 );
        /// 4,5 - top quad nodes
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(1) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(1) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 5 );
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(3) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(3) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 6 );
        assert( points_set.size() == 6 );

        /// reorder faces
        /// bottom quad face ( 0321 )
        this->AssignFaceOrder( 0, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] );
        /// top quad face ( 0425 )
        this->AssignFaceOrder( 1, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 1 ] );
        /// tri face   ( 014 )
        this->AssignFaceOrder( 2, ( point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] == point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] ? point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] : ( point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] == point_tri_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ? point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(0) )[ 1 ] ) ) );
        /// tri face   ( 124 )
        this->AssignFaceOrder( 3, ( point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] != this->GetOriginalFaceId(2) ? point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ) );
        /// tri face  ( 235 )
        this->AssignFaceOrder( 4, ( point_tri_faces.at( this->GetPointCustomOrder(2) )[ 0 ] != this->GetOriginalFaceId(3) ? point_tri_faces.at( this->GetPointCustomOrder(2) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(2) )[ 1 ] ) );
        /// tri face  ( 305 )
        this->AssignFaceOrder( 5, ( point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] != this->GetOriginalFaceId(4) ? point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(3) )[ 1 ] ) );
    }
    else
        csmp_error.Note( csmp::FATAL_ERROR,
                           "CornerPointCell::InitializePoly6Element():", "Unpredictable case for 6 point cell!!");
    CheckFaceOrder();
    CheckNodeOrder();
}



void CornerPointCell::ProcessPoly6Element( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// (1) wedge ( 3 tetrahedrons )
    /// (2) octahedron ( 2 tetrahedrons )

    if( ( meshing_cycle_ < 2 ) && ( cell_type_ == CORNER_POINT_CELL_6_WEDGE ) )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_6_WEDGE ) )
    {
        /// edges and faces are not meshed
        processPrismElement(tetra_mesh,
                            this->grid_, additional_points, additional_edges,
                            this->elements_,this->faces_,this->extra_nodes_,
                            this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
    else if( ( meshing_cycle_ < 3 ) && ( cell_type_ == CORNER_POINT_CELL_6_OCTAHEDRON ) )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
}


/**
    InitializePoly5Element():
    Poly5: sketch

 augmented tetrahedron:


       |              |
       |              |
     4 x              |
       |              |
       |    |         |    |
       |    |         |    |
       x____|_________x    |
     2 |\   |        1|\   |
       | \  |         | \  |
       |  \ |         |  \ |
       |   \x_________|___\x
       |  3 |         |    | 0
            |              |
            |              |
            |              |

 pyramid: ( shouldn't be observed )

                  pillars(0,1,2,3)

                   \\  //
                    \\//
                      x
                    / | \
                   / / \ \
                  /  / \  \
                 /  /   \  \
                /   /    \  \
               /   /     \  \
              /   /       \  \
             /    /       \   \
            /    /         \  \
           x____ /_________x   \
        2 / \   /        1 \\   \
         /   \  /          \ \  \
   pillar(2)  \ /  pillar(1)  \ \
               \x______________\x
              3 /              \ 0
                /               \
     pillar(3)  /               \ pillar(0)

*/

void CornerPointCell
::InitializePoly5Element(  const std::map<csmp::Point<3U>,size_t>&               cell_nodes,
                           const std::vector<std::vector<csmp::Point<3U> > >&    face_ordered_points,
                           const std::vector<uint32_t>&                            tri_faces,
                           const std::vector<uint32_t>&                            quad_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_tri_faces,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_quad_faces,
                           const std::vector<std::vector<csmp::Point<3U> > >&    edge_ordered_points,
                           const std::map<csmp::Point<3U>,std::vector<uint32_t> >& point_edges )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// point sets
    std::set<csmp::Point<3U> > points_set;

    /// node id
    size_t current_id;

    /// hexa cell info
    const size_t quad_fem_nodes( 4U  );

    /// (1) 2 tetrahedrons or tetrahedron + triangle
    /// (2) pyramid ( 2 tetrahedrons ) ( shouldn't be observed )

    /// 1st case: tetrahedron + triangle
    if( ( num_quad_faces_ == 2 ) && ( num_tri_faces_ == 2 ) )
    {
        /// observed and tested
        cell_type_     = CORNER_POINT_CELL_5_AUGMENTED_TETRA;
        cell_category_ = CORNER_POINT_CELL_DIM_3_OVERLAPPING;

        /// renumbering nodes
        /// 0 - side node
        current_id = 0;
        points_set.clear();
        for( size_t nid =0; nid < this->GetNumNodes(); ++nid )
            if( point_edges.at( this->GetPointOriginalOrder( nid ) ).size() == 2 )
            {
                points_set.insert( this->GetPointOriginalOrder( nid ) );
                this->AssignNodeOrder( current_id++, nid );
                break;
            }
        assert( current_id == 1 );
        /// 1,2,3 - bottom quad nodes
        for( size_t nid = 0; nid < 4U; ++nid )
            if( points_set.find( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ nid ] ) != points_set.end() )
            {
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 1)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 2)%quad_fem_nodes ] ) );
                points_set.insert( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] );
                this->AssignNodeOrder( current_id++, cell_nodes.at( face_ordered_points.at( point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] )[ (nid + 3)%quad_fem_nodes ] ) );
                break;
            }
        assert( current_id == 4 );
        /// 4 - top quad node
        for( std::vector<uint32_t>::const_iterator
             eit = point_edges.at( this->GetPointCustomOrder(2) ).begin(); eit != point_edges.at( this->GetPointCustomOrder(2) ).end(); ++eit )
            for( size_t nid = 0; nid < 2U; ++nid )
                if( points_set.find( edge_ordered_points[ *eit ][ nid ] ) == points_set.end() )
                {
                    points_set.insert( edge_ordered_points[ *eit ][ nid ] );
                    this->AssignNodeOrder( current_id++, cell_nodes.at( edge_ordered_points[ *eit ][ nid ] ) );
                    break;
                }
        assert( current_id == 5 );
        assert( points_set.size() == 5 );

        /// reorder faces
        std::set<uint32_t> faces;
        for( uint32_t fid{0U}; fid<6; ++fid )
            faces.insert(fid);
        /// bottom quad face ( 0321 )
        this->AssignFaceOrder( 0, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] );
        faces.erase(this->GetOriginalFaceId(0));
        /// top quad face  ( 0143 )
        this->AssignFaceOrder( 1, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 1 ] );
        faces.erase(this->GetOriginalFaceId(1));
        /// tri face   ( 124 )
        this->AssignFaceOrder( 2, point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] );
        faces.erase(this->GetOriginalFaceId(2));
        /// tri face   ( 234 )
        this->AssignFaceOrder( 3, point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] );
        faces.erase(this->GetOriginalFaceId(3));
        /// edge faces( 01 and 03 )
        size_t fid = 3;
        for(std::set<uint32_t>::const_iterator fit = faces.begin(); fit != faces.end(); ++fit )
            this->AssignFaceOrder( ++fid, *fit );
    }
    // TODO: numbering seems to be broken
    /// 2nd case: pyramid
    else if( ( num_quad_faces_ == 1 ) && ( num_tri_faces_ == 4 ) )
    {
        /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
        cell_type_     = CORNER_POINT_CELL_5_PYRAMID;
        cell_category_ = CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING;

        /// renumbering nodes
        /// 0,1,2,3 - quad face
        current_id = 0;
        for( typename std::vector<csmp::Point<3U> >::const_iterator
             it = face_ordered_points.at( quad_faces[ 0 ] ).begin(); it != face_ordered_points.at( quad_faces[ 0] ).end(); ++it )
        {
            points_set.insert( *it );
            this->AssignNodeOrder( current_id++, cell_nodes.at( *it ) );
        }
        assert( current_id == 4 );
        /// 4 - apex
        for( size_t nid = 0; nid < this->GetNumNodes(); ++nid )
            if( point_edges.at( this->GetPointOriginalOrder( nid ) ).size() == 4 )
            {
                points_set.insert( this->GetPointOriginalOrder( nid ) );
                this->AssignNodeOrder( current_id++, nid );
                break;
            }
        assert( current_id == 5 );
        assert( points_set.size() == 5 );

        /// reorder faces
        std::set<uint32_t> faces;
        for( uint32_t fid{0U}; fid<6; ++fid )
            faces.insert(fid);
        /// bottom quad face ( 0321 )
        this->AssignFaceOrder( 0, point_quad_faces.at( this->GetPointCustomOrder(0) )[ 0 ] );
        faces.erase(this->GetOriginalFaceId(0));
        /// tri face   ( 014 )
        this->AssignFaceOrder( 1, ( point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] == point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] ? point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] : ( point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] == point_tri_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ? point_tri_faces.at( this->GetPointCustomOrder(0) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(0) )[ 1 ] ) ) );
        faces.erase(this->GetOriginalFaceId(1));
        /// tri face   ( 124 )
        this->AssignFaceOrder( 2, ( point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] != this->GetOriginalFaceId(1) ? point_tri_faces.at( this->GetPointCustomOrder(1) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(1) )[ 1 ] ) );
        faces.erase(this->GetOriginalFaceId(2));
        /// tri face   ( 234 )
        this->AssignFaceOrder( 3, ( point_tri_faces.at( this->GetPointCustomOrder(2) )[ 0 ] != this->GetOriginalFaceId(2) ? point_tri_faces.at( this->GetPointCustomOrder(2) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(2) )[ 1 ] ) );
        faces.erase(this->GetOriginalFaceId(3));
        /// tri face   ( 304 )
        this->AssignFaceOrder( 4, ( point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] != this->GetOriginalFaceId(3) ? point_tri_faces.at( this->GetPointCustomOrder(3) )[ 0 ] : point_tri_faces.at( this->GetPointCustomOrder(3) )[ 1 ] ) );
        faces.erase(this->GetOriginalFaceId(4));
        /// point faces ( 4 )
        size_t fid = 4;
        for(std::set<uint32_t>::const_iterator fit = faces.begin(); fit != faces.end(); ++fit )
            this->AssignFaceOrder( ++fid, *fit );
    }
    else
        csmp_error.Note( csmp::FATAL_ERROR,
                           "CornerPointCell::InitializePoly5Element():", "Unpredictable case for 5 point cell!!");
    CheckFaceOrder();
    CheckNodeOrder();
}




void CornerPointCell
::ProcessPoly5Element( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// diagonal nodes
    std::vector<uint32_t> diag_nodes;

    /// (1) 2 tetrahedrons or tetrahedron + triangle
    /// (2) pyramid ( 2 tetrahedrons ) ( shouldn't be observed )

    if( ( meshing_cycle_ < 2 ) && ( cell_type_ == CORNER_POINT_CELL_5_AUGMENTED_TETRA ) )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_5_AUGMENTED_TETRA ) )
    {
        /// edges and faces are not meshed
        processDegenerateHexahedronElement(tetra_mesh,
                              this->grid_, additional_points, additional_edges,
                              this->elements_,this->faces_,this->extra_nodes_,
                              this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
    else if( ( meshing_cycle_ < 2 ) && ( cell_type_ == CORNER_POINT_CELL_5_PYRAMID ) )
    {
        /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_5_PYRAMID ) )
    {
        /// edges and faces are not meshed
        processPyramidElement(tetra_mesh,
                              this->grid_, additional_points, additional_edges,
                              this->elements_,this->faces_,this->extra_nodes_,
                              this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}



/**
    InitializePoly4Element():
    Poly4: sketch

  quadrilateral:

        |            |
        |            |
        |            |
        |            |
      2 x____________x 3
        |\            \
        | \   |        \   |
        |  \  |         \  |
            \ |          \ |
             \|           \|
            1 x____________x 0
              |            |
              |            |
              |            |


  tetrahedron: ( shouldn't be observed )

            \\ pillars (1,4)
             \\
              \\3
                x
                |\
                |\ \
                | \  \
   pillar(0)    | \    \
              \ |  \     \
               \|  \      \
   pillar(2)__1_x___\_______x_2___pillars (2,4)
                 \   \      /\
                  \  \     /  \
                   \  \   /    \
                    \  \ /
                     \ \/
                    0 \x
                       \\
                        \\ pillars (0,1)
*/

void CornerPointCell
::InitializePoly4Element(  const std::vector<uint32_t>& tri_faces,
                           const std::vector<uint32_t>& quad_faces )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// (1) quadrilateral ( 2 triangles )
    /// (2) tetrahedron ( shouldn't be observed )

    /// 1st case: quadrilateral
    if( ( num_quad_faces_ == 1 ) && ( num_tri_faces_ == 0 ) )
    {
        /// observed and tested
        cell_type_     = CORNER_POINT_CELL_4_QUADRILATERAL;
        cell_category_ = CORNER_POINT_CELL_DIM_2;

        size_t current_fid=0;
        for( size_t fid = 0; fid<6; ++fid )
            if( this->GetNumPolygonFaceNodes(fid) == 4 )
                this->AssignFaceOrder( current_fid++, fid );
        assert( current_fid == 2 );
        for( size_t fid = 0; fid<6; ++fid )
            if( this->GetNumPolygonFaceNodes(fid) != 4 )
                this->AssignFaceOrder( current_fid++, fid );
        assert( current_fid == 6 );
        size_t current_nid=0;
        for( size_t nid = 0; nid<4; ++nid )
            this->AssignNodeOrder( current_nid++, this->GetPolygonFaceNodeLocalId( this->GetOriginalFaceId(0), nid ) );
        assert( current_nid == 4 );
    }
    /// 2nd case: tetrahedron
    else if( ( num_quad_faces_ == 0 ) && ( num_tri_faces_ == 4 ) )
    {
        /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
        cell_type_     = CORNER_POINT_CELL_4_TETRAHEDRON;
        cell_category_ = CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING;

        size_t current_fid=0;
        for( size_t fid = 0; fid<6; ++fid )
            if( this->GetNumPolygonFaceNodes(fid) == 3 )
                this->AssignFaceOrder( current_fid++, fid );
        assert( current_fid == 4 );
        for( size_t fid = 0; fid<6; ++fid )
            if( this->GetNumPolygonFaceNodes(fid) != 3 )
                this->AssignFaceOrder( current_fid++, fid );
        assert( current_fid == 6 );
    }
    else
        csmp_error.Note( csmp::FATAL_ERROR,
                           "CornerPointCell::InitializePoly4Element():", "Unpredictable case for 4 point cell!!");
    CheckFaceOrder();
    CheckNodeOrder();
}



void CornerPointCell::ProcessPoly4Element( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// (1) quadrilateral ( 2 triangles )
    /// (2) tetrahedron ( shouldn't be observed )

    if( ( meshing_cycle_ < 2 ) && ( cell_type_ == CORNER_POINT_CELL_4_QUADRILATERAL ) )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_4_QUADRILATERAL ) )
    {
        /// edges and faces are not meshed
        processQuadrilateralElement(tetra_mesh,
                                    this->grid_, additional_points, additional_edges,
                                    this->elements_,this->faces_,this->extra_nodes_,
                                    this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
    /// SHOULDN'T BE OBSERVED ( only if pillars intersect multiple times )
    else if( ( meshing_cycle_ < 2 ) && ( cell_type_ == CORNER_POINT_CELL_4_TETRAHEDRON ) )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_4_TETRAHEDRON ) )
    {
        /// edges and faces are not meshed
        processTetrahedronElement(tetra_mesh,
                                  this->grid_, additional_points, additional_edges,
                                  this->elements_,this->faces_,this->extra_nodes_,
                                  this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}





/**
    InitializePoly3Element():
    Poly3: sketch

  triangle ( shouldn't be observed ):

         |            |
         |            |
         |            |
      0  x____________x 1
         |\          /|
         | \        / |
         |  \      /  |
         |   \    /   |
         |    \  /    |
   pillar(0)   \/   pillar(1)
               x2
              /  \
             /    \
           pillars(2,3)

*/

void CornerPointCell
::InitializePoly3Element( )
{
    /// triangle ( shouldn't be observed)

    cell_type_     = CORNER_POINT_CELL_3_TRIANGLE;
    cell_category_ = CORNER_POINT_CELL_DIM_2;

    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )

    size_t current_fid=0;
    for( size_t fid = 0; fid<6; ++fid )
        if( this->GetNumPolygonFaceNodes(fid) == 3 )
            this->AssignFaceOrder( current_fid++, fid );
    for( size_t fid = 0; fid<6; ++fid )
        if( this->GetNumPolygonFaceNodes(fid) != 3 )
            this->AssignFaceOrder( current_fid++, fid );
    assert( current_fid == 6 );
    CheckFaceOrder();
    CheckNodeOrder();
}




void CornerPointCell
::ProcessPoly3Element( std::map<std::set<uint32_t>,GridNode*>& additional_points, std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
    if( meshing_cycle_ < 2 )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_3_TRIANGLE ) )
    {
        /// edges and faces are not meshed
        processTriangleElement(tetra_mesh,
                               this->grid_, additional_points, additional_edges,
                               this->elements_,this->faces_,this->extra_nodes_,
                               this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}


/**
    InitializePoly2Element():
    Poly2: sketch

  bar ( shouldn't be observed ):

       \   /        \   /
        \ /          \ /
      0  x____________x 1
        / \          / \
       /   \        /   \
   pillars(0,2)   pillars(1,3)

*/

void CornerPointCell::InitializePoly2Element()
{
    /// bar ( shouldn't be observed)
    cell_type_     = CORNER_POINT_CELL_2_BAR;
    cell_category_ = CORNER_POINT_CELL_DIM_1;

    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
    size_t current_fid=0;
    for( size_t fid = 0; fid<6; ++fid )
        if( this->GetNumPolygonFaceNodes(fid) == 2 )
            this->AssignFaceOrder( current_fid++, fid );
    for( size_t fid = 0; fid<6; ++fid )
        if( this->GetNumPolygonFaceNodes(fid) != 2 )
            this->AssignFaceOrder( current_fid++, fid );
    assert( current_fid == 6 );
    CheckFaceOrder();
    CheckNodeOrder();
}



void CornerPointCell::ProcessPoly2Element( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
    if( meshing_cycle_ < 2 )
    {
        ProcessPolyElement(additional_points,additional_edges,well_element);
    }
    else if( ( meshing_cycle_ == 2U )  && ( cell_type_ == CORNER_POINT_CELL_2_BAR ) )
    {
        /// edges and faces are not meshed
        processBarElement(tetra_mesh,
                          this->grid_, additional_points, additional_edges,
                          this->elements_,this->faces_,this->extra_nodes_,
                          this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}

/**
    InitializePoly1Element():
    Poly1: sketch

  point ( shouldn't be observed ):

      pillars(0,1,2,3)

       \\  //
        \\//
          x
        / | \
       / / \ \
      /  / \  \

*/

void CornerPointCell::InitializePoly1Element()
{
    /// points ( shouldn't be observed)
    cell_type_     = CORNER_POINT_CELL_1_POINT;
    cell_category_ = CORNER_POINT_CELL_DIM_0;

    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
}



void CornerPointCell::ProcessPoly1Element( std::map<std::set<uint32_t>,GridNode*>& additional_points,
                                                std::set<std::set<uint32_t> >& additional_edges, bool well_element, bool tetra_mesh )
{
    /// SHOULDN'T BE OBSERVED ( only if pillars intersect )
    if( meshing_cycle_ != -1U )
    {
        processPointElement(tetra_mesh,
                            this->grid_, additional_points, additional_edges,
                            this->elements_,this->faces_,this->extra_nodes_,
                            this->nodes_in_custom_order_,this->face_nodes_in_custom_order_ );
        FinishMeshing();
        return;
    }
}




std::string toString( CORNER_POINT_CELL_TYPE ctype )
  {
      if ( ctype == CORNER_POINT_CELL_UNDEFINED )         return "CORNER_POINT_CELL_UNDEFINED";
      if ( ctype == CORNER_POINT_CELL_8_HEXAHEDRON )      return "CORNER_POINT_CELL_8_HEXAHEDRON";
      if ( ctype == CORNER_POINT_CELL_7_AUGMENTED_PRISM ) return "CORNER_POINT_CELL_7_AUGMENTED_PRISM";
      if ( ctype == CORNER_POINT_CELL_6_WEDGE )           return "CORNER_POINT_CELL_6_WEDGE";
      if ( ctype == CORNER_POINT_CELL_6_OCTAHEDRON  )     return "CORNER_POINT_CELL_5_AUGMENTED_TETRA";
      if ( ctype == CORNER_POINT_CELL_5_AUGMENTED_TETRA ) return "CORNER_POINT_CELL_5_AUGMENTED_TETRA";  /// observed and tested
      if ( ctype == CORNER_POINT_CELL_5_PYRAMID )         return "CORNER_POINT_CELL_5_PYRAMID";          /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
      if ( ctype == CORNER_POINT_CELL_4_TETRAHEDRON )     return "CORNER_POINT_CELL_4_TETRAHEDRON";      /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
      if ( ctype == CORNER_POINT_CELL_4_QUADRILATERAL )   return "CORNER_POINT_CELL_4_QUADRILATERAL";    /// observed and tested
      if ( ctype == CORNER_POINT_CELL_3_TRIANGLE )        return "CORNER_POINT_CELL_3_TRIANGLE";         /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
      if ( ctype == CORNER_POINT_CELL_2_BAR )             return "CORNER_POINT_CELL_2_BAR";              /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
      if ( ctype == CORNER_POINT_CELL_1_POINT )           return "CORNER_POINT_CELL_1_POINT";            /// --- SHOULDN'T BE OBSERVED ( only if pillars intersect )
      return "CORNER_POINT_CELL_UNDEFINED";
  }




std::string toString( CORNER_POINT_CELL_CATEGORY cat )
{
    if ( cat == CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING ) return "CORNER_POINT_CELL_DIM_3_NON_OVERLAPPING";
    if ( cat == CORNER_POINT_CELL_DIM_3_OVERLAPPING ) return "CORNER_POINT_CELL_DIM_3_OVERLAPPING";
    if ( cat == CORNER_POINT_CELL_DIM_2 ) return "CORNER_POINT_CELL_DIM_2";
    if ( cat == CORNER_POINT_CELL_DIM_1  ) return "CORNER_POINT_CELL_DIM_1";
    if ( cat == CORNER_POINT_CELL_DIM_0   ) return "CORNER_POINT_CELL_DIM_0";
    return "UNDEFINED";
}



std::string toString( CORNER_POINT_CELL_FACE_INDEX fidx )
 {
    if ( fidx == CORNER_POINT_CELL_FACE_Xplus ) return "CORNER_POINT_CELL_FACE_Xplus";
    if ( fidx == CORNER_POINT_CELL_FACE_Xminus ) return "CORNER_POINT_CELL_FACE_Xminus";
    if ( fidx == CORNER_POINT_CELL_FACE_Yplus ) return "CORNER_POINT_CELL_FACE_Yplus";
    if ( fidx == CORNER_POINT_CELL_FACE_Yminus ) return "CORNER_POINT_CELL_FACE_Yminus";
    if ( fidx == CORNER_POINT_CELL_FACE_Zplus ) return "CORNER_POINT_CELL_FACE_Zplus";
    if ( fidx == CORNER_POINT_CELL_FACE_Zminus ) return "CORNER_POINT_CELL_FACE_Zminus";
    return "UNDEFINED";
}




/**
    Prints object state to screen.
*/
void CornerPointCell::Out() const
 {
    // printing the base class first
    PolygonCell::Out();
   
    std::cout <<"\nCornerPointCell::Out: type: "<< toString(cell_type_) <<", category: "<< toString(cell_category_);
    std::cout <<"\n\tquadrilateral faces: "<< num_quad_faces_ <<", triangular faces: "<< num_tri_faces_;
    std::cout <<", meshing cycle: "<< meshing_cycle_;
    std::cout <<"\n\n\t pillar nodes (pillar_nodes_.size()):\n";
    for ( auto it=pillar_nodes_.begin(); it!=pillar_nodes_.end(); ++it )
      std::cout << *(*it) <<" ";
    std::cout <<"\n\n\tpillar node order:\n";
    for ( auto it=pillar_nodes_order_.begin(); it!=pillar_nodes_order_.end(); ++it )
      std::cout << (*it) <<" ";
    std::cout <<"\n\n\torigin well face: "<< well_face_org_ <<", destination well face: "<< well_face_dst_;

    std::cout <<"\n\n\twells penetrating the grid cell:\n";
    // vector<pair<csmp::CSMP_FEM_TYPE,vector<GridNode*> > >  wells_;
    for ( auto it=wells_.begin(); it!=wells_.end(); it++ ) {
         std::cout <<"\n"<< parseFiniteElementType((*it).first) <<": ";
         for ( typename std::vector<GridNode*>::const_iterator nit=(*it).second.begin();
              nit!=(*it).second.end(); nit++ )
           std::cout << *(*nit) <<" ";
      }
    std::cout <<"\n\n";
 }

} // end namespace csmp
