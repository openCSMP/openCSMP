#include "CornerPointGrid_UoM.h"
#include "PropertyData.h"
#include "CSMP_highLevelUtilities.h"

#include "ErrorHandler.h"

namespace csmp {

namespace eclipse {

// CORNER POINT GRID

CornerPointGrid::CornerPointGrid()
:NX_(0),
 NY_(0),
 NZ_(0),
 NX_x_NY_(0),
 elements_(0),
 tetra_mesh_(false),
 exclude_inactive_cells_(false)
{
}

CornerPointGrid::CornerPointGrid( const CornerPointGrid& grid )
:NX_             ( grid.NX_ ),
 NY_             ( grid.NY_ ),
 NZ_             ( grid.NZ_ ),
 pillars_        ( grid.pillars_ ),
 embedded_cells_ ( grid.embedded_cells_ ),
 cell_activity_  ( grid.cell_activity_ ),
 faults_data_    ( grid.faults_data_ ),
 well_face_path_ ( grid.well_face_path_ ),
 well_edge_path_ ( grid.well_edge_path_ ),
 axes_           ( grid.axes_ )
{
}

CornerPointGrid& CornerPointGrid::operator=( const CornerPointGrid& grid )
{
    if( this != &grid )
    {
        NX_             = grid.NX_;
        NY_             = grid.NY_;
        NZ_             = grid.NZ_;
        pillars_        = grid.pillars_;
        embedded_cells_ = grid.embedded_cells_;
        cell_activity_  = grid.cell_activity_;
        faults_data_    = grid.faults_data_;
        well_face_path_ = grid.well_face_path_;
        well_edge_path_ = grid.well_edge_path_;
        axes_           = grid.axes_;
    }
    return *this;
}

CornerPointGrid::~CornerPointGrid()
{

}

void CornerPointGrid::Clear()
{
    /// grid
    embedded_cells_.clear();
    cell_activity_.clear();
    for( std::vector<std::vector<Pillar> >::iterator
         it = pillars_.begin(); it != pillars_.end(); ++it )
        (*it).clear();
    pillars_.clear();
    axes_.clear();
    for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::iterator
        it = faults_data_.begin(); it != faults_data_.end(); ++it )
        (*it).second.clear();
    faults_data_.clear();
    for( std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >::iterator
        it = well_face_path_.begin(); it != well_face_path_.end(); ++it )
        (*it).second.clear();
    well_face_path_.clear();
    for( std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >::iterator
        it = well_edge_path_.begin(); it != well_edge_path_.end(); ++it )
        (*it).second.clear();
    well_edge_path_.clear();
}

void CornerPointGrid::TetraMesh( bool flag )
{
    tetra_mesh_ = flag;
}

void CornerPointGrid::ExcludeInactiveCells( bool flag )
{
    exclude_inactive_cells_ = flag;
}

size_t CornerPointGrid::GetNumCells( ) const
{
    return NX_*NY_*NZ_;
}

size_t CornerPointGrid::CellIndex( size_t i, size_t j, size_t k ) const
{
    return i+j*NX_+k*NX_x_NY_;
}

size_t CornerPointGrid::CellIndexI( size_t cell_id ) const
{
    return (cell_id%NX_);
}

size_t CornerPointGrid::CellIndexJ( size_t cell_id ) const
{
    return ((cell_id/NX_)%NY_);
}

size_t CornerPointGrid::CellIndexK( size_t cell_id ) const
{
    return (cell_id/NX_x_NY_);
}

size_t CornerPointGrid::NodeIndexI( size_t cell_id, size_t nid ) const
{
    return 2*CellIndexI( cell_id ) + NodeIndexIncrementI( nid );
}

size_t CornerPointGrid::NodeIndexJ( size_t cell_id, size_t nid ) const
{
    return 2*CellIndexI( cell_id ) + NodeIndexIncrementJ( nid );
}

size_t CornerPointGrid::NodeIndexK( size_t cell_id, size_t nid ) const
{
    return 2*CellIndexI( cell_id ) + NodeIndexIncrementK( nid );
}

size_t CornerPointGrid::NodeIndexIncrementI( size_t nid ) const
{
    return ( nid%2 );
}

size_t CornerPointGrid::NodeIndexIncrementJ( size_t nid ) const
{
    return ( (nid/2)%2 );
}

size_t CornerPointGrid::NodeIndexIncrementK( size_t nid ) const
{
    return ( nid/4 );
}

void CornerPointGrid::AddWellFacePath( const std::string& well_name, const std::vector<size_t>& cell_ids )
{
    addWellPath( NX_, NY_, NZ_, well_name, cell_ids, well_face_path_ );
}

void CornerPointGrid::AddWellFacePath( const std::string& well_name,
                                            const std::vector<size_t>& cell_ids,
                                            const std::vector<std::pair<size_t,size_t> >& face_ids )
{
    addWellPath( well_name, cell_ids, face_ids, well_face_path_ );
}

void CornerPointGrid::AddWellEdgePath( const std::string& well_name,
                                            const std::vector<size_t>& cell_ids,
                                            const std::vector<std::pair<size_t,size_t> >& edge_ids )
{
    addWellPath( well_name, cell_ids, edge_ids, well_edge_path_ );
}

std::vector<std::vector<Pillar> >& CornerPointGrid::GetPillars()
{
    return pillars_;
}

std::vector<size_t>& CornerPointGrid::GetCellActivity()
{
    return cell_activity_;
}

std::map<std::string,std::vector<std::pair<size_t,size_t> > >& CornerPointGrid::GetFaultData()
{
    return faults_data_;
}

std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& CornerPointGrid::GetWellFacePath()
{
    return well_face_path_;
}

std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& CornerPointGrid::GetWellEdgePath()
{
    return well_edge_path_;
}

void CornerPointGrid::AssignDimensionX( size_t NX )
{
    NX_ = NX;
}

void CornerPointGrid::AssignDimensionY( size_t NY )
{
    NY_ = NY;
}

void CornerPointGrid::AssignDimensionZ( size_t NZ )
{
    NZ_ = NZ;
}


void CornerPointGrid::MinMaxCoordinates( PolygonGridManager& pgm, Point<3U>& xyz_min, Point<3U>& xyz_max ) const
 {
    xyz_min = pgm.GetPoint(0);
    xyz_max = pgm.GetPoint(0);
    const size_t num_nodes( pgm.GetNumNodes() );
    for ( size_t nid = 1U; nid < num_nodes; ++nid )
    {
        const csmp::Point<3U>& p = pgm.GetPoint( nid );
        for( size_t i = 0; i < 3U; i++ )
        {
            xyz_min[i] = std::min( xyz_min[i], p[i] );
            xyz_max[i] = std::min( xyz_max[i], p[i] );
        }
    }
}




void CornerPointGrid::DefineAxes()
{
    csmp::Point<3U> i(0.0);
    i[0] = 1.0;
    csmp::Point<3U> j(0.0);
    j[1] = 1.0;
    csmp::Point<3U> k(0.0);
    k[2] = 1.0;
    axes_.clear();
    axes_.push_back( csmp::Point<3U>(0.0) ); // origin point
    axes_.push_back( i );
    axes_.push_back( j );
    axes_.push_back( k );
}



/**
  x stays, y = -z and z=y.
*/
void CornerPointGrid
::ConvertFromReservoirToCSMPcoordinateSystem( csmp::Point<3U>& pt )
{
    const double y( pt[1U] );
    pt[0U] = pt[0U];
    pt[1U] = -pt[2U];
    pt[2U] = y;
}


const csmp::Point<3U>& CornerPointGrid
::GetCellNode( size_t i, size_t j, size_t k, size_t nid ) const
{
    const size_t node_incrementJ(NodeIndexIncrementJ( nid ));
    const size_t node_incrementI(NodeIndexIncrementI( nid ));
    return pillars_[ j + node_incrementJ ][ i + node_incrementI ].GetPoint( i, j, k, nid );
}



void CornerPointGrid
::AssignCellNodeToPillar( size_t i, size_t j, size_t k, size_t nid, const csmp::Point<3U>& pt )
{
    const size_t node_incrementJ(NodeIndexIncrementJ( nid ));
    const size_t node_incrementI(NodeIndexIncrementI( nid ));
    return pillars_[ j + node_incrementJ ][ i + node_incrementI ].AssignPoint( i, j, k, nid, pt );
}


void CornerPointGrid
::InitializeGridSpecs( )
{
    NX_x_NY_  = NX_*NY_;
    elements_ = NX_x_NY_*NZ_;
}





/**
    Attaches the cells to the pillars.
    Repairs cells where degenerate cells arise because beads are missing; situations considered:
 
    1) top nodes of cell are lower that bottom nodes of same cell ( top point = btm point )
    2) bottom nodes of cell are lower than top nodes of neighbouring cell ( top point = btm point )
*/
void CornerPointGrid::ProcessPillars()
{
    // Assign attached cells to pillars

    // internal pillars
    for( size_t i=1; i<NX_; ++i )
        for( size_t j=1; j<NY_; ++j )
        {
            pillars_[j][i].AssignAttachedCellId( i-1, j-1, 0 );
            pillars_[j][i].AssignAttachedCellId( i,   j-1, 1 );
            pillars_[j][i].AssignAttachedCellId( i-1, j,   2 );
            pillars_[j][i].AssignAttachedCellId( i,   j,   3 );
        }
    /// lateral pillars
    for( size_t j=1; j<NY_; ++j )
      {
          pillars_[j][0].AssignAttachedCellId(   0,     j-1, 0 );
          pillars_[j][0].AssignAttachedCellId(   0,     j,   1 );
          pillars_[j][NX_].AssignAttachedCellId( NX_-1, j-1, 0 );
          pillars_[j][NX_].AssignAttachedCellId( NX_-1, j,   1 );
      }
    for( size_t i=1; i<NX_; ++i )
      {
          pillars_[0][i].AssignAttachedCellId(   i-1,   0,     0 );
          pillars_[0][i].AssignAttachedCellId(   i,     0,     1 );
          pillars_[NY_][i].AssignAttachedCellId( i-1,   NY_-1, 0 );
          pillars_[NY_][i].AssignAttachedCellId( i,     NY_-1, 1 );
      }
    /// corner pillars
    pillars_[0][0].AssignAttachedCellId(       0,     0,     0 );
    pillars_[0][NX_].AssignAttachedCellId(     NX_-1, 0,     0 );
    pillars_[NY_][0].AssignAttachedCellId(     0,     NY_-1, 0 );
    pillars_[NY_][NX_].AssignAttachedCellId(   NX_-1, NY_-1, 0 );

// TESTING - in some models, no points are stored on the pillars
// for( size_t i=1; i<NX_; ++i )
//  for( size_t j=1; j<NY_; ++j ) pillars_[i][j].Out();

    /// repair pillar nodes in the following cases:
    /// 1) top nodes of cell are lower that bottom nodes of same cell ( top point = btm point )
    /// 2) bottom nodes of cell are lower than top nodes of neighbouring cell ( top point = btm point )
    const size_t quad_fem_nodes( 4U );
    /// process top nodes below bottom nodes
    bool first_print(true);
    for( size_t k = 0; k < NZ_; ++k )
        for( size_t j = 0; j < NY_; ++j )
            for( size_t i = 0; i < NX_; ++i )
            {
                if ( pillars_[i][j].GetNumPoints() == 0 ) {
                     if ( first_print ) {
// SKM testing
std::cerr <<"\nCornerPointGrid::ProcessPillars: pillars without beads: ";
                          first_print = false;
                       }
// SKM testing
std::cerr <<"pillar("<< i <<","<< j <<") ";
                     continue;
                  }
                /// hexa cell global id

                /// indexing cell nodes
                for( size_t nid = 0; nid < quad_fem_nodes; ++nid ) {
                    /// read top node on current pillar
                    csmp::Point<3U>  pt_top = GetCellNode(i,j,k,nid);
                    /// read bottom node on current pillar
                    csmp::Point<3U>  pt_btm = GetCellNode(i,j,k,nid+4);
                    if( pt_top[3] > pt_btm[3] )
                    {
                        pt_top = pt_btm;
                        AssignCellNodeToPillar(i,j,k,nid,pt_top);
                    }
                }
            }
    /// process bottom nodes below top nodes of neighbour
    bool modified( false );
    for( size_t k = 0; k < (NZ_-1); ++k )
        for( size_t j = 0; j < NY_; ++j )
            for( size_t i = 0; i < NX_; ++i )
            {
                /// hexa cell global id

                /// indexing cell nodes
                for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                {
                    /// read bottom node on current cell
                    csmp::Point<3U>  pt_btm = GetCellNode(i,j,k,nid+4);
                    /// read top node on neighbouring cell
                    csmp::Point<3U>  pt_top = GetCellNode(i,j,k+1,nid);
                    if( pt_btm[3] > pt_top[3] )
                    {
                        pt_top = pt_btm;
                        AssignCellNodeToPillar(i,j,k+1,nid,pt_top);
                        modified = true;
                    }
                }
            }
    /// process top nodes below bottom nodes ( if some nodes were moved )
    if( modified )
    {
        for( size_t k = 0; k < NZ_; ++k )
            for( size_t j = 0; j < NY_; ++j )
                for( size_t i = 0; i < NX_; ++i )
                {
                    /// hexa cell global id

                    /// indexing cell nodes
                    for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                    {
                        /// read top node on current pillar
                        csmp::Point<3U>  pt_top;
                        pt_top = GetCellNode(i,j,k,nid);
                        /// read bottom node on current pillar
                        csmp::Point<3U>  pt_btm;
                        pt_btm = GetCellNode(i,j,k,nid+4);
                        if( pt_top[3] > pt_btm[3] )
                        {
                            pt_top = pt_btm;
                            AssignCellNodeToPillar(i,j,k,nid,pt_top);
                        }
                    }
                }
    }
  
} // end ProcessPillars






/**
    MASTER METHOD for the creation of corner point grids 
    from the original cell-centered grids read by the EclipseInterface.
    
    The PolygonGridManager is used to construct the pillars.
    
    @todo this is confused! - why should the corner point grid be responsible for building the CSMP model? - code should be in the EclipseModel
    
    @todo regions are not recognised properly.
*/
void CornerPointGrid::CreateModel( const std::string&     model_name,
                                   csmp::VSet<3U>&       vset,
                                   csmp::ModelTopology&   model_topology,
                                   std::set<std::string>& regions,
                                   std::set<std::string>& faults,
                                   std::set<std::string>& wells )
 {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    PolygonGridManager pgm;

    // ==============================================================
    // 0. Initialise grid specs and map between pillar and cell nodes
    // ==============================================================
    InitializeGridSpecs();
    ProcessPillars();

    // fe types and cell ids
    csmp::CSMP_FEM_TYPE   hexa_fem_type( csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON );
    std::string           hexa_fem_typename( csmp::parseFiniteElementType( hexa_fem_type ) );
    std::set<std::string> cell_fem_types;
    std::set<std::string> face_fem_types;
    std::set<std::string> edge_fem_types;

    // cell ids
    std::vector<size_t>   cell_elmts;
    std::vector<size_t>   face_elmts;
    std::vector<size_t>   edge_elmts;

    const size_t hexa_fem_nodes( 8U );
    size_t hexa_cell_id( 0 );
    size_t cell_id ( 0 );
    size_t num_cell_nodes;
    std::vector<CornerPointCell> poly( elements_, CornerPointCell(pgm) );
    std::set<size_t> regular_cells;
    std::set<size_t> regular_well_cells;
    std::set<size_t> degenerate_well_cells;
    std::set<size_t> degenerate_nonoverlap_cells;
    std::set<size_t> degenerate_overlap_cells;
    std::set<size_t> degenerate_lowdim_cells;
    std::set<std::set<size_t> > additional_edges;
    std::map<std::set<size_t>,GridNode*> additional_points;
    std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > embedded_cells_vector;

    // if there is no cell activity information assume that all cells are active
    if( cell_activity_.empty() )
        cell_activity_.resize( elements_, 1 );


    // ================================
    // 1. Initialize polyhedral cells
    // ================================
    // Cell index (i,j,k) = i +j*nx + k*nx*ny
    // convert from reservoir coordinate system (x,y,z) to csmp coordinate system (x,-z,y)
    // renumber nodes and eliminate duplicates )
   
    for( size_t k = 0; k < NZ_; ++k )
        for( size_t j = 0; j < NY_; ++j )
            for( size_t i = 0; i < NX_; ++i )
            {
                /// hexa cell global id
                hexa_cell_id = CellIndex(i,j,k);

                if( cell_activity_[ hexa_cell_id ] == 1 || !exclude_inactive_cells_ )
                {
                    /// indexing cell nodes
                    for( size_t nid = 0; nid < hexa_fem_nodes; ++nid )
                    {
                        /// convert from reservoir to csmp coordinate system
                        csmp::Point<3U>  pt = GetCellNode(i,j,k,nid);
                        ConvertFromReservoirToCSMPcoordinateSystem(pt);
                        AssignCellNodeToPillar(i,j,k,nid,pt);
                        poly[ hexa_cell_id ].AssignNode( nid, pt );
                    }

                    /// initialize corner point cell
                    poly[ hexa_cell_id ].InitializeCornerPointCell();
                    num_cell_nodes = poly[ hexa_cell_id ].GetNumNodes();

                    /// add cell to the list of degenerates if it the case
                    if( num_cell_nodes < hexa_fem_nodes )
                    {
                        if( poly[ hexa_cell_id ].IsVolumetricNonOverlapingElement() )
                            degenerate_nonoverlap_cells.insert( hexa_cell_id );
                        else if( poly[ hexa_cell_id ].IsVolumetricOverlapingElement() )
                            degenerate_overlap_cells.insert( hexa_cell_id );
                        else if( poly[ hexa_cell_id ].IsLowDimensionalElement() )
                            degenerate_lowdim_cells.insert( hexa_cell_id );
                        /// resolve mesh conflicts
                        poly[ hexa_cell_id ].ResolveCornerPointCellConflicts( additional_points, additional_edges );
                    }
                    else
                        regular_cells.insert( hexa_cell_id );
                }
            }

    // processing cells with wells
    size_t face_id_org,face_id_dst;
    for( std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >::const_iterator
         it  = well_face_path_.begin(); it != well_face_path_.end(); ++it )
        for( std::vector<std::pair<size_t,std::pair<size_t,size_t> > >::const_iterator
             wit  = (*it).second.begin(); wit != (*it).second.end(); ++wit )
        {
            /// hexa cell global id
            hexa_cell_id = (*wit).first;
            face_id_org  = (*wit).second.first;
            face_id_dst  = (*wit).second.second;
            if( cell_activity_[ hexa_cell_id ] == 1 || !exclude_inactive_cells_ )
            {
                num_cell_nodes = poly[ hexa_cell_id ].GetNumNodes();

                poly[ hexa_cell_id ].AddCornerPointWellCell( face_id_org, face_id_dst );
                /// add cell to the list of degenerates if it the case
                if( num_cell_nodes < hexa_fem_nodes )
                    degenerate_well_cells.insert( hexa_cell_id );
                else
                    regular_well_cells.insert( hexa_cell_id );
            }
        }

    DefineAxes();

    // ================================
    // 2. Assign cell data to VSet
    // ================================
   
    // 2.1: processing non degenerate cells only meaning that there will be no tetrahedra
    // ==================================================================================
   
    if ( degenerate_nonoverlap_cells.empty() &&
         degenerate_overlap_cells.empty() &&
         degenerate_lowdim_cells.empty() &&
         regular_well_cells.empty() &&
         degenerate_well_cells.empty() &&
         !tetra_mesh_ )
      {
          /// cell types
          cell_fem_types.insert( hexa_fem_typename );

          /// assign element data to VSet
          cell_id = 0;
          for( size_t k = 0; k < NZ_; ++k )
              for( size_t j = 0; j < NY_; ++j )
                  for( size_t i = 0; i < NX_; ++i )
                  {
                      /// hexa cell global id
                      hexa_cell_id = CellIndex(i,j,k);

                      if( cell_activity_[ hexa_cell_id ] == 1 || !exclude_inactive_cells_ )
                      {
                          /// fill VSet
                          vset.ResizePlist( cell_id + 1, hexa_fem_nodes );
                          vset.ResizeElementTypes( cell_id + 1 );
                          for( size_t nid = 0U; nid < hexa_fem_nodes; ++nid )
                              vset.Plist( cell_id, nid, poly[ hexa_cell_id ].GetNodeGlobalIdCustomOrder( nid ) );
                          vset.ElementType( cell_id, hexa_fem_type );
                          cell_elmts.push_back( cell_id );

                          embedded_cells_vector.clear();
                          embedded_cells_vector.push_back( std::make_pair( cell_id, hexa_fem_type ) );
                          ++cell_id;

                          embedded_cells_.insert( std::make_pair( hexa_cell_id, embedded_cells_vector ) );
                      }
                  }
      }

    // 2.2: processing degenerate cells
    // ===============================================================================
    // which appear due to collapse of point pairs along pillars.
    // main reason: fault alignment meshing
    // geological reasons:
    // erosion, pinch-out of geological layers,
    // truncated artificial layers ( zero thickness cells )
    else
    {
        // 2.2.1 process well elements
        bool well_element( true );

        // process degenerate cells which contains wells
        for( std::set<size_t>::const_iterator
             it = degenerate_well_cells.begin(); it != degenerate_well_cells.end(); ++it )
          {
              /// hexa cell global id
              hexa_cell_id = *it;
              poly[ hexa_cell_id ].ProcessCornerPointDegenerateCell( additional_points, additional_edges, well_element, tetra_mesh_ );
          }

        // process regular cells which contains wells
        for( std::set<size_t>::const_iterator
             it = regular_well_cells.begin(); it != regular_well_cells.end(); ++it )
          {
              /// hexa cell global id
              hexa_cell_id = *it;
              poly[ hexa_cell_id ].ProcessCornerPointRegularCell( additional_points, additional_edges, well_element, tetra_mesh_ );
          }

        // 2.2.2 process degenerate volumetric elements
        bool well_free_element( false );

        // process degenerate volumetric overlapping elements
        const size_t num_deg_op_cycles( 3U );
        for( size_t meshing_cycle = 0; meshing_cycle < num_deg_op_cycles; ++meshing_cycle )
            for( std::set<size_t>::const_iterator
                 it = degenerate_overlap_cells.begin(); it != degenerate_overlap_cells.end(); ++it )
            {
                /// hexa cell global id
                hexa_cell_id = *it;
                if( ( degenerate_well_cells.find( hexa_cell_id ) == degenerate_well_cells.end() ) && ( poly[ hexa_cell_id ].GetNumNodes() > 4 ) )
                    poly[ hexa_cell_id ].ProcessCornerPointDegenerateVolumetricOverlappingCell( additional_points, additional_edges, well_free_element, tetra_mesh_ );
             }

        // process degenerate volumetric non-overlapping elements
        const size_t num_deg_nop_cycles( 3U );
        for( size_t meshing_cycle = 0; meshing_cycle < num_deg_nop_cycles; ++meshing_cycle )
            for( std::set<size_t>::const_iterator
                 it = degenerate_nonoverlap_cells.begin(); it != degenerate_nonoverlap_cells.end(); ++it )
            {
                /// hexa cell global id
                hexa_cell_id = *it;
                if( ( degenerate_well_cells.find( hexa_cell_id ) == degenerate_well_cells.end() ) && ( poly[ hexa_cell_id ].GetNumNodes() > 4 ) )
                    poly[ hexa_cell_id ].ProcessCornerPointDegenerateVolumetricNonOverlappingCell( additional_points, additional_edges, well_free_element, tetra_mesh_ );
             }

        // 2.2.3 process regular ( hexahedron ) volumetric elements

        // process regular cells
        const size_t num_reg_cycles( 3U );
        for( size_t meshing_cycle = 0; meshing_cycle < num_reg_cycles; ++meshing_cycle )
            for( std::set<size_t>::const_iterator
                 it = regular_cells.begin(); it != regular_cells.end(); ++it )
            {
                /// hexa cell global id
                hexa_cell_id = *it;
                if( regular_well_cells.find( hexa_cell_id ) == regular_well_cells.end() )
                    poly[ hexa_cell_id ].ProcessCornerPointRegularCell( additional_points, additional_edges, well_free_element, tetra_mesh_ );
            }

        // 2.2.4 process degenerate low dimensional elements

        // process degenerate low dimensional elements
        const size_t num_deg_ld_cycles( 3U );
        for( size_t meshing_cycle = 0; meshing_cycle < num_deg_ld_cycles; ++meshing_cycle )
            for( std::set<size_t>::const_iterator
                 it = degenerate_lowdim_cells.begin(); it != degenerate_lowdim_cells.end(); ++it )
            {
                /// hexa cell global id
                hexa_cell_id = *it;
                if( ( degenerate_well_cells.find( hexa_cell_id ) == degenerate_well_cells.end() ) && ( poly[ hexa_cell_id ].GetNumNodes() < 5 ) )
                    poly[ hexa_cell_id ].ProcessCornerPointDegenerateLowDimensionalCell( additional_points, additional_edges, well_free_element, tetra_mesh_ );
             }

        // 2.2.5 write mesh to VSet

        // assign element data to VSet
        cell_id = 0;
        for( size_t k = 0; k < NZ_; ++k )
            for( size_t j = 0; j < NY_; ++j )
                for( size_t i = 0; i < NX_; ++i )
                  {
                     /// hexa cell global id
                     hexa_cell_id = CellIndex(i,j,k);
//std::cerr <<"\nCell ID: "<< hexa_cell_id;
                     if ( cell_activity_[ hexa_cell_id ] == 1 || !exclude_inactive_cells_ )
                     {
//poly[ hexa_cell_id ].Out();

                         /// assign subcells
                         const size_t num_subcells( poly[ hexa_cell_id ].GetNumElements() );
                         embedded_cells_vector.clear();
                         for( size_t eid = 0; eid < num_subcells; ++eid )
                         {
//std::cerr <<" element type "<< poly[ hexa_cell_id ].GetElementType( eid ); // 42=hex, 36=quad, 28=line

                             const size_t num_cell_nodes( poly[ hexa_cell_id ].GetNumElementNodes( eid ) );
                             if( num_cell_nodes > 1 )
                             {
                                 vset.ResizePlist( cell_id + 1, CSMP_ElementSpecifications::NodesPerElementOfType( poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                 vset.ResizeElementTypes( cell_id + 1 );
                                 for( size_t nid = 0U; nid < num_cell_nodes; ++nid )
                                     vset.Plist( cell_id, nid, poly[ hexa_cell_id ].GetElementNodeGlobalId( eid, nid ) );
                                 vset.ElementType( cell_id, poly[ hexa_cell_id ].GetElementType( eid ) );
                                 const size_t cell_dim ( poly[ hexa_cell_id ].GetElementDim( eid ) );
                                 if( cell_dim == 3U )
                                 {
                                     cell_elmts.push_back( cell_id );
                                     cell_fem_types.insert( csmp::parseFiniteElementType( poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     embedded_cells_vector.push_back( std::make_pair( cell_id, poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     ++cell_id;
                                 }
                                 else if( cell_dim == ( 2U ) )
                                 {
                                     face_elmts.push_back( cell_id );
                                     face_fem_types.insert( csmp::parseFiniteElementType( poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     embedded_cells_vector.push_back( std::make_pair( cell_id, poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     ++cell_id;
                                 }
                                 else if( cell_dim == ( 1U ) )
                                 {
                                     edge_elmts.push_back( cell_id );
                                     edge_fem_types.insert( csmp::parseFiniteElementType( poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     embedded_cells_vector.push_back( std::make_pair( cell_id, poly[ hexa_cell_id ].GetElementType( eid ) ) );
                                     ++cell_id;
                                 }
                             }
                         }
                         embedded_cells_.insert( std::make_pair( hexa_cell_id, embedded_cells_vector ) );
                     }
                }
    }
   
// DEBUGGING CODE ONLY
#ifndef NDEBUG
      {
          std::map<size_t,std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > >::const_iterator cit;

          std::vector<size_t> deg_nonoverlap_cells;
          std::set<std::string> deg_nonoverlap_cells_fem_types;
          for( std::set<size_t>::const_iterator
               it = degenerate_nonoverlap_cells.begin(); it != degenerate_nonoverlap_cells.end(); ++it )
          {
              /// hexa cell global id
              hexa_cell_id = *it;
              cit = embedded_cells_.find( hexa_cell_id );
              if( cit != embedded_cells_.end() )
              {
                  for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                       eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                  {
                      if( CSMP_ElementSpecifications::VolumeElement( (*eit).second ) )
                      {
                          deg_nonoverlap_cells.push_back( (*eit).first );
                          deg_nonoverlap_cells_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                      }
                  }
              }
          }
          if( !deg_nonoverlap_cells.empty() )
          {
              model_topology.AddRegion( "DEGENERATE_NONOVERLAP", deg_nonoverlap_cells_fem_types, deg_nonoverlap_cells );
              regions.insert("DEGENERATE_NONOVERLAP");
          }
          std::vector<size_t> deg_overlap_cells;
          std::set<std::string> deg_overlap_cells_fem_types;
          for( std::set<size_t>::const_iterator
               it = degenerate_overlap_cells.begin(); it != degenerate_overlap_cells.end(); ++it )
          {
              /// hexa cell global id
              hexa_cell_id = *it;
              cit = embedded_cells_.find( hexa_cell_id );
              if( cit != embedded_cells_.end() )
              {
                  for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                       eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                  {
                      if( CSMP_ElementSpecifications::VolumeElement( (*eit).second ) )
                      {
                          deg_overlap_cells.push_back( (*eit).first );
                          deg_overlap_cells_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                      }
                  }
              }
          }
          if( !deg_overlap_cells.empty() )
          {
              model_topology.AddRegion( "DEGENERATE_OVERLAP", deg_overlap_cells_fem_types, deg_overlap_cells );
              regions.insert("DEGENERATE_OVERLAP");
          }
      }
#endif

    regular_cells.clear();
    regular_well_cells.clear();
    degenerate_well_cells.clear();
    degenerate_nonoverlap_cells.clear();
    degenerate_overlap_cells.clear();
    degenerate_lowdim_cells.clear();

// DEBUGGING OUTPUT
//std::cerr <<"\nPolygonGridManager: printing the generated grid prior to output to VSet:\n";
//   pgm.Out();
//    std::cerr <<"\nCornerPointGrid::CreateModel: printing the generated cells prior to output to VSet:\n";
//    for ( auto it=poly.begin(); it!=poly.end(); ++it ) {
//         (*it).Out();
//      }

    // 3. Assign node coordinates
    // ===============================================================================
    const size_t num_nodes( pgm.GetNumNodes() );
    vset.ResizeNodes( num_nodes );

    for ( size_t nid=0; nid<num_nodes; ++nid ){
        const csmp::Point<3U>& pt( pgm.GetPoint( nid ) );
        vset.Px( nid, pt[0U] );
        vset.Py( nid, pt[1U] );
        vset.Pz( nid, pt[2U] );
    }
   
/* NA
    else if ( dim == 2U ) {
          for ( size_t nid=0; nid<num_nodes; ++nid ){
              const csmp::Point<3U>& pt1( pgm.GetPoint( nid ) );
              csmp::Point<3U> pt2;
              getCoordinate( axes_[0], axes_[1], axes_[2], pt1, pt2 );
              vset.Px( nid, pt2[0U] );
              vset.Py( nid, pt2[1U] );
              vset.Pz( nid, 0.0 );
          }
      }
    else if ( dim == 1U ) {
          for ( size_t nid=0; nid<num_nodes; ++nid ){
              const csmp::Point<3U>& pt1( pgm.GetPoint( nid ) );
              csmp::Point<3U> pt2;
              getCoordinate( axes_[0], axes_[1], pt1, pt2 );
              vset.Px( nid, pt2[0U] );
              vset.Py( nid, 0.0 );
              vset.Pz( nid, 0.0 );
          }
      }
*/
    // 4. Establishing Model Topology
    // ================================

    // Assign Model name to model topology
    model_topology.ModelName( model_name.c_str() );

    // Process ACTIVE and INACTIVE regions
    EstablishActiveDomain( vset,model_topology,regions,cell_elmts,cell_fem_types );

    // Process node boundary flags
    EstablishBoundaries( vset,model_topology,poly );

    cell_elmts.clear();
    face_elmts.clear();
    edge_elmts.clear();

    EstablishFaultRegions( vset,model_topology,faults,poly );

    EstablishWellRegions( vset,model_topology,wells,poly );

    poly.clear();
    pgm.Clear();

    // assign mesh type (multi-element or single-element type)
    if ( ( cell_fem_types.size() > 1 ) || ( !faults_data_.empty() ))
        vset.HybridElementTypeMesh( true );
    else
        vset.HybridElementTypeMesh( false );
   
   
    // 5. Reporting
    // ================================

    std::cout << "\nCornerPointGrid::CreateModel: VSet resulting from '"<< model_name << "' grid contains: ";
    if( !embedded_cells_.empty() )
      {
          std::map<csmp::CSMP_FEM_TYPE,size_t> mesh_specs;
          std::pair<std::map<csmp::CSMP_FEM_TYPE,size_t>::iterator,bool> nmit;
          std::map<size_t,std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > >::const_iterator cit;
          for( size_t eid = 0; eid < elements_; ++eid )
          {
              cit = embedded_cells_.find( eid );
              if( cit != embedded_cells_.end() )
                  for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                       eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                  {
                      nmit = mesh_specs.insert( std::make_pair( (*eit).second, 0 ) );
                      if( !nmit.second )
                          ++( (*nmit.first).second);
                  }
              else
              {
                  nmit = mesh_specs.insert( std::make_pair( hexa_fem_type, 0 ) );
                  if( !nmit.second )
                      ++( (*nmit.first).second);
              }
          }
          for( std::map<csmp::CSMP_FEM_TYPE,size_t>::const_iterator
               mit = mesh_specs.begin(); mit != mesh_specs.end(); ++ mit )
              std::cout<<"\n" << csmp::parseFiniteElementType( (*mit).first )<< " elements: "<< (*mit).second;
      }
    else std::cout<<"\n" << hexa_fem_typename << " elements: "<< elements_;

    std::cout<<"\n\tTotal number of nodes: "    << vset.Vertices();
    std::cout<<"\n\tTotal number of elements: " << vset.Elements();

    // topology info
    const size_t num_regions( model_topology.ModelRegions() );
    std::cout<<"\n\tTotal number of regions: "    << num_regions;
    std::cout<<"\n\tList of regions: ";
    const size_t max_per_line( 5 );
    std::list<std::string> regions_list;
    model_topology.Out( regions_list );
    size_t rid( 0 );
    for( std::list<std::string>::const_iterator
         rit = regions_list.begin(); rit != regions_list.end(); ++rit )
      {
          std::cout<<"'" << (*rit) << "'";
          ++rid;
          if( ( rid%max_per_line == 0 ) && ( rid != num_regions ) )
              std::cout<< "," << std::endl;
          else if( rid != num_regions )
              std::cout << ", ";
      }
    std::cout<< std::endl;

    if( csmp_error.Verbose() )
        std::cout <<"\nCornerPointGrid::CreateModel: mesh was successfully output to VSet." << std::endl;

} // CreateModel (MASTER - MONSTER METHOD)




/**
    Eliminating inactive cells in VSet.
    
    @todo SKM question: are inactive cells included in the VSet? - they should be eliminated before because they 
    often dominate the entire cell count.
*/
void CornerPointGrid::EstablishActiveDomain( csmp::VSet<3U>& vset,
                                             csmp::ModelTopology& model_topology,
                                             std::set<std::string>& regions,
                                             const std::vector<size_t>& cell_elmts,
                                             const std::set<std::string>& cell_fem_types )
 {
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// active set of cells
    std::set<std::string> active_cell_fem_types;
    std::set<std::string> active_face_fem_types;
    std::set<std::string> active_edge_fem_types;
    std::vector<size_t>   active_cells;
    std::vector<size_t>   active_faces;
    std::vector<size_t>   active_edges;
    /// inactive set of cells
    std::set<std::string> inactive_cell_fem_types;
    std::set<std::string> inactive_face_fem_types;
    std::set<std::string> inactive_edge_fem_types;
    std::vector<size_t>   inactive_cells;
    std::vector<size_t>   inactive_faces;
    std::vector<size_t>   inactive_edge;


    /// correcting active cell id's due to existance of embedded cells
    if( !embedded_cells_.empty() )
    {
        std::map<size_t,std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > >::const_iterator cit;
        const size_t num_elements( cell_activity_.size() );
        for( size_t cid = 0; cid<num_elements; ++cid )
        {
            cit = embedded_cells_.find( cid );
            if( cit != embedded_cells_.end() )
            {
                if( cell_activity_[ cid ] == 1 )
                {
                    for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                         eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                    {
                        if( CSMP_ElementSpecifications::VolumeElement( (*eit).second ) )
                        {
                            active_cells.push_back( (*eit).first );
                            active_cell_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                        else if( CSMP_ElementSpecifications::SurfaceElement( (*eit).second ) )
                        {
                            active_faces.push_back( (*eit).first );
                            active_face_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                        else if( CSMP_ElementSpecifications::LineElement( (*eit).second ) )
                        {
                            active_edges.push_back( (*eit).first );
                            active_edge_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                    }
                }
                else
                {
                    for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                         eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                    {
                        if( CSMP_ElementSpecifications::VolumeElement( (*eit).second ) )
                        {
                            inactive_cells.push_back( (*eit).first );
                            inactive_cell_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                        else if( CSMP_ElementSpecifications::SurfaceElement( (*eit).second ) )
                        {
                            inactive_faces.push_back( (*eit).first );
                            inactive_face_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                        else if( CSMP_ElementSpecifications::LineElement( (*eit).second ) )
                        {
                            inactive_edge.push_back( (*eit).first );
                            inactive_edge_fem_types.insert( csmp::parseFiniteElementType( (*eit).second ) );
                        }
                    }
                }
            }
        }
    }

    if ( !exclude_inactive_cells_ )
      {
          /// all elements
          if ( !cell_elmts.empty() ) {
              model_topology.AddRegion( "ALL_CELLS", cell_fem_types, cell_elmts );
              regions.insert("ALL_CELLS");
          }

          /// adding active and inctive cell regions
          if ( !active_cells.empty() ) {
                model_topology.AddRegion("ALL_ACTIVE_CELLS", active_cell_fem_types, active_cells );
                regions.insert("ALL_ACTIVE_CELLS");
            }
          if ( !inactive_cells.empty() ) {
                model_topology.AddRegion("ALL_INACTIVE_CELLS",  inactive_cell_fem_types, inactive_cells );
                regions.insert("ALL_INACTIVE_CELLS");
            }
      }
    else
      {
          if ( !active_cells.empty() ) {
                model_topology.AddRegion("ALL_CELLS", active_cell_fem_types, active_cells );
                regions.insert("ALL_CELLS");
            }
      }

    if( csmp_error.Verbose() )
        std::cout <<"CornerPointGrid::EstablishActiveDomain: inactive cells were successfully removed.\n";

 } // EstablishActiveDomain






/**
    Identifies the faces of the cells that make up the boundary.
    The cells are then grouped into boundary regions that are stored in the ModelTopology class.
*/
void CornerPointGrid::EstablishBoundaries( csmp::VSet<3U>& vset,
                                           csmp::ModelTopology& model_topology,
                                           const std::vector<CornerPointCell>& poly )
  {
    if ( !exclude_inactive_cells_ ) {
      throw csmp::Exception( ERROR, "CornerPointGrid::EstablishBoundaryRegions",
                  "Cannot create boundary regions if we are excluding inactive cells" );
    }

    for( size_t j = 0; j < NY_; ++j ) {
      for( size_t i = 0; i < NX_; ++i ) {
        size_t kmin = NZ_;
        size_t kmax = 0;
        for( size_t k = 0; k < NZ_; ++k ) {
          size_t hexa_cell_id = CellIndex(i,j,k);
          
          if( cell_activity_[ hexa_cell_id ] == 1 ) {
            kmin = std::min(kmin, k);
            kmax = std::max(kmax, k);
          }
        }

        for( size_t k = 0; k < NZ_; ++k ) {
          /// hexa cell global id
          size_t hexa_cell_id = CellIndex(i,j,k);
          
          if( cell_activity_[ hexa_cell_id ] == 1 ) {
            const uint8_t TOP_FLAG = 1;
            const uint8_t BOTTOM_FLAG = 2;
            const uint8_t IRREGULAR_FLAG = 4;

            uint8_t node_flags[8];

            std::memset(node_flags, 0, sizeof(node_flags));

            if (i + 1 == NX_ || !cell_activity_[ CellIndex(i+1,j,k) ]) {
              node_flags[1] |= IRREGULAR_FLAG;
              node_flags[3] |= IRREGULAR_FLAG;
              node_flags[5] |= IRREGULAR_FLAG;
              node_flags[7] |= IRREGULAR_FLAG;
            }
            if (i == 0 || !cell_activity_[ CellIndex(i-1,j,k) ]) {
              node_flags[0] |= IRREGULAR_FLAG;
              node_flags[2] |= IRREGULAR_FLAG;
              node_flags[4] |= IRREGULAR_FLAG;
              node_flags[6] |= IRREGULAR_FLAG;
            }
            if (j + 1 == NY_ || !cell_activity_[ CellIndex(i,j+1,k) ]) {
              node_flags[2] |= IRREGULAR_FLAG;
              node_flags[3] |= IRREGULAR_FLAG;
              node_flags[6] |= IRREGULAR_FLAG;
              node_flags[7] |= IRREGULAR_FLAG;
            }
            if (j == 0 || !cell_activity_[ CellIndex(i,j-1,k) ]) {
              node_flags[0] |= IRREGULAR_FLAG;
              node_flags[2] |= IRREGULAR_FLAG;
              node_flags[4] |= IRREGULAR_FLAG;
              node_flags[6] |= IRREGULAR_FLAG;
            }
            if ( k == kmax ) {
              node_flags[0] |= BOTTOM_FLAG;
              node_flags[1] |= BOTTOM_FLAG;
              node_flags[2] |= BOTTOM_FLAG;
              node_flags[3] |= BOTTOM_FLAG;
            }
            if ( k == kmin ) {
              node_flags[4] |= TOP_FLAG;
              node_flags[5] |= TOP_FLAG;
              node_flags[6] |= TOP_FLAG;
              node_flags[7] |= TOP_FLAG;
            }

          auto& polygon = poly[ hexa_cell_id ];
          for ( size_t n = 0; n < 8; ++n ) {
            auto nid = polygon.GetPillarNodeGlobalIdOriginalOrder(n);
              switch (node_flags[n]) {
                  case 0:
                  {
                      break;
                  }

                  case 1: // TOP_FLAG
                  case 5: // IRREGULAR_FLAG | TOP_FLAG
                  {
                      vset.AddBFlag( nid, TOP_OUTSIDE );
                      break;
                  }

                  case 2: // BOTTOM_FLAG
                  case 6: // BOTTOM_FLAG | IRREGULAR_FLAG
                  {
                      vset.AddBFlag( nid, BOTTOM_OUTSIDE );
                      break;
                  }

                  case 4: // IRREGULAR_FLAG
                  {
                      vset.AddBFlag( nid, IRREGULAR_OUTSIDE );
                      break;
                  }
                  
                  case 3: // TOP_FLAG | BOTTOM_FLAG
                  case 7: // TOP_FLAG | BOTTOM_FLAG | IRREGULAR_FLAG
                  {
                    std::cerr << "CornerPointGrid::EstablishBoundaries: "
                      << "Possibly erroneous node " << nid << " at "
                      << vset.Px(nid) << ',' << vset.Py(nid) << ',' << vset.Pz(nid) << '\n';
                    break;
                  }
              }
            }
          }
        }
      }
    }
  }


/**
     Identifies the faces of the cells that make up a fault.
     The cells are then grouped into separate fault regions that are stored in the ModelTopology class.
     
     The faces=side walls of the fault are added to the supplied VSet.

     this breaks the record in terms of Method size and nested loops ! (@@)
     
     TODO: this function never modifies the VSet although we tried several cases with faults present; why?
*/
void CornerPointGrid::EstablishFaultRegions( csmp::VSet<3U>& vset,
                                             csmp::ModelTopology& model_topology,
                                             std::set<std::string>& faults,
                                             const std::vector<CornerPointCell>& poly )
{
    /// cell data
    const size_t quad_fem_nodes( 4U );
    size_t hexa_face_id( 0 );
    size_t hexa_cell_id( 0 );
    size_t hexa_cell_neighbor_id( 0 );
    size_t hexa_face_neighbor_id( 0 );
    std::pair<size_t,size_t> face,nface;

    /// opposite faces
    std::vector<size_t> opposite_face( 6U, 0 );
    opposite_face[CORNER_POINT_CELL_FACE_Xminus]  = CORNER_POINT_CELL_FACE_Xplus;
    opposite_face[CORNER_POINT_CELL_FACE_Xplus ]  = CORNER_POINT_CELL_FACE_Xminus;
    opposite_face[CORNER_POINT_CELL_FACE_Yminus]  = CORNER_POINT_CELL_FACE_Yplus;
    opposite_face[CORNER_POINT_CELL_FACE_Yplus ]  = CORNER_POINT_CELL_FACE_Yminus;
    opposite_face[CORNER_POINT_CELL_FACE_Zminus]  = CORNER_POINT_CELL_FACE_Zplus;
    opposite_face[CORNER_POINT_CELL_FACE_Zplus ]  = CORNER_POINT_CELL_FACE_Zminus;

    /// neighbors from each side
    std::vector<long> opposite_cell_offset( 6U, 0 );
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Xminus]  = -1;
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Xplus ]  = +1;
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Yminus]  = -NX_;
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Yplus ]  = +NX_;
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Zminus]  = -NX_x_NY_;
    opposite_cell_offset[CORNER_POINT_CELL_FACE_Zplus ]  = +NX_x_NY_;

    std::vector<std::vector<size_t> > face_nodes( 6U, std::vector<size_t>( 5, 0 ) );
    face_nodes[ CORNER_POINT_CELL_FACE_Xminus ][0] = 0;
    face_nodes[ CORNER_POINT_CELL_FACE_Xminus ][1] = 2;
    face_nodes[ CORNER_POINT_CELL_FACE_Xminus ][2] = 6;
    face_nodes[ CORNER_POINT_CELL_FACE_Xminus ][3] = 4;
    face_nodes[ CORNER_POINT_CELL_FACE_Xminus ][4] = 0;

    face_nodes[ CORNER_POINT_CELL_FACE_Xplus  ][0] = 1;
    face_nodes[ CORNER_POINT_CELL_FACE_Xplus  ][1] = 3;
    face_nodes[ CORNER_POINT_CELL_FACE_Xplus  ][2] = 7;
    face_nodes[ CORNER_POINT_CELL_FACE_Xplus  ][3] = 5;
    face_nodes[ CORNER_POINT_CELL_FACE_Xplus  ][4] = 1;

    face_nodes[ CORNER_POINT_CELL_FACE_Yminus ][0] = 0;
    face_nodes[ CORNER_POINT_CELL_FACE_Yminus ][1] = 1;
    face_nodes[ CORNER_POINT_CELL_FACE_Yminus ][2] = 5;
    face_nodes[ CORNER_POINT_CELL_FACE_Yminus ][3] = 4;
    face_nodes[ CORNER_POINT_CELL_FACE_Yminus ][4] = 0;

    face_nodes[ CORNER_POINT_CELL_FACE_Yplus  ][0] = 2;
    face_nodes[ CORNER_POINT_CELL_FACE_Yplus  ][1] = 3;
    face_nodes[ CORNER_POINT_CELL_FACE_Yplus  ][2] = 7;
    face_nodes[ CORNER_POINT_CELL_FACE_Yplus  ][3] = 6;
    face_nodes[ CORNER_POINT_CELL_FACE_Yplus  ][4] = 2;

    face_nodes[ CORNER_POINT_CELL_FACE_Zminus ][0] = 0;
    face_nodes[ CORNER_POINT_CELL_FACE_Zminus ][1] = 1;
    face_nodes[ CORNER_POINT_CELL_FACE_Zminus ][2] = 3;
    face_nodes[ CORNER_POINT_CELL_FACE_Zminus ][3] = 2;
    face_nodes[ CORNER_POINT_CELL_FACE_Zminus ][4] = 0;

    face_nodes[ CORNER_POINT_CELL_FACE_Zplus  ][0] = 4;
    face_nodes[ CORNER_POINT_CELL_FACE_Zplus  ][1] = 5;
    face_nodes[ CORNER_POINT_CELL_FACE_Zplus  ][2] = 7;
    face_nodes[ CORNER_POINT_CELL_FACE_Zplus  ][3] = 6;
    face_nodes[ CORNER_POINT_CELL_FACE_Zplus  ][4] = 4;

    std::string fault_name,fault_name_a,fault_name_b;
    std::pair<std::map<std::string,std::vector<std::pair<size_t,size_t> > >::iterator,bool> faults_data_it;
    faults.clear();
  
    // loops over the faults assuming that their cells are adjacent to discontinuities = cell interfaces that are not shared between
    // adjacent cells.
    // To deal with these lower-dimensional surfaces, they are grouped in terms of sides.
    // These sides are logged in 2 new separate regions for furthert processing.
    if( !faults_data_.empty() )
    {
// DEBUGGING
//for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
//     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
//    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
//    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
//          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
//      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
// }
//std::cerr <<"\n";

        // map of fault-name keys to vectors of cell-ID/face-ID pairs
        std::map<std::string,std::vector<std::pair<size_t,size_t> > > temp_faults_data;
        for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
             it = faults_data_.begin(); it != faults_data_.end(); ++it )
        {
            fault_name = (*it).first;

            /// in case if side is already defined
            if( fault_name.find("side") != std::string::npos )
            {
                temp_faults_data.insert( std::make_pair( fault_name, (*it).second ) );
            }
            else
            {
                /// add new fault region ( side A )
                fault_name_a  = fault_name;
                fault_name_a += "_sideA";
                faults_data_it  = temp_faults_data.insert( std::make_pair( fault_name_a, std::vector<std::pair<size_t,size_t> >() ) );
                std::vector<std::pair<size_t,size_t> >& fault_elmts_a( (*faults_data_it.first).second );

                /// add new fault region ( side B )
                fault_name_b  = fault_name;
                fault_name_b += "_sideB";
                faults_data_it  = temp_faults_data.insert( std::make_pair( fault_name_b, std::vector<std::pair<size_t,size_t> >() ) );
                std::vector<std::pair<size_t,size_t> >& fault_elmts_b( (*faults_data_it.first).second );

                // filling the cell-ID/face-ID vectors for the identified fault faces
                for( std::vector<std::pair<size_t,size_t> >::const_iterator
                     vit = (*it).second.begin(); vit != (*it).second.end(); ++vit )
                {
                    /// assign data from provided face
                    hexa_cell_id = (*vit).first;
                    if( cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_ )
                    {
                        hexa_face_id = (*vit).second;
                        if( poly[ hexa_cell_id ].GetNumPolygonFaceNodes( hexa_face_id ) > 2U )
                            fault_elmts_a.push_back( *vit );
                    }

                    /// assign data from opposite face
                    hexa_cell_neighbor_id = hexa_cell_id + opposite_cell_offset[ hexa_face_id ];
                    if( cell_activity_[hexa_cell_neighbor_id] == 1 || !exclude_inactive_cells_ )
                    {
                        hexa_face_neighbor_id = opposite_face[ hexa_face_id ];
                        if( poly[ hexa_cell_neighbor_id ].GetNumPolygonFaceNodes( hexa_face_neighbor_id ) > 2U )
                            fault_elmts_b.push_back( std::make_pair( hexa_cell_neighbor_id, hexa_face_neighbor_id ) );
                    }
                }
                if( fault_elmts_a.empty() )
                    temp_faults_data.erase( fault_name_a );
                if( fault_elmts_b.empty() )
                    temp_faults_data.erase( fault_name_b );
            }
        }
        // the old faults dataset is replaced by the new faults_data
        faults_data_.clear();
        faults_data_ = temp_faults_data;
        temp_faults_data.clear();

// DEBUGGING - appears to recreate faults data exactly = duplicating it
//std::cerr <<"\n\nrecreated fault data:\n";
//for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
//     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
//    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
//    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
//          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
//      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
// }
//std::cerr <<"\n";

    }
    else // if no fault data are provided this information is created automatically (how?)
    {
        std::set<size_t> face_nidsA;
        std::set<size_t> face_nidsB;
        std::set<size_t> edge;
        std::set<std::pair<size_t,size_t> > fault_faces;

        /// face for each direction
        std::vector<std::vector<size_t> > direction_faces(3, std::vector<size_t>( 2,0 ) );
        direction_faces[ 0 ][ 0 ] = CORNER_POINT_CELL_FACE_Xplus;
        direction_faces[ 0 ][ 1 ] = CORNER_POINT_CELL_FACE_Xminus;
        direction_faces[ 1 ][ 0 ] = CORNER_POINT_CELL_FACE_Yplus;
        direction_faces[ 1 ][ 1 ] = CORNER_POINT_CELL_FACE_Yminus;
        direction_faces[ 2 ][ 0 ] = CORNER_POINT_CELL_FACE_Zplus;
        direction_faces[ 2 ][ 1 ] = CORNER_POINT_CELL_FACE_Zminus;

        ///                  edge nodes         cell id , face id
        std::multimap<std::set<size_t>, std::pair<size_t,size_t> > nc;

        /// find fault faces in each direction
        for( size_t d = 0 ; d<2U; ++d )
        {
            face.second  = direction_faces[d][0];
            nface.second = direction_faces[d][1];
            for( size_t k = 0; k < NZ_; ++k )
                for( size_t j = 0; j < NY_; ++j )
                    for( size_t i = 0; i < NX_; ++i )
                    {
                        hexa_cell_id = CellIndex(i,j,k);
                        hexa_cell_neighbor_id  = hexa_cell_id;
                        hexa_cell_neighbor_id += opposite_cell_offset[ direction_faces[d][0] ];
                        face.first   = hexa_cell_id;
                        nface.first  = hexa_cell_neighbor_id;

                        /// only internal faces
                        if( ( i == NX_-1 && d == 0 ) || ( j == NY_-1 && d == 1 ) || ( k == NZ_-1 && d == 2 ) )
                            continue;

                        /// only faces between active cells
                        if(   (exclude_inactive_cells_ && cell_activity_[ hexa_cell_id ] == 1 && cell_activity_[ hexa_cell_neighbor_id ] == 1 )
                           || !exclude_inactive_cells_ )
                        {
                            /// comparing potentially neighboring faces
                            face_nidsA.clear();
                            for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                                face_nidsA.insert( poly[ hexa_cell_id ].GetPillarNodeGlobalIdOriginalOrder( face_nodes[face.second][nid] ) );

                            face_nidsB.clear();
                            for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                                face_nidsB.insert( poly[ hexa_cell_neighbor_id ].GetPillarNodeGlobalIdOriginalOrder( face_nodes[nface.second][nid] ) );

                            if( face_nidsA != face_nidsB )
                            {
                                if( face_nidsA.size() > 2U )
                                {
                                    fault_faces.insert( face );
                                    for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                                    {
                                        edge.clear();
                                        edge.insert( poly[face.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[face.second][nid]   ) );
                                        edge.insert( poly[face.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[face.second][nid+1] ) );
                                        if( edge.size() == 2 )
                                            nc.insert( std::make_pair( edge, face ) );
                                    }
                                }
                                if( face_nidsB.size() > 2U )
                                {
                                    fault_faces.insert( nface );
                                    for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                                    {
                                        edge.clear();
                                        edge.insert( poly[nface.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[nface.second][nid]   ) );
                                        edge.insert( poly[nface.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[nface.second][nid+1] ) );
                                        if( edge.size() == 2 )
                                            nc.insert( std::make_pair( edge, nface ) );
                                    }
                                }
                            }
                        }
                    }
        }

        // if no fault data were provided, but they were recreated, potential subsets are identified as above for the faults data (repetition of code)
        if( !fault_faces.empty() )
        {
            /// identify fault subsets
            std::string faultname_prefix = "FAULT";
            size_t      faultname_index  = 1;
            std::pair<std::multimap<std::set<size_t>,std::pair<size_t,size_t> >::iterator,std::multimap<std::set<size_t>,std::pair<size_t,size_t> >::iterator> neighbor_faces;
            std::set<std::pair<size_t,size_t> >::iterator it;
            std::set<std::pair<size_t,size_t> > neighbors;
            bool found_face( false );
            while( !fault_faces.empty() )
            {
                std::ostringstream oss;
                oss << faultname_prefix << faultname_index++;
                fault_name = oss.str();

                /// add new fault region ( side A )
                fault_name_a  = fault_name;
                fault_name_a += "_sideA";
                faults_data_it  = faults_data_.insert( std::make_pair( fault_name_a, std::vector<std::pair<size_t,size_t> >() ) );
                std::vector<std::pair<size_t,size_t> >& fault_elmts_a( (*faults_data_it.first).second );

                /// add new fault region ( side B )
                fault_name_b  = fault_name;
                fault_name_b += "_sideB";
                faults_data_it  = faults_data_.insert( std::make_pair( fault_name_b, std::vector<std::pair<size_t,size_t> >() ) );
                std::vector<std::pair<size_t,size_t> >& fault_elmts_b( (*faults_data_it.first).second );

                /// add first element
                it = fault_faces.begin();
                neighbors.insert( *it );
                fault_elmts_a.push_back( *it );
                fault_faces.erase( it );
                while( !neighbors.empty() )
                {
                    it  = neighbors.begin();
                    face = *it;
                    neighbors.erase( it );

                    for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                    {
                        edge.clear();
                        edge.insert( poly[face.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[face.second][nid]   ) );
                        edge.insert( poly[face.first].GetPillarNodeGlobalIdOriginalOrder( face_nodes[face.second][nid+1] ) );
                        if( edge.size() == 2 )
                        {
                            neighbor_faces = nc.equal_range( edge );
                            if( neighbor_faces.first != nc.end() )
                            {
                                found_face = false;
                                /// first priority: faces of current cell
                                for( std::multimap<std::set<size_t>,std::pair<size_t,size_t> >::iterator
                                     fit = neighbor_faces.first; fit != neighbor_faces.second; ++fit )
                                {
                                    nface = (*fit).second;
                                    if( nface.first == face.first && nface.second != face.second )
                                    {
                                        it = fault_faces.find( nface );
                                        if( it != fault_faces.end() )
                                        {
                                            fault_elmts_a.push_back( nface );
                                            neighbors.insert( nface );
                                            fault_faces.erase( it );
                                            found_face = true;
                                            break;
                                        }
                                    }
                                }
                                /// second priority: faces which does not belong to opposite cell
                                if( !found_face )
                                {
                                    for( std::multimap<std::set<size_t>,std::pair<size_t,size_t> >::iterator
                                         fit = neighbor_faces.first; fit != neighbor_faces.second; ++fit )
                                    {
                                        nface = (*fit).second;
                                        if( (nface.first != face.first) && (nface.first != (face.first + opposite_cell_offset[face.second]) ) )
                                        {
                                            it = fault_faces.find( nface );
                                            if( it != fault_faces.end() )
                                            {
                                                fault_elmts_a.push_back( nface );
                                                neighbors.insert( nface );
                                                fault_faces.erase( it );
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                /// add elements from opposite side
                for( std::vector<std::pair<size_t,size_t> >::const_iterator
                     fit = fault_elmts_a.begin(); fit != fault_elmts_a.end(); ++fit )
                {
                    face  = *fit;
                    nface.first  = face.first + opposite_cell_offset[ face.second ];
                    nface.second = opposite_face[ face.second ];
                    it = fault_faces.find( nface );
                    if( it != fault_faces.end() )
                    {
                        fault_elmts_b.push_back( *it );
                        fault_faces.erase( it );
                    }
                }
            }
        }

        // outputs "autogen" text file to record the fault information in Eclipse format
        if( !faults_data_.empty() )
        {
            /// face map
            std::map<size_t,std::string> facemap;
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Xminus,"X-" ) );
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Yminus,"Y-" ) );
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Zminus,"Z-" ) );
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Xplus, "X+" ) );
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Yplus, "Y+" ) );
            facemap.insert( std::make_pair( CORNER_POINT_CELL_FACE_Zplus, "Z+" ) );

            /// output
            std::string filename = model_topology.ModelName();
            filename += "-autogen-faults.txt";
            std::ofstream fout( filename, std::ofstream::out );
            fout << "-- Format      : ECLIPSE fault data (ASCII)" << std::endl;
            fout << "-- KEYWORD \"FAULTS\" HAS BEEN WRITTEN FROM CSMP++." << std::endl;
            fout << "FAULTS" << std::endl;
            fout << "-- Matrix Faults" << std::endl;
            for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
                 fit = faults_data_.begin(); fit != faults_data_.end(); ++fit )
            {
                fout << "\n\n-- NAME\tIX1\tIX2\tIY1\tIY2\tIZ1\tIZ2\tFACE" << std::endl;
                fault_name = (*fit).first;
                for( std::vector<std::pair<size_t,size_t> >::const_iterator
                     ffit = (*fit).second.begin(); ffit != (*fit).second.end(); ++ffit )
                {
                    hexa_cell_id = (*ffit).first;
                    hexa_face_id = (*ffit).second;
                    const size_t i = (  hexa_cell_id%NX_ );
                    const size_t j = ( (hexa_cell_id/NX_)%NY_ );
                    const size_t k = (  hexa_cell_id/NX_x_NY_ );
                    fout << "'" << fault_name << "'"
                         << "\t" << i << "\t" << i << "\t" << j << "\t" << j << "\t" << k << "\t" << k
                         << "\t" << "'" << facemap[ hexa_face_id ] << "'" << std::endl;
                }
            }
            fout <<'/' << std::endl;
            fout.close();
        }
    }

    /* process matching fault regions
       normally these shoud not match because the faults displace the grid.
       
       checks whether there is an offset in the faults.
       
       TODO: has this code been tested? - in all our examples it was not reached
     
       TODO: the algorithm can be extended to find out groups of mathing elements
    */
    if( !faults_data_.empty() )
    {
        std::string fault_name_prefix;
        std::map<std::string,std::vector<std::pair<size_t,size_t> > >::iterator fit;
        std::map<std::string,std::vector<std::pair<size_t,size_t> > > temp_faults_data;
        for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::iterator
             it = faults_data_.begin(); it != faults_data_.end(); ++it )
        {
            fault_name = (*it).first;

            if( fault_name.find("sideA") != std::string::npos )
            {
                // recovering the original name of the fault
                fault_name_prefix = fault_name;
                fault_name_prefix.erase(fault_name_prefix.find("_sideA"),6);
                fault_name_a = fault_name;
              
                // creates a backup copy of the contents of fault sideA
                std::set<std::set<size_t> > elmts_a;
                std::vector<std::pair<size_t,size_t> >& fault_elmts_a( (*it).second );
                for( std::vector<std::pair<size_t,size_t> >::const_iterator
                     vit = fault_elmts_a.begin(); vit != fault_elmts_a.end(); ++vit )
                  {
                      /// assign data from provided face
                      hexa_cell_id = (*vit).first;
                      hexa_face_id = (*vit).second;

                      std::set<size_t> elmnt;
                      for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                          elmnt.insert( poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder( face_nodes[hexa_face_id][nid] ) );
                      elmts_a.insert( elmnt );
                  }

                // finds the other side of the fault
                fault_name_b  = fault_name_prefix;
                fault_name_b +="_sideB";
                fit  = faults_data_.find( fault_name_b );
                std::set<std::set<size_t> > elmts_b;
                if( fit != faults_data_.end() )
                {
                    std::vector<std::pair<size_t,size_t> >& fault_elmts_b( (*fit).second );
                    for( std::vector<std::pair<size_t,size_t> >::const_iterator
                         vit = fault_elmts_b.begin(); vit != fault_elmts_b.end(); ++vit )
                    {
                        /// assign data from provided face
                        hexa_cell_id = (*vit).first;
                        hexa_face_id = (*vit).second;

                        std::set<size_t> elmnt;
                        for( size_t nid = 0; nid < quad_fem_nodes; ++nid )
                            elmnt.insert( poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder( face_nodes[hexa_face_id][nid] ) );
                        elmts_b.insert( elmnt );
                    }
                }

                // if the 2 fault sides are the same, the separation in sides a and b is removed?
                if( elmts_a == elmts_b )
                {
                    temp_faults_data.insert( std::make_pair( fault_name_prefix, (*it).second ) );
                }
                else
                {
                    temp_faults_data.insert( std::make_pair( fault_name_a, (*it).second ) );
                    temp_faults_data.insert( std::make_pair( fault_name_b, (*fit).second ) );
                }
            }
        }
      
        //  ???? here we go again?
        /// reassign faults data
        faults_data_.clear();
        faults_data_ = temp_faults_data;
        temp_faults_data.clear();
    }


// DEBUGGING - appears to recreate faults data exactly = duplicating it
//std::cerr <<"\n\n twice recreated fault data:\n";
//for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
//     it = faults_data_.begin(); it != faults_data_.end(); ++it ) {
//    std::cerr <<"\n\n"<< (*it).first <<" (cell-ID/face-ID pairs): ";
//    for ( std::vector<std::pair<size_t,size_t> >::const_iterator
//          fit=(*it).second.begin(); fit!=(*it).second.end(); fit++ )
//      std::cerr << (*fit).first <<","<< (*fit).second <<" ";
// }
//std::cerr <<"\n";



    /* writing results to ModelTopology and VSet
        
       - for each fault two regions of surface elements will be created, being called _sideA and _sideB of the 
         the given fault name
         
       - the lower-dimensional fault-side elements are being pushed back into the VSet
    */
    size_t                cell_id( vset.Elements() );
    std::vector<size_t>   fault_elmts;
    std::set<std::string> fault_fem_types;
    size_t                num_sub_faces;
  
    for( std::map<std::string,std::vector<std::pair<size_t,size_t> > >::const_iterator
        it = faults_data_.begin(); it != faults_data_.end(); ++it )
    {
        fault_name = (*it).first;
        fault_elmts.clear();
        fault_fem_types.clear();
        for( std::vector<std::pair<size_t,size_t> >::const_iterator
             fit = (*it).second.begin(); fit != (*it).second.end(); ++fit )
        {
            hexa_cell_id = (*fit).first;
            if( cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_ )
            {
                hexa_face_id  = (*fit).second;
                assert( hexa_face_id <= 5 );
// DEBUGGING: NODE ORDER OF FACES IS DIFFERENT FROM DEFINITION FOR CORNER-POINT CELL
//poly[ hexa_cell_id ].Out();

                // used to always evaluate to zero; now replaced by subclass method
                num_sub_faces = poly[ hexa_cell_id ].GetNumSubFaces( hexa_face_id );

                for( size_t fid = 0; fid<num_sub_faces; ++fid )
                {
                    const size_t num_face_nodes( poly[ hexa_cell_id ].GetNumFaceNodes( hexa_face_id, fid ) );
                    if( num_face_nodes > 2U )
                    {
                        /// fill VSet
                        vset.ResizePlist( cell_id + 1, CSMP_ElementSpecifications::NodesPerElementOfType(
                                                                   poly[ hexa_cell_id ].GetFaceType( hexa_face_id, fid ) ) );
                        vset.ResizeElementTypes( cell_id + 1 );
                        for( size_t nid = 0U; nid < num_face_nodes; ++nid )
                            vset.Plist( cell_id, nid, poly[ hexa_cell_id ].GetFaceNodeGlobalId( hexa_face_id, fid, nid ) );
                        vset.ElementType( cell_id, poly[ hexa_cell_id ].GetFaceType( hexa_face_id, fid ) );
                        /// add fault elements
                        fault_elmts.push_back( cell_id );
                        fault_fem_types.insert( csmp::parseFiniteElementType( poly[ hexa_cell_id ].GetFaceType( hexa_face_id, fid ) ) );
                        embedded_cells_[ hexa_cell_id ].push_back( std::make_pair( cell_id, poly[ hexa_cell_id ].GetFaceType( hexa_face_id, fid ) ) );
                        ++cell_id;
                    }
               }
            }
        }
      
      // if a fault has been identified, its member elements are added as a region to the model topology
      if ( !fault_elmts.empty() )
        {
            model_topology.AddRegion( fault_name.c_str(), fault_fem_types, fault_elmts );
            faults.insert( fault_name.c_str() );
        }
    }
  
} // end







void CornerPointGrid::EstablishWellRegions(  csmp::VSet<3U>& vset,
                                             csmp::ModelTopology& model_topology,
                                             std::set<std::string>& wells,
                                             const std::vector<CornerPointCell>& poly )
{
    size_t                hexa_cell_id;
    size_t                nid1,nid2;
    size_t                cell_id( vset.Elements() );
    std::string           well_name;
    std::vector<size_t>   well_elmts;
    std::set<std::string> well_fem_types;
    std::set<size_t>      edge_nodes;
    const csmp::CSMP_FEM_TYPE edge_fem_type( csmp::ISOPARAMETRIC_LINEAR_BAR );
    std::string           edge_fem_type_name( csmp::parseFiniteElementType( edge_fem_type ) );
    size_t                num_edge_nodes( 2U );

    for( std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >::const_iterator
        it = well_face_path_.begin(); it != well_face_path_.end(); ++it )
    {
        well_name = (*it).first;
        well_elmts.clear();
        well_fem_types.clear();
        for( std::vector<std::pair<size_t,std::pair<size_t,size_t> > >::const_iterator
             fit = (*it).second.begin(); fit != (*it).second.end(); ++fit )
        {
            hexa_cell_id = (*fit).first;
            if( cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_ )
            {
                /// fill VSet
                const size_t num_well_sections( poly[ hexa_cell_id ].GetNumWellSections() );
                for( size_t eid = 0; eid < num_well_sections; ++eid )
                {
                    vset.ResizePlist( cell_id + 1, num_edge_nodes );
                    vset.ResizeElementTypes( cell_id + 1 );
                    for( size_t nid = 0U; nid < num_edge_nodes; ++nid )
                        vset.Plist( cell_id, nid, poly[ hexa_cell_id ].GetWellSectionNodeGlobalId( eid, nid ) );
                    vset.ElementType( cell_id, edge_fem_type );
                    /// add well elements
                    well_elmts.push_back( cell_id );
                    well_fem_types.insert( edge_fem_type_name );
                    embedded_cells_[ hexa_cell_id ].push_back( std::make_pair( cell_id, edge_fem_type ) );
                    ++cell_id;
                }
            }
        }
        if( !well_elmts.empty() )
        {
            model_topology.AddRegion( well_name.c_str(), well_fem_types, well_elmts );
            wells.insert( well_name.c_str() );
        }
    }

    for( std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >::const_iterator
        it = well_edge_path_.begin(); it != well_edge_path_.end(); ++it )
    {
        well_name = (*it).first;
        well_elmts.clear();
        well_fem_types.clear();
        for( std::vector<std::pair<size_t,std::pair<size_t,size_t> > >::const_iterator
             fit = (*it).second.begin(); fit != (*it).second.end(); ++fit )
        {
            hexa_cell_id = (*fit).first;
            if( cell_activity_[hexa_cell_id] == 1 || !exclude_inactive_cells_ )
            {
                nid1 = (*fit).second.first;
                nid2 = (*fit).second.second;
                edge_nodes.clear();
                edge_nodes.insert( poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid1) );
                edge_nodes.insert( poly[hexa_cell_id].GetPillarNodeGlobalIdOriginalOrder(nid2) );
                num_edge_nodes = edge_nodes.size();
                if( num_edge_nodes == 2 )
                {
                    /// fill VSet
                    vset.ResizePlist( cell_id + 1, num_edge_nodes );
                    vset.ResizeElementTypes( cell_id + 1 );
                    vset.Plist( cell_id, 0, poly[ hexa_cell_id ].GetPillarNodeGlobalIdOriginalOrder( nid1 ) );
                    vset.Plist( cell_id, 1, poly[ hexa_cell_id ].GetPillarNodeGlobalIdOriginalOrder( nid2 ) );
                    vset.ElementType( cell_id, edge_fem_type );
                    /// add well elements
                    well_elmts.push_back( cell_id );
                    well_fem_types.insert( edge_fem_type_name );
                    embedded_cells_[ hexa_cell_id ].push_back( std::make_pair( cell_id, edge_fem_type ) );
                    ++cell_id;
                }
            }
        }
        if( !well_elmts.empty() )
        {
            model_topology.AddRegion( well_name.c_str(), well_fem_types, well_elmts );
            wells.insert( well_name.c_str() );
        }
    }
}


template<class VarType>
void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&             vset,
                                           const std::vector<VarType>& prop_data,
                                           const std::string&          prop_name,
                                           const csmp::PLACEMENT&      prop_place ) const
 {
    /// correcting data cell id's due to existance of embedded cells
    VarType var;
    var = 0.;
    size_t current_id = 0;
    std::vector<VarType> cell_data( vset.Elements(), var );
    if( !embedded_cells_.empty() )
    {
        ;
        for( typename std::vector<VarType>::const_iterator
             it = prop_data.begin(); it!= prop_data.end(); ++it, ++current_id )
          {
              std::map<size_t,std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> > >::const_iterator
                cit = embedded_cells_.find( current_id );
              if( cit != embedded_cells_.end() )
                  for( std::vector<std::pair<size_t,csmp::CSMP_FEM_TYPE> >::const_iterator
                       eit = (*cit).second.begin(); eit != (*cit).second.end(); ++eit )
                      cell_data[ (*eit).first  ] = *it;
          }
        if( prop_place == csmp::NODE ) {
            // Add nodal data to vset
            std::vector<VarType> nodal_data;
            // TODO: why node property?
            extrapolateElementToNodeProperty<3U,VarType>( vset, cell_data, nodal_data );
            //csmp::FEM_Data<VarType> property_values( prop_place, nodal_data );
            const size_t array_length = (var.Size() > 9U ) ? var.Size() : 0U;
            assert( array_length == 1 );
            PropertyData property_values( prop_place, VarType::VariableType, 3U, 0U );
            property_values.Reserve( nodal_data.size() );
            for ( const auto& it : nodal_data ) pushBack( property_values, it );
            vset.AddData( prop_name.c_str(), property_values );
        }
        else {
            // Add cell data to vset
            // csmp::FEM_Data<VarType> property_values( prop_place, cell_data );
            PropertyData property_values( prop_place, VarType::VariableType, 3U, 0U );
            property_values.Reserve( cell_data.size() );
            for ( const auto& it : cell_data ) pushBack( property_values, it );
            vset.AddData( prop_name.c_str(), property_values );
        }
    }
    return;
 }
 
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<ScalarVariable>&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<VectorVariable<3U> >&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<TensorVariable<3U> >&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<ArrayVariable>&, const std::string&, const csmp::PLACEMENT& ) const;
template void CornerPointGrid::WritePropertyToVSet( csmp::VSet<3U>&, const std::vector<FlaggedArrayVariable>&, const std::string&, const csmp::PLACEMENT& ) const;





/**
     j = rows
     i = colums
*/
void CornerPointGrid::Resize( size_t i_pillar_max, size_t j_pillar_max )
 {
    // vector<std::vector<Pillar> > pillars_;   ///< (subvertical) coordinate lines defining the pillars
    pillars_.resize( i_pillar_max );
    for ( auto it=pillars_.begin(); it!=pillars_.end(); it++ )
      (*it).resize( j_pillar_max );
 }
 
 
/** accessing the contained pillars
*/
Pillar& CornerPointGrid::operator()( size_t i, size_t j )
 {
    assert( i < pillars_.size() );
    assert( j < pillars_[i].size() );
   
     return pillars_[i][j];
 }
 
 







// BLOCK CENTERED GRID

BlockCenteredGrid::BlockCenteredGrid()
:dx_(3U)
{

}

BlockCenteredGrid::BlockCenteredGrid( const BlockCenteredGrid& grid )
: NX_( grid.NX_ ),
  NY_( grid.NY_ ),
  NZ_( grid.NZ_ ),
  dx_( grid.dx_ ),
  dy_( grid.dy_ ),
  dz_( grid.dz_ ),
  tops_( grid.tops_ )
{
}

BlockCenteredGrid& BlockCenteredGrid::operator=( const BlockCenteredGrid& grid )
{
    if( this != &grid )
    {
        NX_   = grid.NX_;
        NY_   = grid.NY_;
        NZ_   = grid.NZ_;
        dx_   = grid.dx_;
        dy_   = grid.dy_;
        dz_   = grid.dz_;
        tops_ = grid.tops_;
    }
    return *this;
}



BlockCenteredGrid::~BlockCenteredGrid()
{
}



void BlockCenteredGrid::Clear()
{
    /// grid
    tops_.clear();
    dx_.clear();
    dy_.clear();
    dz_.clear();
}



size_t BlockCenteredGrid::GetNumCells() const
{
    return NX_*NY_*NZ_;
}



std::vector<csmp::ScalarVariable>& BlockCenteredGrid::GetCellDepths()
{
    return tops_;
}



std::vector<ScalarVariable>& BlockCenteredGrid::GetCellSizes( size_t i )
{
    if( i == 0 )
        return dx_;
    if( i == 1 )
        return dy_;
    return dz_;
}



void BlockCenteredGrid::AssignDimensionX( size_t NX )
{
    NX_ = NX;
}



void BlockCenteredGrid::AssignDimensionY( size_t NY )
{
    NY_ = NY;
}



void BlockCenteredGrid::AssignDimensionZ( size_t NZ )
{
    NZ_ = NZ;
}




void BlockCenteredGrid::AssignCellCoordinatesToPillars( std::vector<std::vector<Pillar> >& pillars )
{
    const size_t NXY( NX_*NY_ );
    size_t cell_id(0);
    double64 x_offset(0.0);
    double64 y_offset(0.0);
    double64 z_offset(0.0);
    double64 dx(0.0);
    double64 dy(0.0);

    csmp::Point<3U> pt;
    pillars.resize( (NY_ + 1), std::vector<Pillar>( NX_ + 1 ) );

    for( size_t k = 0; k < NZ_; k++ )
    {
        /// top level
        y_offset = 0.0;
        for( size_t j = 0; j < NY_; j++ )
        {
            x_offset = 0.0;
            for( size_t i = 0; i < NX_; i++ )
            {
                cell_id = i + j*NX_ + k*NXY;

                dx = dx_[cell_id]();
                dy = dy_[cell_id]();

                z_offset = tops_[ cell_id ]();

                // assign point coordinates
                ///(0,0,0) NW,top
                pt[0] = x_offset;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][ i ].AddPoint( pt );
                ///(1,0,0) NE,top
                pt[0] = x_offset + dx;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][i+1].AddPoint( pt );
                ///(0,1,0) SW,top
                pt[0] = x_offset;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][ i ].AddPoint( pt );
                ///(1,1,0) SE,top
                pt[0] = x_offset + dx;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][i+1].AddPoint( pt );

                // assign new x offset
                x_offset += dx;
            }
            // assign new y offset ( take dy value from the last cell in x row )
            y_offset += dy;
        }
        /// bottom level
        y_offset = 0.0;
        for( size_t j = 0; j < NY_; j++ )
        {
            x_offset = 0.0;
            for( size_t i = 0; i < NX_; i++ )
            {
                cell_id = i + j*NX_ + k*NXY;

                dx = dx_[cell_id]();
                dy = dy_[cell_id]();
                double64 dz = dz_[cell_id]();

                z_offset = ( ( k != NZ_-1 ) ? tops_[ cell_id + NXY ]() : tops_[ cell_id ]() + dz );

                // assign point coordinates
                ///(0,0,1) NW,btm
                pt[0] = x_offset;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][ i ].AddPoint( pt );
                ///(1,0,1) NE,btm
                pt[0] = x_offset + dx;
                pt[1] = y_offset;
                pt[2] = z_offset;
                pillars[ j ][i+1].AddPoint( pt );
                ///(0,1,1) SW,btm
                pt[0] = x_offset;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][ i ].AddPoint( pt );
                ///(1,1,1) SE,btm
                pt[0] = x_offset + dx;
                pt[1] = y_offset + dy;
                pt[2] = z_offset;
                pillars[j+1][i+1].AddPoint( pt );

                // assign new x offset
                x_offset += dx;
            }
            // assign new y offset ( take dy value from the last cell in x row )
            y_offset += dy;
        }
    }

    for( size_t j = 0; j <= NY_; j++ )
        for( size_t i = 0; i <= NX_; i++ )
            pillars[ j ][ i ].AssignEnds( );
}






// WELLS

/// add well path based on symmetry assumption ( neighbouring cell defines the direction )
/// by default well is assumed to be vertical
void addWellPath( size_t NX, size_t NY, size_t NZ,
                  const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
{
    const size_t NX_x_NY( NX*NY );
    std::map<long,size_t> face_map;
    face_map.insert( std::make_pair( -1,       CORNER_POINT_CELL_FACE_Xminus ) );
    face_map.insert( std::make_pair( +1,       CORNER_POINT_CELL_FACE_Xplus  ) );
    face_map.insert( std::make_pair( -NX,      CORNER_POINT_CELL_FACE_Yminus ) );
    face_map.insert( std::make_pair( +NX,      CORNER_POINT_CELL_FACE_Yplus  ) );
    face_map.insert( std::make_pair( -NX_x_NY, CORNER_POINT_CELL_FACE_Zminus ) );
    face_map.insert( std::make_pair( +NX_x_NY, CORNER_POINT_CELL_FACE_Zplus  ) );

    std::pair<size_t,size_t> direction;
    std::pair<size_t,std::pair<size_t,size_t> > path;
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
    well_path.insert( std::make_pair( well_name, empty_path ) );
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >& wpath( well_path[ well_name ] );
    size_t num_cells( cell_ids.size() );
    for( size_t cid=0; cid<num_cells; ++cid )
    {
        size_t cell_id = cell_ids[ cid ];
        direction.first  = CORNER_POINT_CELL_FACE_Zminus;
        direction.second = CORNER_POINT_CELL_FACE_Zplus;
        path.first  = cell_id;
        path.second = direction;
        wpath.push_back( path );
    }
    if( wpath.size() > 1 )
    {
        size_t nid( wpath.size() - num_cells + 1 ); /// neighbour is a next cell
        for( size_t cid = 0; cid <(num_cells-1); ++cid, ++nid )
        {
            //              neighbor id          cell id
            long face_id  = wpath[ nid ].first - cell_ids[ cid ];
            if( nid-1 == 0 ) wpath[ nid-1 ].second.first = face_map[ -1 * face_id ];
            wpath[ nid-1 ].second.second = face_map[ face_id ];
            wpath[ nid ].second.first    = face_map[ -1*face_id ];
            wpath[ nid ].second.second   = face_map[ face_id ];
        }
    }
}





/// add well path with explicitly specified faces
void addWellPath( const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  const std::vector<std::pair<size_t,size_t> >& face_ids,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
{
    /// temp data
    std::pair<size_t,size_t> direction;
    std::pair<size_t,std::pair<size_t,size_t> > path;
    size_t cell_id;
    size_t num_cells( cell_ids.size() );
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
    well_path.insert( std::make_pair( well_name, empty_path ) );
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >& wpath( well_path[ well_name ] );
    for( size_t cid=0; cid<num_cells; ++cid )
    {
        cell_id = cell_ids[ cid ];
        direction.first  = face_ids[ cid ].first;
        direction.second = face_ids[ cid ].second;
        path.first  = cell_id;
        path.second = direction;
        wpath.push_back( path );
    }
}




/// add well path with explicitly specified faces ( same for all cells )
void addWellPath( const std::string& well_name,
                  const std::vector<size_t>& cell_ids,
                  std::pair<size_t,size_t> face_id,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
{
    /// temp data
    std::pair<size_t,size_t> direction;
    std::pair<size_t,std::pair<size_t,size_t> > path;
    size_t cell_id;
    size_t num_cells( cell_ids.size() );
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
    well_path.insert( std::make_pair( well_name, empty_path ) );
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >& wpath( well_path[ well_name ] );
    for( size_t i=0; i<num_cells; ++i )
    {
        cell_id = cell_ids[ i ];
        direction.first  = face_id.first;
        direction.second = face_id.second;
        path.first  = cell_id;
        path.second = direction;
        wpath.push_back( path );
    }
}




/// add well path with explicitly specified faces ( for single cell )
void addWellPath( const std::string& well_name,
                  size_t cell_id,
                  std::pair<size_t,size_t> face_id,
                  std::map<std::string,std::vector<std::pair<size_t,std::pair<size_t,size_t> > > >& well_path )
{
    /// temp data
    std::pair<size_t,size_t> direction;
    std::pair<size_t,std::pair<size_t,size_t> > path;
    direction.first  = face_id.first;
    direction.second = face_id.second;
    path.first  = cell_id;
    path.second = direction;
    std::vector<std::pair<size_t,std::pair<size_t,size_t> > >   empty_path;
    well_path.insert( std::make_pair( well_name, empty_path ) );
    well_path[ well_name ].push_back( path );
}

} // eclipse

} // end namespace csmp









