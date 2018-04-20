#include "meshingTools_UoM.h"
#include "ErrorHandler.h"

namespace csmp {

namespace eclipse {

// GATHERING CELLS, FACES AND NODES

bool checkCell( const std::vector<GridNode*>& cell )
{
    /// check uniquines
    const size_t num_cell_nodes( cell.size() );
    std::set<csmp::Point<3U> >unique_points;
    for( size_t nid = 0; nid<num_cell_nodes; ++nid )
        unique_points.insert( cell[nid]->GetPoint() );
    assert( unique_points.size() == num_cell_nodes );
    if( unique_points.size() != num_cell_nodes )
        return false;
    return true;
}





// TODO: do we need this monster vector of GridNode finite-element pairs ?
void addHexaCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                  GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4, GridNode* n5, GridNode* n6, GridNode* n7 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON;
    if( isCounterClockWiseOrientation( n0->GetPoint(), n1->GetPoint(), n2->GetPoint(), n4->GetPoint() ) )
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n1 );
        cell.second.push_back( n2 );
        cell.second.push_back( n3 );
        cell.second.push_back( n4 );
        cell.second.push_back( n5 );
        cell.second.push_back( n6 );
        cell.second.push_back( n7 );
    }
    else
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n3 );
        cell.second.push_back( n2 );
        cell.second.push_back( n1 );
        cell.second.push_back( n4 );
        cell.second.push_back( n7 );
        cell.second.push_back( n6 );
        cell.second.push_back( n5 );
    }
    checkCell( cell.second );
    cells.push_back( cell );
}





/**
    SKM improved form of function
*/
void addHexaCellTo( Cell::Vector& cells, GridNode* const n0, GridNode* const n1, GridNode* const n2, GridNode* const n3,
                                         GridNode* const n4, GridNode* const n5, GridNode* const n6, GridNode* const n7 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_HEXAHEDRON;
    if( isCounterClockWiseOrientation( n0->GetPoint(), n1->GetPoint(), n2->GetPoint(), n4->GetPoint() ) )
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n1 );
        cell.second.push_back( n2 );
        cell.second.push_back( n3 );
        cell.second.push_back( n4 );
        cell.second.push_back( n5 );
        cell.second.push_back( n6 );
        cell.second.push_back( n7 );
    }
    else
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n3 );
        cell.second.push_back( n2 );
        cell.second.push_back( n1 );
        cell.second.push_back( n4 );
        cell.second.push_back( n7 );
        cell.second.push_back( n6 );
        cell.second.push_back( n5 );
    }
    checkCell( cell.second );
    cells.push_back( cell );
  
} // end addHexaCell






void addPrismCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                   GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4, GridNode* n5 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_PRISM;
    if( isCounterClockWiseOrientation( n0->GetPoint(), n1->GetPoint(), n2->GetPoint(), n3->GetPoint() ) )
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n1 );
        cell.second.push_back( n2 );
        cell.second.push_back( n3 );
        cell.second.push_back( n4 );
        cell.second.push_back( n5 );
    }
    else
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n2 );
        cell.second.push_back( n1 );
        cell.second.push_back( n3 );
        cell.second.push_back( n5 );
        cell.second.push_back( n4 );
    }
    checkCell( cell.second );
    cells.push_back( cell );
}







void addPyramidCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_PYRAMID;
    if( isCounterClockWiseOrientation( n0->GetPoint(), n1->GetPoint(), n2->GetPoint(), n4->GetPoint() ) )
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n1 );
        cell.second.push_back( n2 );
        cell.second.push_back( n3 );
        cell.second.push_back( n4 );
    }
    else
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n3 );
        cell.second.push_back( n2 );
        cell.second.push_back( n1 );
        cell.second.push_back( n4 );
    }
    checkCell( cell.second );
    cells.push_back( cell );
}








void addTetraCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_TETRAHEDRON;
    if( isCounterClockWiseOrientation( n0->GetPoint(), n1->GetPoint(), n2->GetPoint(), n3->GetPoint() ) )
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n1 );
        cell.second.push_back( n2 );
        cell.second.push_back( n3 );
    }
    else
    {
        cell.second.push_back( n0 );
        cell.second.push_back( n2 );
        cell.second.push_back( n1 );
        cell.second.push_back( n3 );
    }
    checkCell( cell.second );
    cells.push_back( cell );
}







void addQuadCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_QUADRILATERAL;
    cell.second.push_back( n0 );
    cell.second.push_back( n1 );
    cell.second.push_back( n2 );
    cell.second.push_back( n3 );
    checkCell( cell.second );
    cells.push_back( cell );
}







void addTriCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_TRIANGLE;
    cell.second.push_back( n0 );
    cell.second.push_back( n1 );
    cell.second.push_back( n2 );
    checkCell( cell.second );
    cells.push_back( cell );
}








void addBarCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::ISOPARAMETRIC_LINEAR_BAR;
    cell.second.push_back( n0 );
    cell.second.push_back( n1 );
    checkCell( cell.second );
    cells.push_back( cell );
}







void addPointCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > cell;
    cell.first = csmp::POINT_ELEMENT;
    cell.second.push_back( n0 );
    //checkCell( cell.second );
    cells.push_back( cell );
}






void addQuadFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > face;
    face.first = csmp::ISOPARAMETRIC_LINEAR_QUADRILATERAL;
    face.second.push_back( n0 );
    face.second.push_back( n1 );
    face.second.push_back( n2 );
    face.second.push_back( n3 );
    checkCell( face.second );
    faces[fid].push_back( face );
}







void addTriFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0, GridNode* n1, GridNode* n2 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > face;
    face.first = csmp::ISOPARAMETRIC_LINEAR_TRIANGLE;
    face.second.push_back( n0 );
    face.second.push_back( n1 );
    face.second.push_back( n2 );
    checkCell( face.second );
    faces[fid].push_back( face );
}






void addBarFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0, GridNode* n1 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > face;
    face.first = csmp::ISOPARAMETRIC_LINEAR_BAR;
    face.second.push_back( n0 );
    face.second.push_back( n1 );
    checkCell( face.second );
    faces[fid].push_back( face );
}






void addPointFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0 )
{
    std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > face;
    face.first = csmp::POINT_ELEMENT;
    face.second.push_back( n0 );
    //checkCell( face.second );
    faces[fid].push_back( face );
}






void addExtraNode( std::vector<GridNode*>& extra_nodes, GridNode* n )
{
    extra_nodes.push_back( n );
}







// ADDITIONAL EDGES ON QUADRILATERAL FACES


size_t isEdgeExist( std::set<std::set<size_t> >& additional_edges,
                    std::vector<size_t>& diag_nodes,
                    GridNode* pt0,
                    GridNode* pt1,
                    GridNode* pt2,
                    GridNode* pt3 )
{
    std::set<size_t> edge02;
    edge02.insert( pt0->GetIdx() );
    edge02.insert( pt2->GetIdx() );
    const bool edge02_exist( additional_edges.find( edge02 ) != additional_edges.end() );

    std::set<size_t> edge13;
    edge13.insert( pt1->GetIdx() );
    edge13.insert( pt3->GetIdx() );
    const bool edge13_exist( additional_edges.find( edge13 ) != additional_edges.end() );

    diag_nodes.clear();
    if( edge02_exist && edge13_exist )
    {
        return 3;
    }
    else if( edge02_exist )
    {
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( edge13_exist )
    {
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }
    // no edges exist
    return 0;
}






size_t addValidAndPreferablyFirstEdge( std::set<std::set<size_t> >& additional_edges,
                                       std::vector<size_t>& diag_nodes,
                                       GridNode* pt0,
                                       GridNode* pt1,
                                       GridNode* pt2,
                                       GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    std::set<size_t> edge02;
    edge02.insert( pt0->GetIdx() );
    edge02.insert( pt2->GetIdx() );
    const bool edge02_exist( additional_edges.find( edge02 ) != additional_edges.end() );

    std::set<size_t> edge13;
    edge13.insert( pt1->GetIdx() );
    edge13.insert( pt3->GetIdx() );
    const bool edge13_exist( additional_edges.find( edge13 ) != additional_edges.end() );

    diag_nodes.clear();
    if( edge02_exist && edge13_exist )
    {
        return 3;
    }
    else if( edge02_exist )
    {
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( edge13_exist )
    {
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }

    // choose valid edge
    size_t valid_edge = chooseValidAndPreferablyFirstEdge(diag_nodes,pt0,pt1,pt2,pt3);

    if( valid_edge == 1 )
    {
        additional_edges.insert( edge02 );
        return 1;
    }
    else if( valid_edge == 2 )
    {
        additional_edges.insert( edge13 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addValidAndPreferablyFirstEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}







size_t addValidAndPreferablySecondEdge( std::set<std::set<size_t> >& additional_edges,
                                        std::vector<size_t>& diag_nodes,
                                        GridNode* pt0,
                                        GridNode* pt1,
                                        GridNode* pt2,
                                        GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    std::set<size_t> edge02;
    edge02.insert( pt0->GetIdx() );
    edge02.insert( pt2->GetIdx() );
    const bool edge02_exist( additional_edges.find( edge02 ) != additional_edges.end() );

    std::set<size_t> edge13;
    edge13.insert( pt1->GetIdx() );
    edge13.insert( pt3->GetIdx() );
    const bool edge13_exist( additional_edges.find( edge13 ) != additional_edges.end() );

    diag_nodes.clear();
    if( edge02_exist && edge13_exist )
    {
        return 3;
    }
    else if( edge02_exist )
    {
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( edge13_exist )
    {
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }

    // choose valid edge
    size_t valid_edge = chooseValidAndPreferablySecondEdge(diag_nodes,pt0,pt1,pt2,pt3);

    if( valid_edge == 1 )
    {
        additional_edges.insert( edge02 );
        return 1;
    }
    else if( valid_edge == 2 )
    {
        additional_edges.insert( edge13 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addValidAndPreferablySecondEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}








size_t addValidAndShortestOrBiggerSolidAngleEdge( std::set<std::set<size_t> >& additional_edges,
                                std::vector<size_t>& diag_nodes,
                                GridNode* pt0,
                                GridNode* pt1,
                                GridNode* pt2,
                                GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    std::set<size_t> edge02;
    edge02.insert( pt0->GetIdx() );
    edge02.insert( pt2->GetIdx() );
    const bool edge02_exist( additional_edges.find( edge02 ) != additional_edges.end() );

    std::set<size_t> edge13;
    edge13.insert( pt1->GetIdx() );
    edge13.insert( pt3->GetIdx() );
    const bool edge13_exist( additional_edges.find( edge13 ) != additional_edges.end() );

    diag_nodes.clear();
    if( edge02_exist && edge13_exist )
    {
        return 3;
    }
    else if( edge02_exist )
    {
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( edge13_exist )
    {
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }

    // choose valid edge
    size_t valid_edge = chooseValidAndShortestOrBiggerSolidAngleEdge(diag_nodes,pt0,pt1,pt2,pt3);

    if( valid_edge == 1 )
    {
        additional_edges.insert( edge02 );
        return 1;
    }
    else if( valid_edge == 2 )
    {
        additional_edges.insert( edge13 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addValidAndShortestOrBiggerSolidAngleEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}









/// avoiding ear clipping
size_t chooseValidAndPreferablyFirstEdge( std::vector<size_t>& diag_nodes,
                                          GridNode* pt0,
                                          GridNode* pt1,
                                          GridNode* pt2,
                                          GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// points
    const csmp::Point<3U>& n0 = pt0->GetPoint();
    const csmp::Point<3U>& n1 = pt1->GetPoint();
    const csmp::Point<3U>& n2 = pt2->GetPoint();
    const csmp::Point<3U>& n3 = pt3->GetPoint();

    const double angle1  = dihedralDegAngle( n1, n3, n0, n2 );
    const double angle2  = dihedralDegAngle( n0, n2, n1, n3 );
    const bool is_edge02_valid( angle1 != 0.0 );
    const bool is_edge13_valid( angle2 != 0.0 );

    diag_nodes.clear();
    if( is_edge02_valid && is_edge13_valid )
    {
        /// both edges are valid
        if( ( angle1 == 180.0 ) && ( angle2 == 180.0 ) )
        {
            /// choose first edge as desired
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            return 1;
        }
        else if( angle1 > angle2 )
        {
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            return 1;
        }
        else
        {
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            return 2;
        }
    }
    else if( is_edge02_valid )
    {
        /// edge 13 is not valid!!! choosing edge 02.
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( is_edge13_valid )
    {
        /// edge 02 is not valid!!! choosing edge 13.
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "chooseValidAndPreferablyFirstEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}








size_t chooseValidAndPreferablySecondEdge( std::vector<size_t>& diag_nodes,
                                           GridNode* pt0,
                                           GridNode* pt1,
                                           GridNode* pt2,
                                           GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// points
    const csmp::Point<3U>& n0 = pt0->GetPoint();
    const csmp::Point<3U>& n1 = pt1->GetPoint();
    const csmp::Point<3U>& n2 = pt2->GetPoint();
    const csmp::Point<3U>& n3 = pt3->GetPoint();

    const double angle1  = dihedralDegAngle( n1, n3, n0, n2 );
    const double angle2  = dihedralDegAngle( n0, n2, n1, n3 );
    const bool is_edge02_valid( angle1 != 0.0 );
    const bool is_edge13_valid( angle2 != 0.0 );

    diag_nodes.clear();
    if( is_edge02_valid && is_edge13_valid )
    {
        /// both edges are valid
        if( ( angle1 == 180.0 ) && ( angle2 == 180.0 ) )
        {
            /// choose second edge as desired
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            return 2;
        }
        else if( angle1 > angle2 )
        {
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            return 1;
        }
        else
        {
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            return 2;
        }
    }
    else if( is_edge02_valid )
    {
        /// edge 13 is not valid!!! choosing edge 02.
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( is_edge13_valid )
    {
        /// edge 02 is not valid!!! choosing edge 13.
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "chooseValidAndPreferablySecondEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}









size_t chooseValidAndShortestOrBiggerSolidAngleEdge( std::vector<size_t>& diag_nodes,
                                                     GridNode* pt0,
                                                     GridNode* pt1,
                                                     GridNode* pt2,
                                                     GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// points
    const csmp::Point<3U>& n0 = pt0->GetPoint();
    const csmp::Point<3U>& n1 = pt1->GetPoint();
    const csmp::Point<3U>& n2 = pt2->GetPoint();
    const csmp::Point<3U>& n3 = pt3->GetPoint();

    const double angle1  = dihedralDegAngle( n1, n3, n0, n2 );
    const double angle2  = dihedralDegAngle( n0, n2, n1, n3 );
    const bool is_edge02_valid( angle1 != 0.0 );
    const bool is_edge13_valid( angle2 != 0.0 );

    diag_nodes.clear();
    if( is_edge02_valid && is_edge13_valid )
    {
        /// both edges are valid
        if( ( angle1 == 180.0 ) && ( angle2 == 180.0 ) )
        {
            /// choose the shortest one
            csmp::Point<3U> e02 = n2; e02 -= n0;
            csmp::Point<3U> e13 = n3; e13 -= n1;
            if( e02.Length() < e13.Length() )
            {
                diag_nodes.push_back( 0 );
                diag_nodes.push_back( 2 );
                diag_nodes.push_back( 1 );
                diag_nodes.push_back( 3 );
                return 1;
            }
            else
            {
                diag_nodes.push_back( 1 );
                diag_nodes.push_back( 3 );
                diag_nodes.push_back( 0 );
                diag_nodes.push_back( 2 );
                return 2;
            }
        }
        else if( angle1 > angle2 )
        {
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            return 1;
        }
        else
        {
            diag_nodes.push_back( 1 );
            diag_nodes.push_back( 3 );
            diag_nodes.push_back( 0 );
            diag_nodes.push_back( 2 );
            return 2;
        }
    }
    else if( is_edge02_valid )
    {
        /// edge 13 is not valid!!! choosing edge 02.
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        return 1;
    }
    else if( is_edge13_valid )
    {
        /// edge 02 is not valid!!! choosing edge 13.
        diag_nodes.push_back( 1 );
        diag_nodes.push_back( 3 );
        diag_nodes.push_back( 0 );
        diag_nodes.push_back( 2 );
        return 2;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "chooseValidAndShortestOrBiggerSolidAngleEdge()",
                      "Both edges exist!!! Unpredictable case!!! ");
    return 0;
}








// ADDITIONAL POINTS

/**

   @todo SKM method should return an iterator rather than a naked pointer to a map element

*/
bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points,
                   GridNode*& pt,
                   GridNode* const pt0,
                   GridNode* const pt1,
                   GridNode* const pt2,
                   GridNode* const pt3 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );
    face.insert( pt2->GetIdx() );
    face.insert( pt3->GetIdx() );

    typename std::map<std::set<size_t>,GridNode*>::const_iterator
            it = additional_points.find( face );
    if( it != additional_points.end() )
    {
        pt = (*it).second;
        return true;
    }
    return false;
}





/**
    No references to pointers please !
*/
bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points,
                   GridNode*& pt,
                   GridNode* const pt0,
                   GridNode* const pt1,
                   GridNode* const pt2 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );
    face.insert( pt2->GetIdx() );

    typename std::map<std::set<size_t>,GridNode*>::const_iterator
            it = additional_points.find( face );
    if( it != additional_points.end() )
    {
        pt = (*it).second;
        return true;
    }
    return false;
}







bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points,
                   GridNode*& pt,
                   GridNode* const pt0,
                   GridNode* const pt1 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );

    typename std::map<std::set<size_t>,GridNode*>::const_iterator
            it = additional_points.find( face );
    if( it != additional_points.end() )
    {
        pt = (*it).second;
        return true;
    }
    return false;
}









bool containRemeshedFaces( const std::map<std::set<size_t>,GridNode*>& additional_points,
                             const std::vector<GridNode*>& nodes,
                             const std::vector<std::vector<size_t> >& face_nodes )
{
    const size_t num_faces( face_nodes.size() );

    GridNode* gn(nullptr);
    for( size_t fid = 0; fid < num_faces; ++fid )
    {
        const size_t num_face_nodes( face_nodes[fid].size() );
        if( num_face_nodes == 4U )
        {
            if( doesPointExist( additional_points,gn,
                              nodes[face_nodes[fid][0]],
                              nodes[face_nodes[fid][1]],
                              nodes[face_nodes[fid][2]],
                              nodes[face_nodes[fid][3]] ) )
            {
                return true;
            }
        }
        else if( num_face_nodes == 3U )
        {
            if( doesPointExist( additional_points,gn,
                              nodes[face_nodes[fid][0]],
                              nodes[face_nodes[fid][1]],
                              nodes[face_nodes[fid][2]] ) )
            {
                return true;
            }
        }
    }
    return false;
}








bool containRemeshedEdges( const std::map<std::set<size_t>,GridNode*>& additional_points,
                             const std::vector<GridNode*>& nodes,
                             const std::vector<std::vector<size_t> >& face_nodes )
{
    const size_t num_faces( face_nodes.size() );

    GridNode* gn(nullptr);
    for( size_t fid = 0; fid < num_faces; ++fid )
    {
        const size_t num_face_nodes( face_nodes[fid].size() );
        for( size_t nid = 0; nid < num_face_nodes; ++nid )
        {
            if( doesPointExist( additional_points,gn,
                              nodes[face_nodes[fid][nid]],
                              nodes[face_nodes[fid][(nid+1)%num_face_nodes]] ) )
            {
                return true;
            }
        }
    }
    return false;
}







void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                              GridNode* pt,
                              GridNode* pt0,
                              GridNode* pt1,
                              GridNode* pt2,
                              GridNode* pt3 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );
    face.insert( pt2->GetIdx() );
    face.insert( pt3->GetIdx() );
    additional_points.insert( std::make_pair( face, pt ) );
}







void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                              GridNode* pt,
                              GridNode* pt0,
                              GridNode* pt1,
                              GridNode* pt2 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );
    face.insert( pt2->GetIdx() );
    additional_points.insert( std::make_pair( face, pt ) );
}







void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                              GridNode* pt,
                              GridNode* pt0,
                              GridNode* pt1 )
{
    std::set<size_t> face;
    face.insert( pt0->GetIdx() );
    face.insert( pt1->GetIdx() );
    additional_points.insert( std::make_pair( face, pt ) );
}








bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                         GridNode& gn_ref,
                         GridNode* pt0,
                         GridNode* pt1,
                         GridNode* pt2,
                         GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// check whether point exist
    GridNode* gn_ptr(nullptr);
    if( doesPointExist( additional_points, gn_ptr, pt0, pt1, pt2, pt3 ) )
    {
        gn_ref = *gn_ptr;
        return true;
    }

    /// choose valid adge
    std::vector<size_t> diag_nodes;
    if( chooseValidAndShortestOrBiggerSolidAngleEdge(diag_nodes,pt0,pt1,pt2,pt3) != 0 )
    {
        csmp::Point<3U> pt(0.);
        pt += ( diag_nodes[0] == 0 ? pt0->GetPoint() : pt1->GetPoint() );
        pt += ( diag_nodes[1] == 2 ? pt2->GetPoint() : pt3->GetPoint() );
        pt /= 2.0;
        gn_ref.AssignNullIdx();
        gn_ref.AssignPoint( pt );
        return true;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "chooseFaceCentroid()",
                      "The quadrilateral is unacceptable!!! ");
    return false;
}







bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                         GridNode& gn_ref,
                         GridNode* pt0,
                         GridNode* pt1,
                         GridNode* pt2 )
{
    /// check whether point exist
    GridNode* gn_ptr(nullptr);
    if( doesPointExist( additional_points, gn_ptr, pt0, pt1, pt2 ) )
    {
        gn_ref = *gn_ptr;
        return true;
    }

    /// add new point
    csmp::Point<3U> pt(0.);
    pt += pt0->GetPoint();
    pt += pt1->GetPoint();
    pt += pt2->GetPoint();
    pt /= 3.0;
    gn_ref.AssignNullIdx();
    gn_ref.AssignPoint( pt );
    return true;
}







bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points,
                         GridNode& gn_ref,
                         GridNode* pt0,
                         GridNode* pt1 )
{
    /// check whether point exist
    GridNode* gn_ptr(nullptr);
    if( doesPointExist( additional_points, gn_ptr, pt0, pt1 ) )
    {
        gn_ref = *gn_ptr;
        return true;
    }

    /// add new point
    csmp::Point<3U> pt(0.);
    pt += pt0->GetPoint();
    pt += pt1->GetPoint();
    pt /= 2.0;
    gn_ref.AssignNullIdx();
    gn_ref.AssignPoint( pt );
    return true;
}








bool addFaceCentroid( PolygonGridManager* pgm,
                      std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& gn_ptr,
                      GridNode* pt0,
                      GridNode* pt1,
                      GridNode* pt2,
                      GridNode* pt3 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// choose face centroid
    GridNode gn_ref;
    if( chooseFaceCentroid(additional_points,gn_ref,pt0,pt1,pt2,pt3) )
    {
        /// add new node
        if( gn_ref.IsNullIdx() )
            gn_ptr = pgm->AddNode( gn_ref.GetPoint() );
        else
            gn_ptr = pgm->GetNode( gn_ref.GetIdx() );
        std::set<size_t> face;
        face.insert( pt0->GetIdx() );
        face.insert( pt1->GetIdx() );
        face.insert( pt2->GetIdx() );
        face.insert( pt3->GetIdx() );
        additional_points.insert( std::make_pair( face, gn_ptr ) );
        return true;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addNewFaceCentroid()",
                      "The quadrilateral is unacceptable!!! ");
    return false;
}








bool addFaceCentroid( PolygonGridManager* pgm,
                      std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& gn_ptr,
                      GridNode* pt0,
                      GridNode* pt1,
                      GridNode* pt2 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// choose face centroid
    GridNode gn_ref;
    if( chooseFaceCentroid(additional_points,gn_ref,pt0,pt1,pt2) )
    {
        /// add new node
        if( gn_ref.IsNullIdx() )
            gn_ptr = pgm->AddNode( gn_ref.GetPoint() );
        else
            gn_ptr = pgm->GetNode( gn_ref.GetIdx() );
        std::set<size_t> face;
        face.insert( pt0->GetIdx() );
        face.insert( pt1->GetIdx() );
        face.insert( pt2->GetIdx() );
        additional_points.insert( std::make_pair( face, gn_ptr ) );
        return true;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addNewFaceCentroid()",
                      "The triangle is unacceptable!!! ");
    return false;
}








bool addFaceCentroid( PolygonGridManager* pgm,
                      std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& gn_ptr,
                      GridNode* pt0,
                      GridNode* pt1 )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// choose face centroid
    GridNode gn_ref;
    if( chooseFaceCentroid(additional_points,gn_ref,pt0,pt1) )
    {
        /// add new node
        if( gn_ref.IsNullIdx() )
            gn_ptr = pgm->AddNode( gn_ref.GetPoint() );
        else
            gn_ptr = pgm->GetNode( gn_ref.GetIdx() );
        std::set<size_t> face;
        face.insert( pt0->GetIdx() );
        face.insert( pt1->GetIdx() );
        additional_points.insert( std::make_pair( face, gn_ptr ) );
        return true;
    }
    // both edges are not valid
    csmp_error.notice(csmp::ERROR,
                      "addNewFaceCentroid()",
                      "The triangle is unacceptable!!! ");
    return false;
}










// EXTRA USEFUL MESHING TOOLS

bool addOverlappingCellFaceCentroids( PolygonGridManager* pgm,
                                      std::map<std::set<size_t>,GridNode*>& additional_points,
                                      std::set<std::set<size_t> >& additional_edges,
                                      GridNode*& cgn,
                                      const std::vector<GridNode*>& nodes,
                                      const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_faces( face_nodes.size() );
    const size_t num_face_nodes( num_faces != 0 ? face_nodes[0].size() : 0 );

    if( num_face_nodes == 4 )
    {
        /// find one face centroid
        std::vector<size_t> diag_nodes;
        size_t valid_edge;
        long found_point( -1 );
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( doesPointExist( additional_points,cgn,
                              nodes[face_nodes[fid][0]],
                              nodes[face_nodes[fid][1]],
                              nodes[face_nodes[fid][2]],
                              nodes[face_nodes[fid][3]] ) )
            {
                found_point = fid;
                break;
            }
        }
        if( found_point == -1 )
        {
            found_point = 0;
            size_t fid = 0;
            valid_edge = isEdgeExist( additional_edges,diag_nodes,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]],
                                      nodes[face_nodes[fid][2]],
                                      nodes[face_nodes[fid][3]] );
            if( valid_edge == 3 )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addOverlappingCellFaceCentroids()",
                                  "Face contains two diaginals, which can lead to non-conformity!!! ");
                return false;
            }
            else if( valid_edge == 1 || valid_edge == 2 )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addOverlappingCellFaceCentroids()",
                                  "Face contains diaginal, cannot add face centroid!!! ");
                return false;
            }
            else if( valid_edge == 0 )
            {
                if( !addFaceCentroid( pgm,additional_points,cgn,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]],
                                      nodes[face_nodes[fid][2]],
                                      nodes[face_nodes[fid][3]] ) )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addOverlappingCellFaceCentroids()",
                                      "Cannot define appropriate quad-face centroid!!! ");
                    return false;
                }
            }
        }

        /// assign same face centroid to overlapping faces
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( fid != found_point )
            {
                valid_edge = isEdgeExist( additional_edges,diag_nodes,
                                          nodes[face_nodes[fid][0]],
                                          nodes[face_nodes[fid][1]],
                                          nodes[face_nodes[fid][2]],
                                          nodes[face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addOverlappingCellFaceCentroids()",
                                      "Face contains two diaginals, which can lead to non-conformity!!! ");
                    return false;
                }
                else if( valid_edge == 1 || valid_edge == 2 )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addOverlappingCellFaceCentroids()",
                                      "Face contains diaginal, cannot add face centroid!!! ");
                    return false;
                }
                else if( valid_edge == 0 )
                {
                    addExistingFaceCentroid( additional_points,cgn,
                                             nodes[face_nodes[fid][0]],
                                             nodes[face_nodes[fid][1]],
                                             nodes[face_nodes[fid][2]],
                                             nodes[face_nodes[fid][3]] );
                }
            }
        }
    }
    else if( num_face_nodes == 3 )
    {
        /// find one face centroid
        long found_point( -1 );
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( doesPointExist( additional_points,cgn,
                              nodes[face_nodes[fid][0]],
                              nodes[face_nodes[fid][1]],
                              nodes[face_nodes[fid][2]] ) )
            {
                found_point = fid;
                break;
            }
        }
        if( found_point == -1 )
        {
            found_point = 0;
            size_t fid = 0;
            if( !addFaceCentroid( pgm,additional_points,cgn,
                                  nodes[face_nodes[fid][0]],
                                  nodes[face_nodes[fid][1]],
                                  nodes[face_nodes[fid][2]] ) )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addOverlappingCellFaceCentroids()",
                                  "Cannot define appropriate tri-face centroid!!! ");
                return false;
            }
        }

        /// assign same face centroid to overlapping faces
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( fid != found_point )
            {
                addExistingFaceCentroid( additional_points,cgn,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]],
                                      nodes[face_nodes[fid][2]] );
            }
        }
    }
    else if( num_face_nodes == 2 )
    {
        /// find one face centroid
        long found_point( -1 );
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( doesPointExist( additional_points,cgn,
                              nodes[face_nodes[fid][0]],
                              nodes[face_nodes[fid][1]]) )
            {
                found_point = fid;
                break;
            }
        }
        if( found_point == -1 )
        {
            found_point = 0;
            size_t fid = 0;
            if( !addFaceCentroid( pgm,additional_points,cgn,
                                  nodes[face_nodes[fid][0]],
                                  nodes[face_nodes[fid][1]] ) )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addOverlappingCellFaceCentroids()",
                                  "Cannot define appropriate bar-face centroid!!! ");
                return false;
            }
        }

        /// assign same face centroid to overlapping faces
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            if( fid != found_point )
            {
                addExistingFaceCentroid( additional_points,cgn,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]] );
            }
        }
    }
    else if( num_face_nodes == 1 )
        cgn = nodes[0];
    return true;
}










bool addCellFaceCentroids( PolygonGridManager* pgm,
                           std::map<std::set<size_t>,GridNode*>& additional_points,
                           std::set<std::set<size_t> >& additional_edges,
                           const std::vector<GridNode*>& nodes,
                           const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_faces( face_nodes.size() );

    /// add cell face centroid points
    size_t valid_edge;
    std::vector<size_t> diag_nodes;
    GridNode* gn;
    for( size_t fid = 0; fid < num_faces; ++fid )
    {
        const size_t num_face_nodes( face_nodes[fid].size() );
        if( num_face_nodes == 4U )
        {
            if( !doesPointExist( additional_points,gn,
                               nodes[face_nodes[fid][0]],
                               nodes[face_nodes[fid][1]],
                               nodes[face_nodes[fid][2]],
                               nodes[face_nodes[fid][3]] ) )
            {
                valid_edge = isEdgeExist( additional_edges,diag_nodes,
                                          nodes[face_nodes[fid][0]],
                                          nodes[face_nodes[fid][1]],
                                          nodes[face_nodes[fid][2]],
                                          nodes[face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellFaceCentroids()",
                                      "Face contains two diaginals, which can lead to non-conformity!!! ");
                    return false;
                }
                else if( valid_edge == 1 || valid_edge == 2 )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellFaceCentroids()",
                                      "Face contains diaginal, cannot add face centroid!!! ");
                    return false;
                }
                else if( valid_edge == 0 )
                {
                    if( !addFaceCentroid( pgm,additional_points,gn,
                                          nodes[face_nodes[fid][0]],
                                          nodes[face_nodes[fid][1]],
                                          nodes[face_nodes[fid][2]],
                                          nodes[face_nodes[fid][3]] ) )
                    {
                        csmp_error.notice(csmp::ERROR,
                                          "addCellFaceCentroids()",
                                          "Cannot define appropriate quad-face centroid!!! ");
                        return false;
                    }
                }
            }
        }
        else if( num_face_nodes == 3U )
        {
            if( !doesPointExist( additional_points,gn,
                               nodes[face_nodes[fid][0]],
                               nodes[face_nodes[fid][1]],
                               nodes[face_nodes[fid][2]] ) )
            {
                if( !addFaceCentroid( pgm,additional_points,gn,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]],
                                      nodes[face_nodes[fid][2]] ) )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellFaceCentroids()",
                                      "Cannot define appropriate tri-face centroid!!! ");
                    return false;
                }
            }
        }
        else if( num_face_nodes == 2U )
        {
            if( !doesPointExist( additional_points,gn,
                               nodes[face_nodes[fid][0]],
                               nodes[face_nodes[fid][1]] ) )
            {
                if( !addFaceCentroid( pgm,additional_points,gn,
                                      nodes[face_nodes[fid][0]],
                                      nodes[face_nodes[fid][1]] ) )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellFaceCentroids()",
                                      "Cannot define appropriate bar-face centroid!!! ");
                    return false;
                }
            }
        }
    }
    return true;
}










/** 
    Adding cell centroids based on topology:
    hexahedron, prism: opposite quadrilateral face nodes
    pyramid: quad face nodes and apex node
    tetrahedron, quadrilateral, triangle, bar: nodes
    
    @todo remove reference to pointer in argument list, it's simply too dangerous
*/
bool addCellCentroid( PolygonGridManager* pgm,
                      std::map<std::set<size_t>,GridNode*>& additional_points,
                      std::set<std::set<size_t> >& additional_edges,
                      GridNode*& cgn,
                      const std::vector<GridNode*>& nodes,
                      const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_faces( face_nodes.size() );
    const size_t num_nodes( nodes.size() );

    std::vector<std::vector<size_t> > quad_face_nodes;
    for( size_t fid = 0; fid<num_faces; ++fid )
    {
        const size_t num_face_nodes( face_nodes.size() );
        if( num_face_nodes == 4 )
        {
            std::vector<size_t> fnodes;
            for( size_t nid = 0; nid < num_face_nodes; ++nid )
                fnodes.push_back( face_nodes[fid][nid] );
            quad_face_nodes.push_back(fnodes);
        }
    }

    /// add cell centroid point
    const size_t num_quad_faces( quad_face_nodes.size() );
    std::vector<csmp::Point<3U> > axis;
    if( num_quad_faces > 1 )
    {
        /// hexahedron, prism cases
        size_t valid_edge;
        std::vector<size_t>  diag_nodes;
        GridNode* gn_ptr(0);
        GridNode  gn_ref;
        csmp::Point<3U>     pt;
        for( size_t fid = 0; fid < num_quad_faces; fid += ( num_quad_faces - 1 ) )
        {
            if( doesPointExist( additional_points,gn_ptr,
                              nodes[quad_face_nodes[fid][0]],
                              nodes[quad_face_nodes[fid][1]],
                              nodes[quad_face_nodes[fid][2]],
                              nodes[quad_face_nodes[fid][3]] ) )
            {
                axis.push_back( gn_ptr->GetPoint() );
            }
            else
            {
                valid_edge = isEdgeExist( additional_edges,diag_nodes,
                                          nodes[quad_face_nodes[fid][0]],
                                          nodes[quad_face_nodes[fid][1]],
                                          nodes[quad_face_nodes[fid][2]],
                                          nodes[quad_face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellCentroid()",
                                      "Face contains two diaginals, which can lead to non-conformity!!! ");
                    return false;
                }
                else if( valid_edge == 1 || valid_edge == 2 )
                {
                    pt  = 0.0;
                    pt += nodes[quad_face_nodes[fid][diag_nodes[0]]]->GetPoint();
                    pt += nodes[quad_face_nodes[fid][diag_nodes[1]]]->GetPoint();
                    pt /= 2.0;
                    axis.push_back( pt );
                }
                else if( !chooseFaceCentroid( additional_points,gn_ref,
                                              nodes[quad_face_nodes[fid][0]],
                                              nodes[quad_face_nodes[fid][1]],
                                              nodes[quad_face_nodes[fid][2]],
                                              nodes[quad_face_nodes[fid][3]] ) )
                {
                    csmp_error.notice(csmp::ERROR,
                                      "addCellCentroid()",
                                      "Cannot define appropriate quad-face centroid!!! ");
                    return false;
                }
                else
                {
                    axis.push_back( gn_ref.GetPoint() );
                }
            }
        }

        /// assign centroid as the middle point of axis
        pt  = 0.0;
        const size_t num_axis_nodes( axis.size() );
        for( size_t i = 0; i < num_axis_nodes; ++i )
            pt += axis[i];
        pt /= static_cast<double>( num_axis_nodes );
        cgn = pgm->AddNode( pt );
    }
    else if( num_quad_faces == 1 )
    {
        /// pyramid, quadrilateral cases

        size_t valid_edge;
        std::vector<size_t>  diag_nodes;
        GridNode* gn_ptr(0);
        GridNode  gn_ref;
        csmp::Point<3U>     pt;
        size_t fid = 0;
        if( doesPointExist( additional_points,gn_ptr,
                          nodes[quad_face_nodes[fid][0]],
                          nodes[quad_face_nodes[fid][1]],
                          nodes[quad_face_nodes[fid][2]],
                          nodes[quad_face_nodes[fid][3]] ) )
        {
            axis.push_back( gn_ptr->GetPoint() );
        }
        else
        {
            valid_edge = isEdgeExist( additional_edges,diag_nodes,
                                      nodes[quad_face_nodes[fid][0]],
                                      nodes[quad_face_nodes[fid][1]],
                                      nodes[quad_face_nodes[fid][2]],
                                      nodes[quad_face_nodes[fid][3]] );
            if( valid_edge == 3 )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addCellCentroid()",
                                  "Face contains two diaginals, which can lead to non-conformity!!! ");
                return false;
            }
            else if( valid_edge == 1 || valid_edge == 2 )
            {
                pt  = 0.0;
                pt += nodes[quad_face_nodes[fid][diag_nodes[0]]]->GetPoint();
                pt += nodes[quad_face_nodes[fid][diag_nodes[1]]]->GetPoint();
                pt /= 2.0;
                axis.push_back( pt );
            }
            else if( !chooseFaceCentroid( additional_points,gn_ref,
                                          nodes[quad_face_nodes[fid][0]],
                                          nodes[quad_face_nodes[fid][1]],
                                          nodes[quad_face_nodes[fid][2]],
                                          nodes[quad_face_nodes[fid][3]] ) )
            {
                csmp_error.notice(csmp::ERROR,
                                  "addCellCentroid()",
                                  "Cannot define appropriate quad-face centroid!!! ");
                return false;
            }
            else
            {
                axis.push_back( gn_ref.GetPoint() );
            }
        }

        if( num_nodes != 4U )
        {
            for( size_t nid = 0; nid < num_nodes; ++nid )
                if( nodes[nid]->GetIdx() != nodes[quad_face_nodes[fid][0]]->GetIdx() &&
                    nodes[nid]->GetIdx() != nodes[quad_face_nodes[fid][1]]->GetIdx() &&
                    nodes[nid]->GetIdx() != nodes[quad_face_nodes[fid][2]]->GetIdx() &&
                    nodes[nid]->GetIdx() != nodes[quad_face_nodes[fid][3]]->GetIdx() )
                {
                    axis.push_back( nodes[nid]->GetPoint() );
                    break;
                }
        }

        /// assign centroid as the middle point of axis
        pt  = 0.0;
        const size_t num_axis_nodes( axis.size() );
        for( size_t i = 0; i < num_axis_nodes; ++i )
            pt += axis[i];
        pt /= static_cast<double>( num_axis_nodes );
        cgn = pgm->AddNode( pt );
    }
    else if( num_nodes > 1 )
    {
        /// tetrahedron, triangle, bar cases
        csmp::Point<3U>    pt;
        pt = 0.0;
        for( size_t nid = 0; nid<num_nodes; ++nid )
            pt += nodes[nid]->GetPoint();
        pt /= static_cast<double>( num_nodes );
        cgn = pgm->AddNode( pt );
    }
    else if( num_nodes == 1 )
        cgn = nodes[0];
    return true;
}











// PROCESSING CELLS


/**
    processing cells (with\without) cell centroids, (with\without) face centroids, (with\without) edge centroids
    cases without cell centroids can be proceeded only for restricted set of cell types
*/
bool processCellWithCentroids( bool add_edge_centroids,
                               bool add_face_centroids,
                               bool add_cell_centroids,
                               bool is_volumetric_element,
                               bool tetra_mesh,
                               PolygonGridManager* pgm,
                               std::map<std::set<size_t>,GridNode*>& additional_points,
                               std::set<std::set<size_t> >& additional_edges,
                               std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                               std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                               std::vector<GridNode* >& extra_nodes,
                               const std::vector<GridNode*>& nodes,
                               const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_nodes( nodes.size() );
    const size_t num_faces( face_nodes.size() );

    GridNode* cgn(nullptr);

    /// 0.0 process degenerate cases ( < 2 nodes )
    if( num_nodes == 0 )
    {
        return true;
    }
    else if( num_nodes == 1 )
    {
        addPointCell(cells,nodes[0]);
        return true;
    }
    else if( num_nodes == 2 )
    {
        bool edge_centroid_exist( false );
        if( doesPointExist( additional_points,cgn,
                          nodes[0],
                          nodes[1] ) )
        {
            edge_centroid_exist = true;
        }
        else if( add_edge_centroids || add_face_centroids || add_cell_centroids )
        {
            addFaceCentroid( pgm,additional_points,cgn,
                             nodes[0],
                             nodes[1] );
            edge_centroid_exist = true;
        }
        if( edge_centroid_exist )
        {
            addBarCell( cells, nodes[0], cgn );
            addBarCell( cells, cgn, nodes[1] );
        }
        else
        {
            addBarCell( cells, nodes[0], nodes[1] );
        }
        return true;
    }

    /// 1.0 add cell centroid if required
    if( add_cell_centroids )
    {
        addCellCentroid(pgm,additional_points,additional_edges,cgn,nodes,face_nodes);
        addExtraNode( extra_nodes, cgn );
    }

    /// 2.0 add face centroids if required
    GridNode* gn(nullptr);
    std::vector<GridNode*> face_cgn( num_faces, NULL );
    std::vector<GridNode*> quad_face_gn;
    std::vector<GridNode*> tri_face_gn;
    bool face_centroid_exist( false );
    for( size_t fid = 0; fid<num_faces; ++fid )
    {
        const size_t num_face_nodes( face_nodes[fid].size() );
        face_centroid_exist = false;
        if( num_face_nodes == 4 )
        {
            face_centroid_exist = doesPointExist( additional_points,gn,
                                                nodes[face_nodes[fid][0]],
                                                nodes[face_nodes[fid][1]],
                                                nodes[face_nodes[fid][2]],
                                                nodes[face_nodes[fid][3]] );
            if( !face_centroid_exist && add_face_centroids )
            {
                addFaceCentroid( pgm,additional_points,gn,
                                 nodes[face_nodes[fid][0]],
                                 nodes[face_nodes[fid][1]],
                                 nodes[face_nodes[fid][2]],
                                 nodes[face_nodes[fid][3]] );
                face_centroid_exist = true;
            }
            quad_face_gn.push_back( gn );
        }
        else if( num_face_nodes == 3 )
        {
            face_centroid_exist = doesPointExist( additional_points,gn,
                                                nodes[face_nodes[fid][0]],
                                                nodes[face_nodes[fid][1]],
                                                nodes[face_nodes[fid][2]] );
            if( !face_centroid_exist && add_face_centroids )
            {
                addFaceCentroid( pgm,additional_points,gn,
                                 nodes[face_nodes[fid][0]],
                                 nodes[face_nodes[fid][1]],
                                 nodes[face_nodes[fid][2]] );
                face_centroid_exist = true;
            }
            tri_face_gn.push_back( gn );
        }
        if( face_centroid_exist )
        {
            face_cgn[fid] = gn;
            addExtraNode( extra_nodes, gn );
        }
    }

    /// face categories
    const size_t num_quad_faces( quad_face_gn.size() );
    const size_t num_tri_faces ( tri_face_gn.size() );

    /// 3.0 choose appropriate connecting point: cell centroid or face centroid
    if( !add_cell_centroids )
    {
        bool found_central_point( false );
        for( size_t fid = 0; fid<num_quad_faces; ++fid )
        {
            if( quad_face_gn[fid] != NULL )
            {
                cgn  = quad_face_gn[fid];
                found_central_point = true;
                break;
            }
        }
        if( !found_central_point )
        {
            for( size_t fid = 0; fid<num_tri_faces; ++fid )
            {
                if( tri_face_gn[fid] != NULL )
                {
                    cgn  = tri_face_gn[fid];
                    found_central_point = true;
                    break;
                }
            }
            if( !found_central_point )
                return false;
        }
    }

    /// 4.0 process cell with existing face centroids
    if( add_face_centroids )
    {
        /// considering cells with face centroids and possible edge centroids
        /// distinguish cases between:
        /// volumetric and planar elements ( num_faces > 1 || num_faces == 1 )
        /// elements with cell centroid and without ( cgn != face_cgn || cgn == face_cgn )
        std::set<size_t> edge;
        std::map<std::set<size_t>,GridNode*> edge_nodes;
        std::pair<typename std::map<std::set<size_t>,GridNode*>::iterator,bool> eit;
        bool edge_centroid_exist( false );
        bool planar_cell_exist( false );
        for( size_t fid = 0; fid<num_faces; ++fid )
        {
            const size_t num_face_nodes( face_nodes[fid].size() );
            if( (num_face_nodes == 4) || (num_face_nodes == 3) )
            {
                for( size_t nid = 0; nid < num_face_nodes; ++nid )
                {
                    edge.clear();
                    edge.insert( nodes[face_nodes[fid][nid]]->GetIdx() );
                    edge.insert( nodes[face_nodes[fid][(nid+1)%num_face_nodes]]->GetIdx() );
                    edge_centroid_exist = false;
                    if( edge_nodes.find( edge ) != edge_nodes.end() )
                    {
                        edge_centroid_exist = true;
                    }
                    else if( doesPointExist( additional_points,gn,
                                           nodes[face_nodes[fid][nid]],
                                           nodes[face_nodes[fid][(nid+1)%num_face_nodes]] ) )
                    {
                        edge_centroid_exist = true;
                    }
                    else if( add_edge_centroids )
                    {
                        addFaceCentroid( pgm,additional_points,gn,
                                         nodes[face_nodes[fid][nid]],
                                         nodes[face_nodes[fid][(nid+1)%num_face_nodes]] );
                        edge_centroid_exist = true;
                    }
                    if( edge_centroid_exist )
                    {
                        /// edge with centroid
                        eit = edge_nodes.insert( std::make_pair( edge, gn ) );
                        if( eit.second == true )
                            addExtraNode( extra_nodes, gn );
                        GridNode* en( (*eit.first).second );

                        /// volumetric element with cell centroid
                        if( is_volumetric_element )
                        {
                            if( face_cgn[fid]->GetPoint() != cgn->GetPoint() )
                            {
                                /// tetra0 cell
                                addTetraCell( cells, nodes[face_nodes[fid][nid]],en,face_cgn[fid],cgn );
                                /// tetra1 cell
                                addTetraCell( cells, nodes[face_nodes[fid][(nid+1)%num_face_nodes]],en,face_cgn[fid],cgn );
                            }
                        }
                        /// planar element
                        else if( !planar_cell_exist )
                        {
                            /// tri0 cell
                            addTriCell( cells, nodes[face_nodes[fid][nid]],en,face_cgn[fid] );
                            /// tri1 cell
                            addTriCell( cells, nodes[face_nodes[fid][(nid+1)%num_face_nodes]],en,face_cgn[fid] );
                        }
                        /// tri0 face
                        addTriFace(faces,fid, nodes[face_nodes[fid][nid]],en,face_cgn[fid] );
                        /// tri1 face
                        addTriFace(faces,fid, nodes[face_nodes[fid][(nid+1)%num_face_nodes]],en,face_cgn[fid] );
                    }
                    else
                    {
                        /// edge without centroid

                        /// volumetric element with cell centroid
                        if( is_volumetric_element )
                        {
                            if( face_cgn[fid]->GetPoint() != cgn->GetPoint() )
                                /// tetra cell
                                addTetraCell( cells, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],face_cgn[fid],cgn );
                        }
                        /// planar element
                        else if( !planar_cell_exist )
                        {
                            /// tri cell
                            addTriCell( cells, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],cgn );
                        }
                        /// tri face
                        addTriFace( faces,fid, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],face_cgn[fid] );
                    }
                }
                if( !is_volumetric_element )
                    planar_cell_exist = true;
            }
            else if( ( num_face_nodes == 2 ) && !is_volumetric_element )
            {
                for( size_t nid = 0; nid < num_face_nodes; ++nid )
                {
                    edge.clear();
                    edge.insert( nodes[face_nodes[fid][nid]]->GetIdx() );
                    edge.insert( nodes[face_nodes[fid][(nid+1)%num_face_nodes]]->GetIdx() );
                    edge_centroid_exist = false;
                    if( edge_nodes.find( edge ) != edge_nodes.end() )
                    {
                        edge_centroid_exist = true;
                    }
                    else if( doesPointExist( additional_points,gn,
                                           nodes[face_nodes[fid][nid]],
                                           nodes[face_nodes[fid][(nid+1)%num_face_nodes]] ) )
                    {
                        edge_centroid_exist = true;
                    }
                    else if( add_edge_centroids )
                    {
                        addFaceCentroid( pgm,additional_points,gn,
                                         nodes[face_nodes[fid][nid]],
                                         nodes[face_nodes[fid][(nid+1)%num_face_nodes]] );
                        edge_centroid_exist = true;
                    }
                    if( edge_centroid_exist )
                    {
                        /// edge with centroid
                        eit = edge_nodes.insert( std::make_pair( edge, gn ) );
                        if( eit.second == true )
                            addExtraNode( extra_nodes, gn );
                        GridNode* en( (*eit.first).second );

                        /// bar0 face
                        addBarFace(faces,fid, nodes[face_nodes[fid][nid]],en );
                        /// bar1 face
                        addBarFace(faces,fid, nodes[face_nodes[fid][(nid+1)%num_face_nodes]],en );
                    }
                    else
                    {
                        /// bar face
                        addBarFace( faces,fid, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]] );
                    }
                }
            }
        }
    }
    /// 5.0 process cell without all face centroids
    else
    {
        /// considering cells with only possible face centroids
        /// distinguish cases between:
        /// volumetric and planar elements ( num_faces > 1 || num_faces == 1 )
        /// elements with cell centroid and without ( cgn != face_cgn || cgn == face_cgn )
        size_t valid_edge;
        std::vector<size_t> diag_nodes;
        bool planar_cell_exist( false );
        for( size_t fid = 0; fid < num_faces; ++fid )
        {
            const size_t num_face_nodes( face_nodes[fid].size() );
            if( (num_face_nodes == 4) || (num_face_nodes == 3) )
            {
                /// face with centroid
                if( face_cgn[fid] != NULL )
                {
                    /// volumetric element with cell centroid
                    if( is_volumetric_element )
                    {
                        if( face_cgn[fid]->GetPoint() != cgn->GetPoint() )
                        {
                            for( size_t nid = 0; nid < num_face_nodes; ++nid )
                            {
                                /// tetra cell
                                addTetraCell( cells, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],face_cgn[fid],cgn );
                            }
                        }
                    }
                    /// planar element
                    else if( !planar_cell_exist )
                    {
                        for( size_t nid = 0; nid < num_face_nodes; ++nid )
                        {
                            /// tri cell
                            addTriCell( cells, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],face_cgn[fid] );
                        }
                    }
                    for( size_t nid = 0; nid < num_face_nodes; ++nid )
                    {
                        /// tri face
                        addTriFace( faces,fid, nodes[face_nodes[fid][nid]],nodes[face_nodes[fid][(nid+1)%num_face_nodes]],face_cgn[fid] );
                    }
                }
                /// quadrilateral face without centroid
                else if( num_face_nodes == 4 )
                {
                    valid_edge = isEdgeExist( additional_edges, diag_nodes,
                                              nodes[face_nodes[fid][0]],
                                              nodes[face_nodes[fid][1]],
                                              nodes[face_nodes[fid][2]],
                                              nodes[face_nodes[fid][3]] );
                    if( valid_edge == 3 )
                    {
                        csmp_error.notice( csmp::FATAL_ERROR,
                                           "processCell()",
                                           "Two diagonals exist!!! Unpredictable case!!!" );
                        return false;
                    }
                    else if( valid_edge == 0 && !tetra_mesh )
                    {
                        /// volumetric element with cell centroid
                        if( is_volumetric_element )
                        {
                            /// pyramid cell
                            addPyramidCell( cells, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]],cgn );
                        }
                        /// planar element
                        else if( !planar_cell_exist )
                        {
                            /// quad cell
                            addQuadCell( cells, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]] );
                        }
                        /// quad face
                        addQuadFace( faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]] );
                    }
                    else
                    {
                        if( valid_edge == 0 )
                        {
                            valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes,
                                                                                    nodes[face_nodes[fid][0]],
                                                                                    nodes[face_nodes[fid][1]],
                                                                                    nodes[face_nodes[fid][2]],
                                                                                    nodes[face_nodes[fid][3]] );
                            if( valid_edge == 3 )
                            {
                                csmp_error.notice( csmp::FATAL_ERROR,
                                                   "processCell()",
                                                   "Two diagonals exist!!! Unpredictable case!!!" );
                                return false;
                            }
                            else if( valid_edge == 0 )
                            {
                                csmp_error.notice( csmp::FATAL_ERROR,
                                                   "processCell()",
                                                   "Both diagonals are invalid!!! Unpredictable case!!!" );
                                return false;
                            }
                        }
                        /// volumetric element with cell centroid
                        if( is_volumetric_element )
                        {
                            /// tetra0 cell
                            addTetraCell( cells, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[2]]],cgn );
                            /// tetra1 cell
                            addTetraCell( cells, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[3]]],cgn );
                        }
                        /// planar element
                        else if( !planar_cell_exist )
                        {
                            /// tri0 cell
                            addTriCell( cells, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[2]]] );
                            /// tri1 cell
                            addTriCell( cells, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[3]]] );
                        }
                        /// tri0 face
                        addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[2]]] );
                        /// tri1 face
                        addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[0]]],nodes[face_nodes[fid][diag_nodes[1]]],nodes[face_nodes[fid][diag_nodes[3]]] );
                    }
                }
                /// triangular face without centroid
                else if( num_face_nodes == 3 )
                {
                    /// volumetric element with cell centroid
                    if( is_volumetric_element )
                    {
                        /// tetra cell
                        addTetraCell( cells, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],cgn );
                    }
                    /// planar element
                    else if( !planar_cell_exist )
                    {
                        /// tri cell
                        addTriCell( cells, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
                    }
                    /// tri face
                    addTriFace( faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
                }
                if( !is_volumetric_element )
                    planar_cell_exist = true;
            }
            else if( ( num_face_nodes == 2 ) && !is_volumetric_element )
            {
                /// bar face
                addBarFace( faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]] );
            }
        }
    }
    return true;
}










/**
    process degenerate hexahedron ( expected node and face numbering )

     5 x_____________x 4
       |\            |
       | \           |
       |  \          |
       | 6 x         |
       |   |         |
       x___|_________x
     2  \  |        1 \
         \ |           \
          \|            \
           x_____________\x
            3              0

Faces:  bottom quad face ( 0123 )
        right quad face  ( 1254 )
        left quad face   ( 2365 )
        top quad face    ( 0654 )
        right tri face   ( 014  )
        left tri face    ( 063  )
*/

bool processDegenerateHexahedronElement( bool tetra_mesh,
                                         PolygonGridManager* pgm,
                                         std::map<std::set<size_t>,GridNode*>& additional_points,
                                         std::set<std::set<size_t> >& additional_edges,
                                         std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                                         std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                                         std::vector<GridNode*>& extra_nodes,
                                         const std::vector<GridNode*>& nodes,
                                         const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_quad_faces( 4 );
    const size_t cell_faces( 6 );
    const size_t cell_nodes( 7 );

    /// build up topology: face and edge nodes
    size_t num_diags( 0 );
    std::vector<std::vector<size_t> > diag_nodes( cell_faces );
    std::vector<size_t>               point_diags(cell_nodes,0);

    /// count number of actually existing diagonals and
    /// assign to each vertex number of exisiting diagonals
    size_t valid_edge( 0 );
    num_diags = 0;
    for( size_t fid = 0; fid < num_quad_faces; ++fid )
    {
        valid_edge = isEdgeExist( additional_edges, diag_nodes[fid],
                                  nodes[face_nodes[fid][0]],
                                  nodes[face_nodes[fid][1]],
                                  nodes[face_nodes[fid][2]],
                                  nodes[face_nodes[fid][3]] );
        if( valid_edge == 3 )
        {
            csmp_error.notice( csmp::FATAL_ERROR,
                               "processDegenerateHexahedronElement()",
                               "Two diagonals exist!!! Unpredictable case!!!" );
            return false;
        }
        else if( valid_edge != 0 )
        {
            ++num_diags;
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ]);
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ]);
        }
    }

    if( num_diags < num_quad_faces )
    {
        /// inserting remaining diagonals
        for( size_t fid = 0; fid < num_quad_faces; ++fid )
        {
            if( diag_nodes[fid].empty() )
            {
                valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes[fid],
                                                                        nodes[face_nodes[fid][0]],
                                                                        nodes[face_nodes[fid][1]],
                                                                        nodes[face_nodes[fid][2]],
                                                                        nodes[face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processDegenerateHexahedronElement()",
                                       "Two diagonals exist!!! Unpredictable case!!!" );
                    return false;
                }
                else if( valid_edge == 0 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processDegenerateHexahedronElement()",
                                       "Both diagonals are invalid!!! Unpredictable case!!!" );
                    return false;
                }
                else
                {
                    ++num_diags;
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ]);
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ]);
                }
            }
        }
    }

    /// define if we deal with undividable case
    /// check if topology is correct
    size_t sum( 0 );
    for( size_t nid = 0; nid < cell_nodes; ++nid )
        sum += point_diags[nid];
    assert( sum == 2*num_quad_faces );
    assert( num_diags == num_quad_faces );

    /// dividable cases are:
    /// symmetric diagonals on bottom and top faces: (0,2),(0,5) or (1,3),(4,6)
    /// 3 diagonals comming though nodes 2 or 5
    bool undividable( true );
    size_t dividable_case( 0 );
    if( ( face_nodes[0][diag_nodes[0][0]] == 0 || face_nodes[0][diag_nodes[0][1]] == 0 ) && ( face_nodes[3][diag_nodes[3][0]] == 0 || face_nodes[3][diag_nodes[3][1]] == 0 ) )
    {
        dividable_case = 0;
        undividable = false;
    }
    else if( ( face_nodes[0][diag_nodes[0][0]] == 1 || face_nodes[0][diag_nodes[0][1]] == 1 ) && ( face_nodes[3][diag_nodes[3][0]] == 4 || face_nodes[3][diag_nodes[3][1]] == 4 ) )
    {
        dividable_case = 1;
        undividable = false;
    }
    else if( point_diags[2] == 3 || point_diags[5] == 3 )
    {
        dividable_case = 2;
        undividable = false;
    }
    else
        undividable = true;

    if( !undividable )
    {
        /// we are considering 3 configurations

        /// 1) diagonals 02 and 05 ( 4 tetrahedra )
        if( dividable_case == 0 )
        {
            /// tetrahedron0
            addTetraCell(cells,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][2]]],nodes[0]);

            /// tetrahedron1
            addTetraCell(cells,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][3]]],nodes[0]);

            /// tetrahedron2
            addTetraCell(cells,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][2]]],nodes[0]);

            /// tetrahedron3
            addTetraCell(cells,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][3]]],nodes[0]);
        }
        /// 2) diagonals (1,3) and (4,6) (5 tetrahedra )
        else if( dividable_case == 1 )
        {
            std::vector<size_t> nodemap;
            if( face_nodes[1][diag_nodes[1][0]] == 1 || face_nodes[1][diag_nodes[1][1]] == 1 )
            {
                nodemap.push_back(1);
                nodemap.push_back(6);
                nodemap.push_back(4);
                nodemap.push_back(3);
                nodemap.push_back(5);
            }
            else
            {
                nodemap.push_back(4);
                nodemap.push_back(3);
                nodemap.push_back(1);
                nodemap.push_back(6);
                nodemap.push_back(2);
            }

            /// tetrahedron0
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[0]);

            /// tetrahedron1
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[3]],nodes[0]);

            /// tetrahedron2
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[4]]);

            /// tetrahedron3
            addTetraCell(cells,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][2]]],nodes[nodemap[0]]);

            /// tetrahedron4
            addTetraCell(cells,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][3]]],nodes[nodemap[0]]);
        }
        /// 3) diagonals 02 and 46 or 05 and 13 (4 tetrahedra )
        else if( dividable_case == 3 )
        {
            std::vector<size_t> nodemap;
            if( face_nodes[0][diag_nodes[0][0]] == 0 || face_nodes[0][diag_nodes[0][1]] == 0 )
            {
                nodemap.push_back(0);
                nodemap.push_back(2);
                nodemap.push_back(4);
                nodemap.push_back(6);
                nodemap.push_back(1);
                nodemap.push_back(3);
            }
            else
            {
                nodemap.push_back(0);
                nodemap.push_back(5);
                nodemap.push_back(1);
                nodemap.push_back(3);
                nodemap.push_back(4);
                nodemap.push_back(6);
            }

            /// tetrahedron0
            addTetraCell(cells,nodes[2],nodes[5],nodes[nodemap[2]],nodes[nodemap[3]]);

            /// tetrahedron1
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[3]]);

            /// tetrahedron2
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[4]]);

            /// tetrahedron3
            addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[3]],nodes[nodemap[5]]);
        }

        /// add faces
        /// face (0,1,2,3)
        addTriFace(faces,0,nodes[face_nodes[0][diag_nodes[0][0]]],nodes[face_nodes[0][diag_nodes[0][1]]],nodes[face_nodes[0][diag_nodes[0][2]]]);
        addTriFace(faces,0,nodes[face_nodes[0][diag_nodes[0][0]]],nodes[face_nodes[0][diag_nodes[0][1]]],nodes[face_nodes[0][diag_nodes[0][3]]]);
        /// face (1,2,5,4)
        addTriFace(faces,1,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][2]]]);
        addTriFace(faces,1,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][3]]]);
        /// face (2,3,6,5)
        addTriFace(faces,2,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][2]]]);
        addTriFace(faces,2,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][3]]]);
        /// face (0,4,5,6)
        addTriFace(faces,3,nodes[face_nodes[3][diag_nodes[3][0]]],nodes[face_nodes[3][diag_nodes[3][1]]],nodes[face_nodes[3][diag_nodes[3][2]]]);
        addTriFace(faces,3,nodes[face_nodes[3][diag_nodes[3][0]]],nodes[face_nodes[3][diag_nodes[3][1]]],nodes[face_nodes[3][diag_nodes[3][3]]]);

        return true;
    }

    /// in case if one deals with undividable case
    if( undividable )
    {
        /// add cell centroid
        GridNode* cgn;
        addCellCentroid(pgm,additional_points,additional_edges,cgn,nodes,face_nodes);
        addExtraNode( extra_nodes, cgn );

        for( size_t fid = 0; fid<num_quad_faces; ++fid )
        {
            assert( !diag_nodes.empty() );

            /// tetrahedron0
            addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]],cgn);

            /// tetrahedron1
            addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]],cgn);

            /// add faces
            addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]]);
            addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]]);
        }
        /// tetrahedron ( 0,1,4,7 )
        addTetraCell(cells,nodes[0],nodes[1],nodes[4],cgn);
        /// tetrahedron ( 0,3,6,7 )
        addTetraCell(cells,nodes[0],nodes[3],nodes[6],cgn);
        /// triangle ( 0,1,4 )
        addTriFace(faces,4, nodes[0],nodes[1],nodes[4]);
        /// triangle ( 0,3,6 )
        addTriFace(faces,5, nodes[0],nodes[3],nodes[6]);

        return true;
    }
    return false;
}











/**
    process hexahedron ( expected node and face numbering )

     7 x_____________x 6
       |\            |\
       | \           | \
       |  \          |  \
       | 4 x_________|___x 5
       |   |         |   |
       x___|_________x   |
      3 \  |        2 \  |
         \ |           \ |
          \|            \|
           x_____________x
            0              1

Faces:  bottom quad face ( 0321 )
        left quad face   ( 0473 )
        front quad face  ( 0154 )
        right quad face  ( 1265 )
        back quad face   ( 2376 )
        top quad face    ( 4567 )

certain pieces of algorithm are based on:
[ 1 ] J.Dompierre,P.Labbe,M-G. Vallet,R.Camarereo. How to subdivide pyramids,prisms aand hexahedra into tetrahedra.
[ 2 ] G.Albertelli,R.A.Crawfis. Efficient subdivision of finite-element dataset into consistent tetrahedra.

*/

bool processHexahedronElement( bool tetra_mesh,
                              PolygonGridManager* pgm,
                              std::map<std::set<size_t>,GridNode*>& additional_points,
                              std::set<std::set<size_t> >& additional_edges,
                              std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                              std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                              std::vector<GridNode*>& extra_nodes,
                              const std::vector<GridNode*>& nodes,
                              const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_quad_faces( 6U );
    const size_t num_hexa_fem_nodes( 8U );

    /// build up topology: face and edge nodes
    size_t num_diags( 0 );
    std::vector<std::vector<size_t> > diag_nodes(  num_quad_faces );
    std::vector<size_t>               point_diags( num_hexa_fem_nodes,0);

    /// count number of actually existing diagonals and
    /// assign to each vertex number of exisiting diagonals
    size_t valid_edge( 0 );
    num_diags = 0;
    for( size_t fid = 0; fid < num_quad_faces; ++fid )
    {
        valid_edge = isEdgeExist( additional_edges, diag_nodes[fid],
                                  nodes[face_nodes[fid][0]],
                                  nodes[face_nodes[fid][1]],
                                  nodes[face_nodes[fid][2]],
                                  nodes[face_nodes[fid][3]] );
        if( valid_edge == 3 )
        {
            csmp_error.notice( csmp::FATAL_ERROR,
                               "processHexahedronElement()",
                               "Two diagonals exist!!! Unpredictable case!!!" );
            return false;
        }
        else if( valid_edge != 0 )
        {
            ++num_diags;
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ]);
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ]);
        }
    }

    if( num_diags == 0 && !tetra_mesh )
    {
        /// hexahedron
        addHexaCell(cells,nodes[0],nodes[1],nodes[2],nodes[3],nodes[4],nodes[5],nodes[6],nodes[7]);
        /// add faces
        for( size_t fid = 0; fid<num_quad_faces; ++fid )
        {
            if( face_nodes[fid].size() == 4 )
                addQuadFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]] );
        }
        return true;
    }
    else if( ( num_diags < num_quad_faces ) && tetra_mesh )
    {
        /// inserting remaining diagonals
        for( size_t fid = 0; fid < num_quad_faces; ++fid )
        {
            if( diag_nodes[fid].empty() )
            {
                valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes[fid],
                                                                        nodes[face_nodes[fid][0]],
                                                                        nodes[face_nodes[fid][1]],
                                                                        nodes[face_nodes[fid][2]],
                                                                        nodes[face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processHexahedronElement()",
                                       "Two diagonals exist!!! Unpredictable case!!!" );
                    return false;
                }
                else if( valid_edge == 0 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processHexahedronElement()",
                                       "Both diagonals are invalid!!! Unpredictable case!!!" );
                    return false;
                }
                else
                {
                    ++num_diags;
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ] );
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ] );
                }
            }
        }
    }

    bool undividable( true );
    if( tetra_mesh )
    {
        /// check if topology is correct
        size_t sum( 0 );
        for( size_t nid = 0; nid < num_hexa_fem_nodes; ++nid )
            sum += point_diags[nid];
        assert( sum == 2*num_quad_faces );
        assert( num_diags == num_quad_faces );

        /// dividable cases are:
        /// 3 diagonals comming trough one of the points
        /// undividable: 2,2,2,2,1,1,1,1
        std::vector<size_t> three_diag_nodes;
        for( size_t nid = 0; nid < num_hexa_fem_nodes; ++nid )
            if( point_diags[nid] == 3U )
                three_diag_nodes.push_back( nid );
        if( three_diag_nodes.empty() )
            undividable = true;
        else
            undividable = false;

        if( !undividable )
        {
            /// put three diagonal node at 0 position
            std::vector<size_t> nodemap0;
            for( size_t nid = 0; nid < num_hexa_fem_nodes; ++nid )
                nodemap0.push_back( nid );
            std::vector<size_t> facemap0;
            for( size_t fid = 0; fid < num_quad_faces; ++fid )
                facemap0.push_back( fid );
            switch( three_diag_nodes[0] )
            {
            case 1U:
                facemap0[0] = 2;
                facemap0[1] = 3;
                facemap0[2] = 0;
                facemap0[3] = 1;
                facemap0[4] = 5;
                facemap0[5] = 4;

                nodemap0[0] = 1;
                nodemap0[1] = 0;
                nodemap0[2] = 4;
                nodemap0[3] = 5;
                nodemap0[4] = 2;
                nodemap0[5] = 3;
                nodemap0[6] = 7;
                nodemap0[7] = 6;
                break;
            case 2U:
                facemap0[0] = 3;
                facemap0[1] = 4;
                facemap0[2] = 0;
                facemap0[3] = 2;
                facemap0[4] = 5;
                facemap0[5] = 1;

                nodemap0[0] = 2;
                nodemap0[1] = 1;
                nodemap0[2] = 5;
                nodemap0[3] = 6;
                nodemap0[4] = 3;
                nodemap0[5] = 0;
                nodemap0[6] = 4;
                nodemap0[7] = 7;
                break;
            case 3U:
                facemap0[0] = 0;
                facemap0[1] = 4;
                facemap0[2] = 1;
                facemap0[3] = 2;
                facemap0[4] = 3;
                facemap0[5] = 5;

                nodemap0[0] = 3;
                nodemap0[1] = 0;
                nodemap0[2] = 1;
                nodemap0[3] = 2;
                nodemap0[4] = 7;
                nodemap0[5] = 4;
                nodemap0[6] = 5;
                nodemap0[7] = 6;
                break;
            case 4U:
                facemap0[0] = 1;
                facemap0[1] = 5;
                facemap0[2] = 2;
                facemap0[3] = 0;
                facemap0[4] = 4;
                facemap0[5] = 3;

                nodemap0[0] = 4;
                nodemap0[1] = 0;
                nodemap0[2] = 3;
                nodemap0[3] = 7;
                nodemap0[4] = 5;
                nodemap0[5] = 1;
                nodemap0[6] = 2;
                nodemap0[7] = 6;
                break;
            case 5U:
                facemap0[0] = 2;
                facemap0[1] = 5;
                facemap0[2] = 3;
                facemap0[3] = 0;
                facemap0[4] = 1;
                facemap0[5] = 4;

                nodemap0[0] = 5;
                nodemap0[1] = 1;
                nodemap0[2] = 0;
                nodemap0[3] = 4;
                nodemap0[4] = 6;
                nodemap0[5] = 2;
                nodemap0[6] = 3;
                nodemap0[7] = 7;
                break;
            case 6U:
                facemap0[0] = 3;
                facemap0[1] = 5;
                facemap0[2] = 4;
                facemap0[3] = 0;
                facemap0[4] = 2;
                facemap0[5] = 1;

                nodemap0[0] = 6;
                nodemap0[1] = 2;
                nodemap0[2] = 1;
                nodemap0[3] = 5;
                nodemap0[4] = 7;
                nodemap0[5] = 3;
                nodemap0[6] = 0;
                nodemap0[7] = 4;
                break;
            case 7U:
                facemap0[0] = 4;
                facemap0[1] = 5;
                facemap0[2] = 1;
                facemap0[3] = 0;
                facemap0[4] = 3;
                facemap0[5] = 2;

                nodemap0[0] = 7;
                nodemap0[1] = 3;
                nodemap0[2] = 2;
                nodemap0[3] = 6;
                nodemap0[4] = 4;
                nodemap0[5] = 0;
                nodemap0[6] = 1;
                nodemap0[7] = 5;
                break;
            default:
                break;
            }

            /// in order to get canonical numbering lets perform suitable rotation.
            /// in case of vertex 6 having 1 diagonal we want face with this diagonal on the right
            /// in case of vertex 6 having 2 diagonals we want face adjacent to those 2 on the rigth as well
            /// in case of vertex 6 hvaing 0 or 3 diagonals we keep the positions of faces
            std::vector<size_t> fcode(3,0);
            for( size_t eid = 0; eid < 3U; ++eid )
                fcode[ eid ] = ( ( face_nodes[ facemap0[eid+3] ][diag_nodes[ facemap0[eid+3] ][0]] == nodemap0[6] || face_nodes[ facemap0[eid+3] ][diag_nodes[ facemap0[eid+3] ][1]] == nodemap0[6] ) ? 1 : 0 );

            std::vector<size_t> nodemap;
            for( size_t nid = 0; nid < num_hexa_fem_nodes; ++nid )
                nodemap.push_back( nodemap0[nid] );
            std::vector<size_t> facemap;
            for( size_t fid = 0; fid < num_quad_faces; ++fid )
                facemap.push_back( facemap0[fid] );
            if(    ( fcode[0] == 0 && fcode[1] == 0 && fcode[2] == 1 )
                || ( fcode[0] == 1 && fcode[1] == 1 && fcode[2] == 0 ) )
            {
                /// 120 rotation over the diagonal 06
                facemap[0] = facemap0[2];
                facemap[1] = facemap0[0];
                facemap[2] = facemap0[1];
                facemap[3] = facemap0[5];
                facemap[4] = facemap0[3];
                facemap[5] = facemap0[4];
                nodemap[0] = nodemap0[0];
                nodemap[1] = nodemap0[4];
                nodemap[2] = nodemap0[5];
                nodemap[3] = nodemap0[1];
                nodemap[4] = nodemap0[3];
                nodemap[5] = nodemap0[7];
                nodemap[6] = nodemap0[6];
                nodemap[7] = nodemap0[2];
            }
            else if(    ( fcode[0] == 0 && fcode[1] == 1 && fcode[2] == 0 )
                     || ( fcode[0] == 1 && fcode[1] == 0 && fcode[2] == 1 ) )
            {
                /// 240 rotation over the diagonal 06
                facemap[0] = facemap0[1];
                facemap[1] = facemap0[2];
                facemap[2] = facemap0[0];
                facemap[3] = facemap0[4];
                facemap[4] = facemap0[5];
                facemap[5] = facemap0[3];
                nodemap[0] = nodemap0[0];
                nodemap[1] = nodemap0[3];
                nodemap[2] = nodemap0[7];
                nodemap[3] = nodemap0[4];
                nodemap[4] = nodemap0[1];
                nodemap[5] = nodemap0[2];
                nodemap[6] = nodemap0[6];
                nodemap[7] = nodemap0[5];
            }
            //facemap0.clear();
            //nodemap0.clear();

            /// After all the necessary rotations we performed we consider 4 configurations

            /// 1) no diagonal that goes through vertex 6 ( 5 tetrahedra )
            if( point_diags[ nodemap[ 6 ] ] == 0 )
            {
                assert( point_diags[ nodemap[ 0 ] ] == 3 );
                assert( point_diags[ nodemap[ 1 ] ] == 0 );
                assert( point_diags[ nodemap[ 2 ] ] == 3 );
                assert( point_diags[ nodemap[ 3 ] ] == 0 );
                assert( point_diags[ nodemap[ 4 ] ] == 0 );
                assert( point_diags[ nodemap[ 5 ] ] == 3 );
                assert( point_diags[ nodemap[ 6 ] ] == 0 );
                assert( point_diags[ nodemap[ 7 ] ] == 3 );

                /// tetrahedron0 ( 0,1,2,5 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[5]]);

                /// tetrahedron1 ( 0,2,7,5 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[5]]);

                /// tetrahedron2 ( 0,2,3,7 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]],nodes[nodemap[7]]);

                /// tetrahedron3 ( 0,5,7,4 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[4]]);

                /// tetrahedron4 ( 2,7,5,6 )
                addTetraCell(cells,nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[5]],nodes[nodemap[6]]);

                /// add faces
                /// bottom face ( 0321 )
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]]);
                /// left face   ( 3047 )
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]]);
                /// front face  ( 0154 )
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[4]]);
                /// right face  ( 1265 )
                addTriFace(faces,facemap[3], nodes[nodemap[2]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[3], nodes[nodemap[2]],nodes[nodemap[5]],nodes[nodemap[6]]);
                /// back face   ( 2376 )
                addTriFace(faces,facemap[4], nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[4], nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[6]]);
                /// top face    ( 4567 )
                addTriFace(faces,facemap[5], nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[4]]);
                addTriFace(faces,facemap[5], nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[6]]);
            }
            /// 2) one diagonal that goes through vertex 6 ( 6 tetrahedra )
            else if( point_diags[ nodemap[ 6 ] ] == 1 )
            {
                assert( point_diags[ nodemap[ 0 ] ] == 3 );
                assert( point_diags[ nodemap[ 1 ] ] == 1 );
                assert( point_diags[ nodemap[ 2 ] ] == 2 );
                assert( point_diags[ nodemap[ 3 ] ] == 0 );
                assert( point_diags[ nodemap[ 4 ] ] == 0 );
                assert( point_diags[ nodemap[ 5 ] ] == 2 );
                assert( point_diags[ nodemap[ 6 ] ] == 1 );
                assert( point_diags[ nodemap[ 7 ] ] == 3 );

                /// tetrahedron0 ( 0,5,7,4 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[4]]);

                /// tetrahedron1 ( 0,1,7,5 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[7]],nodes[nodemap[5]]);

                /// tetrahedron2 ( 1,6,7,5 )
                addTetraCell(cells,nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[7]],nodes[nodemap[5]]);

                /// tetrahedron ( 0,7,2,3 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[2]],nodes[nodemap[3]]);

                /// tetrahedron3 ( 0,7,1,2 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[1]],nodes[nodemap[2]]);

                /// tetrahedron5 ( 1,7,6,2 )
                addTetraCell(cells,nodes[nodemap[1]],nodes[nodemap[7]],nodes[nodemap[6]],nodes[nodemap[2]]);

                /// add faces
                /// bottom face ( 0321 )
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]]);
                /// left face   ( 3047 )
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]]);
                /// front face  ( 0154 )
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[4]]);
                /// right face  ( 1265 )
                addTriFace(faces,facemap[3], nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[2]]);
                addTriFace(faces,facemap[3], nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[5]]);
                /// back face   ( 2376 )
                addTriFace(faces,facemap[4], nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[4], nodes[nodemap[2]],nodes[nodemap[7]],nodes[nodemap[6]]);
                /// top face    ( 4567 )
                addTriFace(faces,facemap[5], nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[4]]);
                addTriFace(faces,facemap[5], nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[6]]);
            }
            /// 3) two diagonals that goes through vertex 6 ( 6 tetrahedra )
            else if( point_diags[ nodemap[ 6 ] ] == 2 )
            {
                assert( point_diags[ nodemap[ 0 ] ] == 3 );
                assert( point_diags[ nodemap[ 1 ] ] == 0 );
                assert( point_diags[ nodemap[ 2 ] ] == 2 );
                assert( point_diags[ nodemap[ 3 ] ] == 1 );
                assert( point_diags[ nodemap[ 4 ] ] == 1 );
                assert( point_diags[ nodemap[ 5 ] ] == 2 );
                assert( point_diags[ nodemap[ 6 ] ] == 2 );
                assert( point_diags[ nodemap[ 7 ] ] == 1 );

                /// tetrahedron0 ( 0,4,5,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[4]],nodes[nodemap[5]],nodes[nodemap[6]]);

                /// tetrahedron1 ( 0,3,7,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[3]],nodes[nodemap[7]],nodes[nodemap[6]]);

                /// tetrahedron2 ( 0,7,4,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]],nodes[nodemap[6]]);

                /// tetrahedron ( 0,1,2,5 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[5]]);

                /// tetrahedron3 ( 0,3,6,2 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[3]],nodes[nodemap[6]],nodes[nodemap[2]]);

                /// tetrahedron5 ( 0,6,5,2 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[6]],nodes[nodemap[5]],nodes[nodemap[2]]);

                /// add faces
                /// bottom face ( 0321 )
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]]);
                /// left face   ( 3047 )
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]]);
                /// front face  ( 0154 )
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[4]]);
                /// right face  ( 1265 )
                addTriFace(faces,facemap[3], nodes[nodemap[2]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[3], nodes[nodemap[2]],nodes[nodemap[5]],nodes[nodemap[6]]);
                /// back face   ( 2376 )
                addTriFace(faces,facemap[4], nodes[nodemap[3]],nodes[nodemap[6]],nodes[nodemap[2]]);
                addTriFace(faces,facemap[4], nodes[nodemap[3]],nodes[nodemap[6]],nodes[nodemap[7]]);
                /// top face    ( 4567 )
                addTriFace(faces,facemap[5], nodes[nodemap[4]],nodes[nodemap[6]],nodes[nodemap[5]]);
                addTriFace(faces,facemap[5], nodes[nodemap[4]],nodes[nodemap[6]],nodes[nodemap[7]]);
            }
            /// 4) three diagonals that goes through vertex 6 ( 6 tetrahedra )
            else if( point_diags[ nodemap[ 6 ] ] == 3 )
            {
                assert( point_diags[ nodemap[ 0 ] ] == 3 );
                assert( point_diags[ nodemap[ 1 ] ] == 1 );
                assert( point_diags[ nodemap[ 2 ] ] == 1 );
                assert( point_diags[ nodemap[ 3 ] ] == 1 );
                assert( point_diags[ nodemap[ 4 ] ] == 1 );
                assert( point_diags[ nodemap[ 5 ] ] == 1 );
                assert( point_diags[ nodemap[ 6 ] ] == 3 );
                assert( point_diags[ nodemap[ 7 ] ] == 1 );

                /// tetrahedron0 ( 0,2,3,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]],nodes[nodemap[6]]);

                /// tetrahedron1 ( 0,3,7,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[3]],nodes[nodemap[7]],nodes[nodemap[6]]);

                /// tetrahedron2 ( 0,7,4,6 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]],nodes[nodemap[6]]);

                /// tetrahedron ( 0,5,6,4 )
                addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[6]],nodes[nodemap[4]]);

                /// tetrahedron3 ( 1,5,6,0 )
                addTetraCell(cells,nodes[nodemap[1]],nodes[nodemap[5]],nodes[nodemap[6]],nodes[nodemap[0]]);

                /// tetrahedron5 ( 1,6,2,0 )
                addTetraCell(cells,nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[2]],nodes[nodemap[0]]);

                /// add faces
                /// bottom face ( 0321 )
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[0], nodes[nodemap[0]],nodes[nodemap[2]],nodes[nodemap[3]]);
                /// left face   ( 3047 )
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[3]]);
                addTriFace(faces,facemap[1], nodes[nodemap[0]],nodes[nodemap[7]],nodes[nodemap[4]]);
                /// front face  ( 0154 )
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[1]]);
                addTriFace(faces,facemap[2], nodes[nodemap[0]],nodes[nodemap[5]],nodes[nodemap[4]]);
                /// right face  ( 1265 )
                addTriFace(faces,facemap[3], nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[2]]);
                addTriFace(faces,facemap[3], nodes[nodemap[1]],nodes[nodemap[6]],nodes[nodemap[5]]);
                /// back face   ( 2376 )
                addTriFace(faces,facemap[4], nodes[nodemap[3]],nodes[nodemap[6]],nodes[nodemap[2]]);
                addTriFace(faces,facemap[4], nodes[nodemap[3]],nodes[nodemap[6]],nodes[nodemap[7]]);
                /// top face    ( 4567 )
                addTriFace(faces,facemap[5], nodes[nodemap[4]],nodes[nodemap[6]],nodes[nodemap[5]]);
                addTriFace(faces,facemap[5], nodes[nodemap[4]],nodes[nodemap[6]],nodes[nodemap[7]]);
            }
            return true;
        }
    }

    /// in case if one deals with undividable case
    if( undividable )
    {
        /// add cell centroid
        GridNode* cgn;
        addCellCentroid(pgm,additional_points,additional_edges,cgn,nodes,face_nodes);
        addExtraNode( extra_nodes, cgn );

        for( size_t fid = 0; fid<num_quad_faces; ++fid )
        {
            if( !diag_nodes[ fid ].empty() )
            {
                /// tetrahedron0
                addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]],cgn );
                /// tetrahedron1
                addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]],cgn );

                /// add faces
                addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]] );
                addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]] );
            }
            else
            {
                /// pyramid
                addPyramidCell(cells,nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]],cgn );
            }
        }
        return true;
    }
    return false;
}













/**
  process prism ( expected node and face numbering )

                5
                x
               /|\
              / | \
             /  |  \
            /   |   \x 2
           /    |  / |
          /     | /  |
         /      |/   |
        /       /    |
       /       /|    |
      /       / |    |
   3 x_______/__x 4  |
      \     /    \   |
       \   /      \  |
        \ /        \ |
         x__________x|
        0            1

Faces:  left quad face    ( 0253 )
        bottom quad face  ( 0143 )
        right quad face   ( 1452 )
        front tri face    ( 012  )
        back tri face     ( 354  )
*/

bool processPrismElement( bool tetra_mesh,
                          PolygonGridManager* pgm,
                          std::map<std::set<size_t>,GridNode*>& additional_points,
                          std::set<std::set<size_t> >& additional_edges,
                          std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                          std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                          std::vector<GridNode*>& extra_nodes,
                          const std::vector<GridNode*>& nodes,
                          const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    const size_t num_quad_faces( 3 );
    const size_t cell_faces( 5 );
    const size_t cell_nodes( 6 );

    /// build up topology: face and edge nodes
    size_t num_diags( 0 );
    std::vector<std::vector<size_t> > diag_nodes( cell_faces );
    std::vector<size_t>               point_diags(cell_nodes,0);

    /// count number of actually existing diagonals and
    /// assign to each vertex number of exisiting diagonals
    size_t valid_edge( 0 );
    num_diags = 0;
    for( size_t fid = 0; fid < num_quad_faces; ++fid )
    {
        valid_edge = isEdgeExist( additional_edges, diag_nodes[fid],
                                  nodes[face_nodes[fid][0]],
                                  nodes[face_nodes[fid][1]],
                                  nodes[face_nodes[fid][2]],
                                  nodes[face_nodes[fid][3]] );
        if( valid_edge == 3 )
        {
            csmp_error.notice( csmp::FATAL_ERROR,
                               "processPrismElement()",
                               "Two diagonals exist!!! Unpredictable case!!!" );
            return false;
        }
        else if( valid_edge != 0 )
        {
            ++num_diags;
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ]);
            ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ]);
        }
    }

    if( num_diags == 0 && !tetra_mesh )
    {
        /// prism
        addPrismCell(cells,nodes[0],nodes[1],nodes[2],nodes[3],nodes[4],nodes[5]);
        /// add faces
        for( size_t fid = 0; fid<cell_faces; ++fid )
        {
            if( face_nodes[fid].size() == 4 )
                addQuadFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]],nodes[face_nodes[fid][3]] );
            else if( face_nodes[fid].size() == 4 )
                addTriFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
        }
        return true;
    }
    else if( num_diags < num_quad_faces )
    {
        /// inserting remaining diagonals
        for( size_t fid = 0; fid < num_quad_faces; ++fid )
        {
            if( diag_nodes[fid].empty() )
            {
                valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes[fid],
                                                                        nodes[face_nodes[fid][0]],
                                                                        nodes[face_nodes[fid][1]],
                                                                        nodes[face_nodes[fid][2]],
                                                                        nodes[face_nodes[fid][3]] );
                if( valid_edge == 3 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processPrismElement()",
                                       "Two diagonals exist!!! Unpredictable case!!!" );
                    return false;
                }
                else if( valid_edge == 0 )
                {
                    csmp_error.notice( csmp::FATAL_ERROR,
                                       "processPrismElement()",
                                       "Both diagonals are invalid!!! Unpredictable case!!!" );
                    return false;
                }
                else
                {
                    ++num_diags;
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][0]] ]);
                    ++(point_diags[ face_nodes[fid][diag_nodes[fid][1]] ]);
                }
            }
        }
    }

    /// check if topology is correct
    size_t sum( 0 );
    for( size_t nid = 0; nid < cell_nodes; ++nid )
        sum += point_diags[nid];
    assert( sum == 2*num_quad_faces );
    assert( num_diags == num_quad_faces );

    /// dividable cases are:
    /// 2 diagonals comming trough one of the points
    bool undividable( true );
    size_t two_diag_point(0);
    /// considering only nodes from first triangular face
    for( size_t nid = 0; nid < 3U; ++nid )
    {
        if( point_diags[nid] == 2 )
        {
            undividable = false;
            two_diag_point = nid;
        }
    }

    if( !undividable )
    {
        std::vector<size_t> nodemap;
        if( two_diag_point == 0 )
        {
            nodemap.push_back( 0 );
            nodemap.push_back( 3 );
            nodemap.push_back( 4 );
            nodemap.push_back( 5 );
            nodemap.push_back( face_nodes[2][diag_nodes[2][0]] );
            nodemap.push_back( face_nodes[2][diag_nodes[2][1]] );
            nodemap.push_back( face_nodes[2][diag_nodes[2][2]] );
            nodemap.push_back( face_nodes[2][diag_nodes[2][3]] );

            /// add faces
            /// face (0,2,5,3)
            addTriFace(faces,0,nodes[0],nodes[5],nodes[2]);
            addTriFace(faces,0,nodes[0],nodes[5],nodes[3]);
            /// face (0,1,4,3)
            addTriFace(faces,1,nodes[0],nodes[4],nodes[1]);
            addTriFace(faces,1,nodes[0],nodes[4],nodes[3]);
            /// face (1,2,5,4)
            addTriFace(faces,2,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][2]]]);
            addTriFace(faces,2,nodes[face_nodes[2][diag_nodes[2][0]]],nodes[face_nodes[2][diag_nodes[2][1]]],nodes[face_nodes[2][diag_nodes[2][3]]]);
        }
        else if( two_diag_point == 1 )
        {
            nodemap.push_back( 1 );
            nodemap.push_back( 4 );
            nodemap.push_back( 3 );
            nodemap.push_back( 5 );
            nodemap.push_back( face_nodes[0][diag_nodes[0][0]] );
            nodemap.push_back( face_nodes[0][diag_nodes[0][1]] );
            nodemap.push_back( face_nodes[0][diag_nodes[0][2]] );
            nodemap.push_back( face_nodes[0][diag_nodes[0][3]] );

            /// add faces
            /// face (0,1,4,3)
            addTriFace(faces,1,nodes[1],nodes[3],nodes[0]);
            addTriFace(faces,1,nodes[1],nodes[3],nodes[4]);
            /// face (1,2,5,4)
            addTriFace(faces,2,nodes[1],nodes[5],nodes[2]);
            addTriFace(faces,2,nodes[1],nodes[5],nodes[4]);
            /// face (0,2,5,3)
            addTriFace(faces,0,nodes[face_nodes[0][diag_nodes[0][0]]],nodes[face_nodes[0][diag_nodes[0][1]]],nodes[face_nodes[0][diag_nodes[0][2]]]);
            addTriFace(faces,0,nodes[face_nodes[0][diag_nodes[0][0]]],nodes[face_nodes[0][diag_nodes[0][1]]],nodes[face_nodes[0][diag_nodes[0][3]]]);
        }
        else
        {
            nodemap.push_back( 2 );
            nodemap.push_back( 5 );
            nodemap.push_back( 3 );
            nodemap.push_back( 4 );
            nodemap.push_back( face_nodes[1][diag_nodes[1][0]] );
            nodemap.push_back( face_nodes[1][diag_nodes[1][1]] );
            nodemap.push_back( face_nodes[1][diag_nodes[1][2]] );
            nodemap.push_back( face_nodes[1][diag_nodes[1][3]] );

            /// add faces
            /// face (0,2,5,3)
            addTriFace(faces,0,nodes[2],nodes[3],nodes[0]);
            addTriFace(faces,0,nodes[2],nodes[3],nodes[5]);
            /// face (1,2,5,4)
            addTriFace(faces,2,nodes[2],nodes[4],nodes[1]);
            addTriFace(faces,2,nodes[2],nodes[4],nodes[5]);
            /// face (0,1,4,3)
            addTriFace(faces,1,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][2]]]);
            addTriFace(faces,1,nodes[face_nodes[1][diag_nodes[1][0]]],nodes[face_nodes[1][diag_nodes[1][1]]],nodes[face_nodes[1][diag_nodes[1][3]]]);
        }

        /// tetrahedron0
        addTetraCell(cells,nodes[nodemap[0]],nodes[nodemap[1]],nodes[nodemap[2]],nodes[nodemap[3]]);
        /// tetrahedron1
        addTetraCell(cells,nodes[nodemap[4]],nodes[nodemap[5]],nodes[nodemap[6]],nodes[nodemap[0]]);
        /// tetrahedron2
        addTetraCell(cells,nodes[nodemap[4]],nodes[nodemap[5]],nodes[nodemap[7]],nodes[nodemap[0]]);

        return true;
    }

    /// in case if one deals with undividable case
    if( undividable )
    {
        /// add cell centroid
        GridNode* cgn;
        addCellCentroid(pgm,additional_points,additional_edges,cgn,nodes,face_nodes);
        addExtraNode( extra_nodes, cgn );

        for( size_t fid = 0; fid<num_quad_faces; ++fid )
        {
            assert( !diag_nodes[fid].empty() );

            /// tetrahedron0
            addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]],cgn);
            /// tetrahedron1
            addTetraCell(cells,nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]],cgn);

            /// triangle0
            addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][2]]]);
            /// triangle1
            addTriFace(faces,fid, nodes[face_nodes[fid][diag_nodes[fid][0]]],nodes[face_nodes[fid][diag_nodes[fid][1]]],nodes[face_nodes[fid][diag_nodes[fid][3]]]);
        }
        /// tetrahedron0 ( 0,1,2,6 )
        addTetraCell(cells,nodes[0],nodes[1],nodes[2],cgn);
        /// tetrahedron1 ( 3,4,5,6 )
        addTetraCell(cells,nodes[3],nodes[4],nodes[5],cgn);
        /// triangle0 ( 0,1,2 )
        addTriFace(faces,3, nodes[0],nodes[1],nodes[2]);
        /// triangle1 ( 3,4,5 )
        addTriFace(faces,4, nodes[3],nodes[4],nodes[5]);

        return true;
    }
    return false;
}











/**
    process pyramid ( expected node and face numbering )

                      4
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
          3 \   /         2 \   \
             \  /            \  \
              \ /             \ \
               \x______________\x
              0                 1

Faces:  bottom quad face ( 0321 )
        left tri face    ( 043  )
        front tri face   ( 014  )
        right tri face   ( 124  )
        back tri face    ( 234  )
 */

bool processPyramidElement( bool tetra_mesh,
                           PolygonGridManager* pgm,
                           std::map<std::set<size_t>,GridNode*>& additional_points,
                           std::set<std::set<size_t> >& additional_edges,
                           std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                           std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                           std::vector<GridNode*>& extra_nodes,
                           const std::vector<GridNode*>& nodes,
                           const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// diagonal nodes
    const size_t num_faces( face_nodes.size() );
    std::vector<size_t> diag_nodes;
    size_t valid_edge( 0 );

    valid_edge = isEdgeExist( additional_edges, diag_nodes,
                              nodes[face_nodes[0][0]],
                              nodes[face_nodes[0][1]],
                              nodes[face_nodes[0][2]],
                              nodes[face_nodes[0][3]] );
    if( valid_edge == 3 )
    {
        csmp_error.notice( csmp::FATAL_ERROR,
                           "processPyramidElement()",
                           "Two diagonals exist!!! Unpredictable case!!!" );
        return false;
    }
    else
    {
         if( valid_edge == 0 )
         {
             if( !tetra_mesh )
             {
                 /// quadrilateral
                 addPyramidCell(cells, nodes[0],nodes[3],nodes[2],nodes[1],nodes[4]);

                 /// add faces
                 for( size_t fid = 0; fid<num_faces; ++fid )
                 {
                     if( face_nodes[fid].size() == 4 )
                         addQuadFace(faces,fid, nodes[face_nodes[0][0]],nodes[face_nodes[0][1]],nodes[face_nodes[0][2]],nodes[face_nodes[0][3]] );
                     else if( face_nodes[fid].size() == 3 )
                         addTriFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
                 }
                 return true;
             }
             else
             {
                 /// face0: diagonal (1,3),(0,2)
                 valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes,
                                                                         nodes[face_nodes[0][0]],
                                                                         nodes[face_nodes[0][1]],
                                                                         nodes[face_nodes[0][2]],
                                                                         nodes[face_nodes[0][3]] );
                 if( valid_edge == 3 )
                 {
                     csmp_error.notice( csmp::FATAL_ERROR,
                                        "processPyramidElement()",
                                        "Two diagonals exist!!! Unpredictable case!!!" );
                     return false;
                 }
                 else if( valid_edge == 0 )
                 {
                     csmp_error.notice( csmp::FATAL_ERROR,
                                        "processPyramidElement()",
                                        "Both diagonals are invalid!!! Unpredictable case!!!" );
                     return false;
                 }
             }
         }
         /// tetra0
         addTetraCell(cells, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[2]]],nodes[4]);
         /// tetra1
         addTetraCell(cells, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[3]]],nodes[4]);

         /// add faces
         for( size_t fid = 0; fid<num_faces; ++fid )
         {
             if( face_nodes[fid].size() == 4 )
             {
                 /// triangle0
                 addTriFace(faces,fid, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[2]]] );
                 /// triangle1
                 addTriFace(faces,fid, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[3]]] );
             }
             else if( face_nodes[fid].size() == 3 )
                 addTriFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
         }
         return true;
    }
    return false;
}








/**

  tetrahedron:

                3
                x
                |\
                |\ \
                | \  \
                | \    \
                |  \     \
                |  \      \
              0 x___\_______x 2
                 \   \      /
                  \  \     /
                   \  \   /
                    \  \ /
                     \ \/
                      \x
                       1        return 1;

Faces:  tri face ( 021 )
        tri face ( 013 )
        tri face ( 123 )
        tri face ( 032 )
*/

bool processTetrahedronElement( bool tetra_mesh,
                               PolygonGridManager* pgm,
                               std::map<std::set<size_t>,GridNode*>& additional_points,
                               std::set<std::set<size_t> >& additional_edges,
                               std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                               std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                               std::vector<GridNode*>& extra_nodes,
                               const std::vector<GridNode*>& nodes,
                               const std::vector<std::vector<size_t> >& face_nodes )
{
    addTetraCell(cells,nodes[0],nodes[1],nodes[2],nodes[3]);
    /// add faces
    const size_t num_faces( face_nodes.size() );
    for( size_t fid = 0; fid<num_faces; ++fid )
    {
        if( face_nodes[fid].size() == 3 )
            addTriFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
    }
    return true;
}









/**
    process quadrilateral ( expected node and face numbering )

           3               2
           x_______________x
            \               \
             \               \
              \               \nodes[0],nodes[3],nodes[2],nodes[1]
               \x______________\x
              0                 1

Faces:  quad face ( 0321 )

 */

bool processQuadrilateralElement( bool tetra_mesh,
                                 PolygonGridManager* pgm,
                                 std::map<std::set<size_t>,GridNode*>& additional_points,
                                 std::set<std::set<size_t> >& additional_edges,
                                 std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                                 std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                                 std::vector<GridNode*>& extra_nodes,
                                 const std::vector<GridNode*>& nodes,
                                 const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    /// diagonal nodes
    const size_t num_faces( face_nodes.size() );
    std::vector<size_t> diag_nodes;
    size_t valid_edge( 0 );

    valid_edge = isEdgeExist( additional_edges, diag_nodes,
                              nodes[face_nodes[0][0]],
                              nodes[face_nodes[0][1]],
                              nodes[face_nodes[0][2]],
                              nodes[face_nodes[0][3]] );
    if( valid_edge == 3 )
    {
        csmp_error.notice( csmp::FATAL_ERROR,
                           "processQuadrilateralElement()",
                           "Two diagonals exist!!! Unpredictable case!!!" );
        return false;
    }
    else
    {
         if( valid_edge == 0 )
         {
             if( !tetra_mesh )
             {
                 /// quadrilateral
                 addQuadCell(cells, nodes[0],nodes[3],nodes[2],nodes[1]);

                 /// add faces
                 for( size_t fid = 0; fid<num_faces; ++fid )
                 {
                     if( face_nodes[fid].size() == 4 )
                         addQuadFace(faces,fid, nodes[face_nodes[0][0]],nodes[face_nodes[0][1]],nodes[face_nodes[0][2]],nodes[face_nodes[0][3]]);
                     else if( face_nodes[fid].size() == 2 )
                         addBarFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]] );
                 }
                 return true;
             }
             else
             {
                 /// face0: diagonal (1,3),(0,2)
                 valid_edge = addValidAndShortestOrBiggerSolidAngleEdge( additional_edges, diag_nodes,
                                                                         nodes[face_nodes[0][0]],
                                                                         nodes[face_nodes[0][1]],
                                                                         nodes[face_nodes[0][2]],
                                                                         nodes[face_nodes[0][3]] );
                 if( valid_edge == 3 )
                 {
                     csmp_error.notice( csmp::FATAL_ERROR,
                                        "processQuadrilateralElement()",
                                        "Two diagonals exist!!! Unpredictable case!!!" );
                     return false;
                 }
                 else if( valid_edge == 0 )
                 {
                     csmp_error.notice( csmp::FATAL_ERROR,
                                        "processQuadrilateralElement()",
                                        "Both diagonals are invalid!!! Unpredictable case!!!" );
                     return false;
                 }
             }
         }
         /// triangle0
         addTriCell(cells, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[2]]] );
         /// triangle1
         addTriCell(cells, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[3]]] );

         /// add faces
         for( size_t fid = 0; fid<num_faces; ++fid )
         {
             if( face_nodes[fid].size() == 4 )
             {
                 /// triangle0
                 addTriFace(faces,fid, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[2]]] );
                 /// triangle1
                 addTriFace(faces,fid, nodes[face_nodes[0][diag_nodes[0]]],nodes[face_nodes[0][diag_nodes[1]]],nodes[face_nodes[0][diag_nodes[3]]] );
             }
             else if( face_nodes[fid].size() == 2 )
                 addBarFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]] );
         }
         return true;
    }
    return false;
}








/**

  triangle:

      0  x____________x 1
          \          /
           \        /
            \      /
             \    /
              \  /
               \/
               x 2

Faces:  tri face ( 012 )

*/

bool processTriangleElement( bool tetra_mesh,
                            PolygonGridManager* pgm,
                            std::map<std::set<size_t>,GridNode*>& additional_points,
                            std::set<std::set<size_t> >& additional_edges,
                            std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                            std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                            std::vector<GridNode*>& extra_nodes,
                            const std::vector<GridNode*>& nodes,
                            const std::vector<std::vector<size_t> >& face_nodes )
{
    addTriCell(cells,nodes[0],nodes[1],nodes[2]);
    /// add faces
    const size_t num_faces( face_nodes.size() );
    for( size_t fid = 0; fid<num_faces; ++fid )
    {
        if( face_nodes[fid].size() == 3 )
            addTriFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]],nodes[face_nodes[fid][2]] );
        else if( face_nodes[fid].size() == 2 )
            addBarFace(faces,fid, nodes[face_nodes[fid][0]],nodes[face_nodes[fid][1]] );
    }
    return true;
}









/**

  bar:
      0  x____________x 1

Faces: line face ( 01 )

*/

bool processBarElement( bool tetra_mesh,
                       PolygonGridManager* pgm,
                       std::map<std::set<size_t>,GridNode*>& additional_points,
                       std::set<std::set<size_t> >& additional_edges,
                       std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                       std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                       std::vector<GridNode*>& extra_nodes,
                       const std::vector<GridNode*>& nodes,
                       const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    addBarCell(cells,nodes[0],nodes[1]);
    return true;
}









/**
  point:
            x 0

Faces: point face ( 0 )

*/

bool processPointElement( bool tetra_mesh,
                         PolygonGridManager* pgm,
                         std::map<std::set<size_t>,GridNode*>& additional_points,
                         std::set<std::set<size_t> >& additional_edges,
                         std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                         std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                         std::vector<GridNode*>& extra_nodes,
                         const std::vector<GridNode*>& nodes,
                         const std::vector<std::vector<size_t> >& face_nodes )
{
    csmp::ErrorHandler& csmp_error( csmp::ErrorHandler::Instance() );

    addPointCell(cells,nodes[0]);
    return true;
}

} // eclipse

} // end namespace csmp









