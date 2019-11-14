#ifndef MESHING_TOOLS_UOM_H
#define MESHING_TOOLS_UOM_H

#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"
#include "Point.h"
#include "PolygonGrid.h"
#include "GeometricCalculations.h"


/**
    Collection of functions that process correctly numbered finite-element cells from the supplied point vectors.
*/


// TODO: clean up function arguments and make these global functions member of classes
// TODO: typedef some of the nested containers

// these functions are used by the Eclipse interface
namespace csmp {

/// vector of pairs of FE-type - GridNode pointer vectors
 struct Cell {
    typedef std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >  Vector;
 };

/// gathering cells

bool checkCell( const std::vector<GridNode*>& cell );

/// alternative form of addHexaCell; @todo create struct to group swath of arguments

void addHexaCellTo( Cell::Vector&, GridNode* const, GridNode* const, GridNode* const, GridNode* const,
                    GridNode* const, GridNode* const, GridNode* const, GridNode* const );


void addHexaCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                  GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4, GridNode* n5, GridNode* n6, GridNode* n7 );


void addPrismCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                   GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4, GridNode* n5 );


void addPyramidCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                     GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3, GridNode* n4 );


void addTetraCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                   GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 );


void addQuadCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 );


void addTriCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1, GridNode* n2 );


void addBarCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0, GridNode* n1 );


void addPointCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells, GridNode* n0 );

/// gathering faces

void addQuadFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                  size_t fid, GridNode* n0, GridNode* n1, GridNode* n2, GridNode* n3 );


void addTriFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces,
                 size_t fid, GridNode* n0, GridNode* n1, GridNode* n2 );


void addBarFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0, GridNode* n1 );


void addPointFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > >& faces, size_t fid, GridNode* n0 );

/// gathering extra nodes

void addExtraNode( std::vector<GridNode*>& extra_nodes, GridNode* n );

/// additional edges


size_t isEdgeExist( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t chooseValidAndPreferablyFirstEdge( std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t chooseValidAndPreferablySecondEdge( std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t chooseValidAndShortestOrBiggerSolidAngleEdge( std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t addValidAndPreferablyFirstEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t addValidAndPreferablySecondEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


size_t addValidAndShortestOrBiggerSolidAngleEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );

/// rename and check whether 'additional_points' already contains a point set that is equal to the supplied points; if so the corresponging grid node is returned into point; also only 3D version is needed
 // TODO: replace *& by an const_iterator to the add-points map, or return pair<iterator,bool> to get desired output
bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points,
                   GridNode*& pt, GridNode* const pt0, GridNode* const pt1, GridNode* const pt2, GridNode* const pt3 );



bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points, GridNode*& pt, GridNode* const pt0, GridNode* const pt1, GridNode* const pt2 );


bool doesPointExist( const std::map<std::set<size_t>,GridNode*>& additional_points, GridNode*& pt, GridNode* const pt0, GridNode* const pt1 );


bool containRemeshedFaces( const std::map<std::set<size_t >,GridNode*>& additional_points, const std::vector<GridNode*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );


bool containRemeshedEdges( const std::map<std::set<size_t >,GridNode*>& additional_points, const std::vector<GridNode*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );


bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode& pt, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode& pt, GridNode* pt0, GridNode* pt1, GridNode* pt2 );


bool chooseFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode& pt, GridNode* pt0, GridNode* pt1 );


void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode* pt, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode* pt, GridNode* pt0, GridNode* pt1, GridNode* pt2 );


void addExistingFaceCentroid( std::map<std::set<size_t>,GridNode*>& additional_points, GridNode* pt, GridNode* pt0, GridNode* pt1 );


 // TODO: check whether all these methods want to be part of the PolygonGridManager ?
bool addFaceCentroid( PolygonGridManager* pgm, std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& pt, GridNode* pt0, GridNode* pt1, GridNode* pt2, GridNode* pt3 );


bool addFaceCentroid( PolygonGridManager* pgm, std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& pt, GridNode* pt0, GridNode* pt1, GridNode* pt2 );


bool addFaceCentroid( PolygonGridManager* pgm, std::map<std::set<size_t>,GridNode*>& additional_points,
                      GridNode*& pt, GridNode* pt0, GridNode* pt1 );

/// extra useful meshing techiques

bool addCellCentroid( PolygonGridManager* pgm, std::map<std::set<size_t >,GridNode*>& additional_points,
                      std::set<std::set<size_t > >& additional_edges,
                      GridNode*& cgn, const std::vector<GridNode*>& nodes,
                      const std::vector<std::vector<size_t> >& quad_face_nodes );


bool addCellFaceCentroids( PolygonGridManager* pgm, std::map<std::set<size_t >,GridNode*>& additional_points,
                           std::set<std::set<size_t > >& additional_edges,
                           const std::vector<GridNode*>& nodes,
                           const std::vector<std::vector<size_t> >& face_nodes );


bool addOverlappingCellFaceCentroids( PolygonGridManager* pgm, std::map<std::set<size_t >,GridNode*>& additional_points,
                                      std::set<std::set<size_t > >& additional_edges,
                                      GridNode*& cgn, const std::vector<GridNode*>& nodes,
                                      const std::vector<std::vector<size_t> >& face_nodes );

/// processing certain elements

// TODO: give these typedefs class scope, put it into CornerPointCell
/// faces[i][j] is a 6 x 5 matrix of faces defined by face-type and corresponding node pointers as entries (type, value pairs)
typedef std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > > CellFaces;
// ???
typedef std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > > Cells;


// TODO: reduce this crazy number of arguments, these functions should be member functions of the CornerPointCell or the PolygonCell (what is the difference between these?)
// in fact, these classes already take the pgm as a constructor argument !!!
// break up - far too many options make debugging impossible
bool processCellWithCentroids( bool add_edge_centroid, bool add_face_centroid, bool add_cell_centroid, bool is_volumetric_element, bool tetra_mesh,
                               PolygonGridManager* pgm,
                               std::map<std::set<size_t >,GridNode*>& additional_points,
                               std::set<std::set<size_t > >& additional_edges,
                               std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<GridNode*> > >& cells,
                               CellFaces& faces,
                               std::vector<GridNode* >& extra_nodes,
                               const std::vector<GridNode*>& nodes,
                               const std::vector<std::vector<size_t> >& face_nodes );

// TODO: tetra mesh is an option that will produce tetrahedra from the cell, but this should rather be a seperate method
bool processDegenerateHexahedronElement( bool tetra_mesh,
                                         PolygonGridManager* pgm,
                                         std::map<std::set<size_t>,GridNode*>& additional_points,
                                         std::set<std::set<size_t> >& additional_edges,
                                         Cells& cells,
                                         CellFaces& faces,
                                         std::vector<GridNode* >& extra_nodes,
                                         const std::vector<GridNode*>& nodes,
                                         const std::vector<std::vector<size_t> >& face_nodes );


bool processHexahedronElement( bool tetra_mesh,
                               PolygonGridManager* pgm,
                               std::map<std::set<size_t>,GridNode*>& additional_points,
                               std::set<std::set<size_t> >& additional_edges,
                               Cells& cells,
                               CellFaces& faces,
                               std::vector<GridNode* >& extra_nodes,
                               const std::vector<GridNode*>& nodes,
                               const std::vector<std::vector<size_t> >& face_nodes );


bool processPrismElement( bool tetra_mesh,
                          PolygonGridManager* pgm,
                          std::map<std::set<size_t>,GridNode*>& additional_points,
                          std::set<std::set<size_t> >& additional_edges,
                          Cells& cells,
                          CellFaces& faces,
                          std::vector<GridNode* >& extra_nodes,
                          const std::vector<GridNode*>& nodes,
                          const std::vector<std::vector<size_t> >& face_nodes );


bool processPyramidElement( bool tetra_mesh,
                            PolygonGridManager* pgm,
                            std::map<std::set<size_t>,GridNode*>& additional_points,
                            std::set<std::set<size_t> >& additional_edges,
                            Cells& cells,
                            CellFaces& faces,
                            std::vector<GridNode* >& extra_nodes,
                            const std::vector<GridNode*>& nodes,
                            const std::vector<std::vector<size_t> >& face_nodes );


bool processTetrahedronElement( bool tetra_mesh,
                                PolygonGridManager* pgm,
                                std::map<std::set<size_t>,GridNode*>& additional_points,
                                std::set<std::set<size_t> >& additional_edges,
                                Cells& cells,
                                CellFaces& faces,
                                std::vector<GridNode* >& extra_nodes,
                                const std::vector<GridNode*>& nodes,
                                const std::vector<std::vector<size_t> >& face_nodes );


bool processQuadrilateralElement( bool tetra_mesh,
                                  PolygonGridManager* pgm,
                                  std::map<std::set<size_t>,GridNode*>& additional_points,
                                  std::set<std::set<size_t> >& additional_edges,
                                  Cells& cells,
                                  CellFaces& faces,
                                  std::vector<GridNode* >& extra_nodes,
                                  const std::vector<GridNode*>& nodes,
                                  const std::vector<std::vector<size_t> >& face_nodes );


bool processTriangleElement( bool tetra_mesh,
                             PolygonGridManager* pgm,
                             std::map<std::set<size_t>,GridNode*>& additional_points,
                             std::set<std::set<size_t> >& additional_edges,
                             Cells& cells,
                             CellFaces& faces,
                             std::vector<GridNode* >& extra_nodes,
                             const std::vector<GridNode*>& nodes,
                             const std::vector<std::vector<size_t> >& face_nodes );


bool processBarElement( bool tetra_mesh,
                        PolygonGridManager* pgm,
                        std::map<std::set<size_t>,GridNode*>& additional_points,
                        std::set<std::set<size_t> >& additional_edges,
                        Cells& cells,
                        CellFaces& faces,
                        std::vector<GridNode* >& extra_nodes,
                        const std::vector<GridNode*>& nodes,
                        const std::vector<std::vector<size_t> >& face_nodes );


bool processPointElement( bool tetra_mesh,
                          PolygonGridManager* pgm,
                          std::map<std::set<size_t>,GridNode*>& additional_points,
                          std::set<std::set<size_t> >& additional_edges,
                          Cells& cells,
                          CellFaces& faces,
                          std::vector<GridNode* >& extra_nodes,
                          const std::vector<GridNode*>& nodes,
                          const std::vector<std::vector<size_t> >& face_nodes );

/**

@author R. Manasipov
@date 2015

Common gridding utilities relying on a suite of classes in other directories: GridNode, PoygonGridManager

@todo SKM: turn code into C++ which will also make it a lot more efficient
@todo SKM: refactor GeometricCalculations collection of functions to make design cleared

*/

}// end namespace csmp

#endif

