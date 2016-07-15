#ifndef MESHING_TOOLS_H
#define MESHING_TOOLS_H

#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"
#include "Point.h"
#include "PolygonGrid.h"
#include "GeometricCalculations.h"

// TODO: clean up function arguments and make these global functions member of classes

// these functions are used by the Eclipse interface
namespace csmp {

/// vector of pairs of FE-type - GridNode pointer vectors
template<size_t dim> struct Cell {
    typedef std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >  Vector;
 };

/// gathering cells
template<size_t dim>
bool checkCell( const std::vector<csmp::GridNode<dim>*>& cell );

/// alternative form of addHexaCell; @todo create struct to group swath of arguments
template<size_t dim>
void addHexaCellTo( typename Cell<dim>::Vector&, GridNode<dim>* const, GridNode<dim>* const, GridNode<dim>* const, GridNode<dim>* const,
                                                 GridNode<dim>* const, GridNode<dim>* const, GridNode<dim>* const, GridNode<dim>* const );

template<size_t dim>
void addHexaCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3, csmp::GridNode<dim>* n4, csmp::GridNode<dim>* n5, csmp::GridNode<dim>* n6, csmp::GridNode<dim>* n7 );

template<size_t dim>
void addPrismCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3, csmp::GridNode<dim>* n4, csmp::GridNode<dim>* n5 );

template<size_t dim>
void addPyramidCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3, csmp::GridNode<dim>* n4 );

template<size_t dim>
void addTetraCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3 );

template<size_t dim>
void addQuadCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3 );

template<size_t dim>
void addTriCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2 );

template<size_t dim>
void addBarCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1 );

template<size_t dim>
void addPointCell( std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, csmp::GridNode<dim>* n0 );

/// gathering faces
template<size_t dim>
void addQuadFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, size_t fid, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2, csmp::GridNode<dim>* n3 );

template<size_t dim>
void addTriFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, size_t fid, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1, csmp::GridNode<dim>* n2 );

template<size_t dim>
void addBarFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, size_t fid, csmp::GridNode<dim>* n0, csmp::GridNode<dim>* n1 );

template<size_t dim>
void addPointFace( std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, size_t fid, csmp::GridNode<dim>* n0 );

/// gathering extra nodes
template<size_t dim>
void addExtraNode( std::vector<csmp::GridNode<dim>*>& extra_nodes, csmp::GridNode<dim>* n );

/// additional edges

template<size_t dim>
size_t isEdgeExist( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t chooseValidAndPreferablyFirstEdge( std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t chooseValidAndPreferablySecondEdge( std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t chooseValidAndShortestOrBiggerSolidAngleEdge( std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t addValidAndPreferablyFirstEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t addValidAndPreferablySecondEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
size_t addValidAndShortestOrBiggerSolidAngleEdge( std::set<std::set<size_t> >& additional_edges, std::vector<size_t>& ids, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

/// rename and check whether 'additional_points' already contains a point set that is equal to the supplied points; if so the corresponging grid node is returned into point; also only 3D version is needed
template<size_t dim> // TODO: replace *& by an const_iterator to the add-points map, or return pair<iterator,bool> to get desired output
bool isPointExist( const std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points,
                   csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* const pt0, csmp::GridNode<dim>* const pt1, csmp::GridNode<dim>* const pt2, csmp::GridNode<dim>* const pt3 );


template<size_t dim>
bool isPointExist( const std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* const pt0, csmp::GridNode<dim>* const pt1, csmp::GridNode<dim>* const pt2 );

template<size_t dim>
bool isPointExist( const std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* const pt0, csmp::GridNode<dim>* const pt1 );

template<size_t dim>
bool isContainRemeshedFaces( std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool isContainRemeshedEdges( std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool chooseFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
bool chooseFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2 );

template<size_t dim>
bool chooseFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1 );

template<size_t dim>
void addExistingFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>* pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
void addExistingFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>* pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2 );

template<size_t dim>
void addExistingFaceCentroid( std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>* pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1 );

template<size_t dim>
bool addFaceCentroid( PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2, csmp::GridNode<dim>* pt3 );

template<size_t dim>
bool addFaceCentroid( PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1, csmp::GridNode<dim>* pt2 );

template<size_t dim>
bool addFaceCentroid( PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, csmp::GridNode<dim>*& pt, csmp::GridNode<dim>* pt0, csmp::GridNode<dim>* pt1 );

/// extra useful meshing techiques
template<size_t dim>
bool addCellCentroid( PolygonGridManager<dim>* pgm, std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t > >& additional_edges, csmp::GridNode<dim>*& cgn, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& quad_face_nodes );

template<size_t dim>
bool addCellFaceCentroids( PolygonGridManager<dim>* pgm, std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t > >& additional_edges, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool addOverlappingCellFaceCentroids( PolygonGridManager<dim>* pgm, std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t > >& additional_edges, csmp::GridNode<dim>*& cgn, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

/// processing certain elements
template<size_t dim>
bool processCellWithCentroids( bool add_edge_centroid, bool add_face_centroid, bool add_cell_centroid, bool is_volumetric_element, bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t >,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t > >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processDegenerateHexahedronElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processHexahedronElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processPrismElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processPyramidElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processTetrahedronElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processQuadrilateralElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells, std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces, std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes, const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processTriangleElement( bool tetra_mesh, PolygonGridManager<dim>* pgm,
                             std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points, std::set<std::set<size_t> >& additional_edges,
                             std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells,
                             std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces,
                             std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes,
                             const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processBarElement( bool tetra_mesh, PolygonGridManager<dim>* pgm, std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points,
                        std::set<std::set<size_t> >& additional_edges, std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells,
                        std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces,
                        std::vector<csmp::GridNode<dim>* >& extra_nodes, const std::vector<csmp::GridNode<dim>*>& nodes,
                        const std::vector<std::vector<size_t> >& face_nodes );

template<size_t dim>
bool processPointElement( bool tetra_mesh, PolygonGridManager<dim>* pgm,
                          std::map<std::set<size_t>,csmp::GridNode<dim>*>& additional_points,
                          std::set<std::set<size_t> >& additional_edges,
                          std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > >& cells,
                          std::vector<std::vector<std::pair<csmp::CSMP_FEM_TYPE,std::vector<csmp::GridNode<dim>*> > > >& faces,
                          std::vector<csmp::GridNode<dim>* >& extra_nodes,
                          const std::vector<csmp::GridNode<dim>*>& nodes,
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

