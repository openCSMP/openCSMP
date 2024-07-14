#include "vsetMakers.h"
#include "Box.h"
#include "CSMP_definitions.h"
#include "CSMP_global_enumerations.h"
#include "TensorVariable.h"

#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearQuadrilateral.h"




///this variable controls the perturbation of the node when it is skewed, the smaller it is, the smaller the perturbation will be
#define PERTURBATION 0.0001

using namespace std;

namespace csmp {

/**
    creates Mesh with 4 linear quadrilateral elements
*/
VSet<2U> create_Quadrilateral_VSet()
{
    const size_t iNodes(9);
    const size_t iNrOfElements(4);
    
  	IsoparametricLinearQuadrilateral iso_quadrilateral;
  	VSet< 2>  vset( iso_quadrilateral.Nodes(),
                    iso_quadrilateral.Neighbors(),
                    iso_quadrilateral.ElementType(), iNodes, iNrOfElements );
  	
  	//define nodes
  	deque<double> px(iNodes);
  	deque<double> py(iNodes);
  	deque<double> pz(iNodes);
  
    // xy plane, left-right numbering from bottom y to top y
  	px[0]=0.;py[0]=0.;pz[0]=0.;
  	px[1]=1.;py[1]=0.;pz[1]=0.;
  	px[2]=2.;py[2]=0.;pz[2]=0.;
  	px[3]=0.;py[3]=1.;pz[3]=0.;
  	px[4]=1.;py[4]=1.;pz[4]=0.;
  	px[5]=2.;py[5]=1.;pz[5]=0.;
  	px[6]=0.;py[6]=2.;pz[6]=0.;
  	px[7]=1.;py[7]=2.;pz[7]=0.;
  	px[8]=2.;py[8]=2.;pz[8]=0.;
  	  
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    
    //define elements
    deque<vector<size_t> > deqElements(iNrOfElements);
    vector<size_t> vecNodes(4);
    //element 0
    vecNodes[0]=0;
    vecNodes[1]=1;
    vecNodes[2]=4;
    vecNodes[3]=3;
    deqElements[0]=vecNodes;
    //element 1
    vecNodes[0]=1;
    vecNodes[1]=2;
    vecNodes[2]=5;
    vecNodes[3]=4;
    deqElements[1]=vecNodes;
    //element 2
    vecNodes[0]=4;
    vecNodes[1]=5;
    vecNodes[2]=8;
    vecNodes[3]=7;
    deqElements[2]=vecNodes;
    //element 3
    vecNodes[0]=3;
    vecNodes[1]=4;
    vecNodes[2]=7;
    vecNodes[3]=6;
    deqElements[3]=vecNodes;
    
    vset.AddPlist( deqElements.begin(),deqElements.end());

    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    vector<int64_t> vecNeighbors(4);
    //element 0
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=2;
    vecNeighbors[2]=4;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[0]=vecNeighbors;
    //element 1
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=3;
    vecNeighbors[3]=1;
    deqElementNeighbors[1]=vecNeighbors;
    //element 2
    vecNeighbors[0]=2;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=4;
    deqElementNeighbors[2]=vecNeighbors;
    //element 3
    vecNeighbors[0]=1;
    vecNeighbors[1]=3;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[3]=vecNeighbors;
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//define boundaries
  	vset.BFlag( 0, CNR1 );
    vset.BFlag( 1, BOTTOM_OUTSIDE );
    vset.BFlag( 2, CNR2 );
    vset.BFlag( 3, LEFT_OUTSIDE );
    vset.BFlag( 5, RIGHT_OUTSIDE );
    vset.BFlag( 6, CNR4 );
	  vset.BFlag( 7, TOP_OUTSIDE );
    vset.BFlag( 8, CNR3 );
    
    vset.EstablishZeroBasedNumbering();
    vset.Out();
    
    return vset;
    
} // end create_QuadraticQUadrilateral_VSet








void create_1Square_VSet(VSet<2U>& vset, double length_of_sides, bool bSkewed )
{
  	IsoparametricLinearQuadrilateral iso_quad;
  	
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int8_t> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_QUADRILATERAL;
 
    const size_t   n_nodes(4), n_elmts{1};
    const uint32_t nodes_per_elmt{4}, nbors_per_elmt{4};
  	vset.Resize( nodes_per_elmt, nbors_per_elmt, ISOPARAMETRIC_LINEAR_QUADRILATERAL, n_nodes, n_elmts );
 
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(n_nodes);
  	deque<double> py(n_nodes);
  	deque<double> pz(n_nodes);
  	
  	px[0]=0;                py[0]=0;                  pz[0]=0;
  	px[1]=length_of_sides;  py[1]=0;                  pz[1]=0;
  	px[2]=length_of_sides;  py[2]=length_of_sides;    pz[2]=0;
  	px[3]=0;                py[3]=length_of_sides;    pz[3]=0;
  	
  	if( bSkewed )
  	 for ( unsigned int i = 0; i < 8; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
    //define quadrilateral elements (elements 0->26), assign nodes per element
    deque< vector<size_t> > deqElements(1);
    deqElements[0].resize(n_nodes);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(4);
  	deqElementNeighbors[0][0]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][1]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][2]= TOP_OUTSIDE;
  	deqElementNeighbors[0][3]= LEFT_OUTSIDE;
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 0, CNR1);
  	vset.BFlag( 1, CNR2);
  	vset.BFlag( 2, CNR3);
  	vset.BFlag( 3, CNR4);
  	
    vset.EstablishZeroBasedNumbering();
    //vset.Out();
    
} // end create_1Square_VSet





/**
    Model TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object at the box boundary.
    Model is rectangle shaped
*/
ModelTopology  create_SimplePolyElement2DModel( VSet<2U>& vset )
 {
    //--------------------------ELEMENT TYPES
  	//add element types
    const CSMP_FEM_TYPE T(ISOPARAMETRIC_LINEAR_TRIANGLE), Q(ISOPARAMETRIC_LINEAR_QUADRILATERAL), P(ISOPARAMETRIC_LINEAR_BAR);
    deque<int8_t> vecElementTypes = { T,T,Q,P,       // elements
                                      P,P,P,P,P,P }; // faces
    const int n_cells{10};
  	assert( vecElementTypes.size() == n_cells );
    const size_t     n_nodes(6); // number of nodes
  	deque<uint32_t>  npes(n_cells,2);  // default: number of nodes per element
    deque<uint32_t>  epes(n_cells,2);  // default: nbors per element
    for ( size_t i{0U}; i<vecElementTypes.size(); ++i ) {
         if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_TRIANGLE ) {
              npes[i] = 3;
              epes[i] = 3;
           }
         else if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) {
              npes[i] = 4;
              epes[i] = 4;
           }
      }
    const int n_faces{6}, n_interfaces{0};
  	vset.Resize( vecElementTypes, npes, epes, n_nodes, n_faces, n_interfaces );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES (6)
  	//define node coordinates
  	deque<double> px = { 0.,1.,2.,0.,1.,2. };
    assert( px.size() == n_nodes );
  	deque<double> py = { 1.,1.,1., 0.,0.,0. };
    assert( py.size() == n_nodes );
  	deque<double> pz(n_nodes,0.);
  	  	  	
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------NODE BOUNDARY FLAGS
    BOX_BOUNDARY B{BOTTOM}, R{RIGHT}, U{TOP}, L{LEFT};
    vector<int8_t> bflags = { CNR1, B, CNR2, CNR3, U, CNR4 };
    assert( bflags.size() == n_nodes );
    vset.AddBFlags( bflags.begin(), bflags.end() );


  	//--------------------------ELEMENTS
  	// define nodes per element
    deque< vector<size_t> > plist = { {0,3,1}, {1,3,4}, {1,4,5,2}, {1,3}, // 4 elements
                                       {3,4}, {4,5}, {5,2}, {2,1}, {1,0}, {0,3} }; // 17 faces
    assert( plist.size() == n_cells );
    vset.AddPlist( plist.begin(), plist.end() );

     //---------------------------------NEIGHBORS
    //define neighbors per element
    deque<vector<int64_t> > pfverts = { {1,U,L}, {B,2,0}, {1,B,R,U}, {CNR1,U}, // element neighbors
                                        // face neighbors: { face-nbors, connected high-dim elmts, local face ids of high dim elmts }
                                        {5,9,1,B,0,B}, {6,4,2,B,1,B}, {7,5,2,R,2,R}, {8,6,2,U,3,U}, {9,7,0,U,1,U}, {4,8,0,L,2,L} };
    assert( pfverts.size() == n_cells );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );

    // creating a matching model topology
    ModelTopology mesh_topology( "MODEL_TINY", true );
    
    // surface elements fall in 2 domains "lower" and "upper"
    mesh_topology.AddDomain( "TRIA", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE"},
                              vector<size_t>{ 0 } );
    mesh_topology.AddDomain( "MIXED", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE", "ISOPARAMETRIC_LINEAR_QUADRILATERAL"},
                              vector<size_t>{ 1,2 } );
    mesh_topology.AddDomain( "LINE", set<string>{"ISOPARAMETRIC_LINEAR_BAR"},
                              vector<size_t>{ 3 } );
    // boundaries
    mesh_topology.AddDomain( "BOTTOM", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{4,5} );
    mesh_topology.AddDomain( "RIGHT",  set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{6} );
    mesh_topology.AddDomain( "TOP",    set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{7,8} );
    mesh_topology.AddDomain( "LEFT",   set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{9} );

    assert( mesh_topology.Cells() == vset.Elements() + vset.Faces() + vset.Interfaces() );
    assert( mesh_topology.Cells() == vset.Cells() );

    // adding corresponding materials to VSet
    const size_t n_elements{4};
    vector<int32_t> pmtrl = { 1,2,2,3 };
    assert( pmtrl.size() == n_elements );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    // adding node and element numbers for comparisons
    PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );
    // permeability
    PropertyData perm( ELEMENT, SCALAR, 2U );
    perm.Reserve( vset.Elements() );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("TRIA"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-15 ) );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("MIXED"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-14 ) );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("LINE"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-12 ) );
    vset.AddData( "permeability", perm );

    vset.Out();
    
    return mesh_topology;
    
 } // end create_SimplePolyElement2DModel
 
 




void create_TrianglePatch_VSet( VSet<2U>& vset )
{    
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int8_t> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_TRIANGLE;
  	
    const size_t   n_nodes(9), n_elmts{10};
    const uint32_t nodes_per_elmt{3}, nbors_per_elmt{3};
  	vset.Resize( nodes_per_elmt, nbors_per_elmt, ISOPARAMETRIC_LINEAR_TRIANGLE, n_nodes, n_elmts );
    
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(n_nodes);
  	deque<double> py(n_nodes);
  	deque<double> pz(n_nodes);
  	
  	px[0]=0.;                 py[0]=100.;               pz[0]=0.;
  	px[1]=52.73842909594504;  py[1]=100.;               pz[1]=0.;
  	px[2]=100.;               py[2]=100.;               pz[2]=0.;
  	px[3]=53.51413499621665;  py[3]=76.29281878271485;  pz[3]=0.;
  	px[4]=100.;               py[4]=48.94918579813998;  pz[4]=0.;
  	px[5]=26.17050201164178;  py[5]=23.35089108917626;  pz[5]=0.;
  	px[6]=79.5002826553162;   py[6]=28.97475886614556;  pz[6]=0.;
  	px[7]=0.;                 py[7]=0.;                 pz[7]=0.;
  	px[8]=100.;               py[8]=0.;                 pz[8]=0.;
  	  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
  	//--------------------------ELEMENTS
  	// define nodes per element
  	vector<size_t> node_dummy(3);
    deque< vector<size_t> > deqElements(10,node_dummy);
    deqElements[0][0]= 0;
    deqElements[0][1]= 3;
    deqElements[0][2]= 1;
    deqElements[1][0]= 1;
    deqElements[1][1]= 3;
    deqElements[1][2]= 2;
    deqElements[2][0]= 0;
    deqElements[2][1]= 5;
    deqElements[2][2]= 3;
    deqElements[3][0]= 3;
    deqElements[3][1]= 5;
    deqElements[3][2]= 6;
    deqElements[4][0]= 3;
    deqElements[4][1]= 6;
    deqElements[4][2]= 4;
    deqElements[5][0]= 2;
    deqElements[5][1]= 3;
    deqElements[5][2]= 4;
    deqElements[6][0]= 0;
    deqElements[6][1]= 7;
    deqElements[6][2]= 5;
    deqElements[7][0]= 5;
    deqElements[7][1]= 7;
    deqElements[7][2]= 8;
    deqElements[8][0]= 5;
    deqElements[8][1]= 8;
    deqElements[8][2]= 6;
    deqElements[9][0]= 4;
    deqElements[9][1]= 6;
    deqElements[9][2]= 8;

    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors per element
  	vector<int64_t> nbor_dummy(3);
    deque<vector<int64_t> > deqElementNeighbors(10,nbor_dummy);
    deqElementNeighbors[0][0]= 1;
    deqElementNeighbors[0][1]= TOP_OUTSIDE;
    deqElementNeighbors[0][2]= 2;
    deqElementNeighbors[1][0]= 5;
    deqElementNeighbors[1][1]= TOP_OUTSIDE;
    deqElementNeighbors[1][2]= 0;
    deqElementNeighbors[2][0]= 3;
    deqElementNeighbors[2][1]= 0;
    deqElementNeighbors[2][2]= 6;
    deqElementNeighbors[3][0]= 8;
    deqElementNeighbors[3][1]= 4;
    deqElementNeighbors[3][2]= 2;
    deqElementNeighbors[4][0]= 9;
    deqElementNeighbors[4][1]= 5;
    deqElementNeighbors[4][2]= 3;
    deqElementNeighbors[5][0]= 4;
    deqElementNeighbors[5][1]= RIGHT_OUTSIDE;
    deqElementNeighbors[5][2]= 1;
    deqElementNeighbors[6][0]= 7;
    deqElementNeighbors[6][1]= 2;
    deqElementNeighbors[6][2]= LEFT_OUTSIDE;
    deqElementNeighbors[7][0]= BOTTOM_OUTSIDE;
    deqElementNeighbors[7][1]= 8;
    deqElementNeighbors[7][2]= 6;
    deqElementNeighbors[8][0]= 9;
    deqElementNeighbors[8][1]= 3;
    deqElementNeighbors[8][2]= 7;
    deqElementNeighbors[9][0]= 8;
    deqElementNeighbors[9][1]= RIGHT_OUTSIDE;
    deqElementNeighbors[9][2]= 4;
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 7, CNR1);
  	vset.BFlag( 8, CNR2);
  	vset.BFlag( 2, CNR3);
  	vset.BFlag( 0, CNR4);
  	vset.BFlag( 1, TOP_OUTSIDE );
  	vset.BFlag( 4, RIGHT_OUTSIDE);
   
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    vset.Out();
    
} // end create_TrianglePatch_VSet






/**
        Test mixed element (triangle, quadrilateral and line elements) dataset for testing
        VData etc. functionality for the clean-up of 'pfverts', line-element connectivity etc.
                
        Also returns corresponding model topology and property data in the form of "element number" and "node number" for testing.
        
        @attention VSet does not contain any Face objects. It therefore needs to be used in conjunction with Element to Face conversion
        to create Boundary objects.
        
        @author SKM
        @date 1/10/2021
*/
ModelTopology create_MeshPatchWithLineElements_VSet( VSet<2U>& vset )
{
    //--------------------------ELEMENT TYPES
  	//add element types
    const CSMP_FEM_TYPE T(ISOPARAMETRIC_LINEAR_TRIANGLE), Q(ISOPARAMETRIC_LINEAR_QUADRILATERAL), B(ISOPARAMETRIC_LINEAR_BAR);
    deque<int8_t> vecElementTypes = { T,T,T,T,Q,T,T,T,T,T,Q,Q,Q,T,T,T,T,Q,T,T,T,T,T,T,T,
                                      B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B,B };
  	
    const size_t     nodes(22);   // number of nodes
  	deque<uint32_t>  npes(45,3);  // default: number of nodes per triangle
    deque<uint32_t>  epes(45,3);  // default: neighbor elements per triangle
    for ( size_t i{0U}; i<vecElementTypes.size(); ++i ) {
         if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_QUADRILATERAL ) {
              npes[i] = 4;
              epes[i] = 4;
           }
         // default: else if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_TRIANGLE ) {
         //       npes[i] = 3;
         //       epes[i] = 3;
         //   }
         else if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_BAR ) {
              npes[i] = 2;
              epes[i] = 2;
           }
      }
    
  	vset.Resize( vecElementTypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define node coordinates
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes,0.);
  	deque<int8_t> bflags(nodes,NOT), gflags(nodes,MESH_VERTEX);
  	
  	px[0]=0.;    py[0]=5.;    bflags[0] = CNR4;       gflags[0] = EXTERIOR_POINT;
  	px[1]=2.;    py[1]=5.;    bflags[1] = TOP;        gflags[1] = EXTERIOR_LINE;
  	px[2]=5.;    py[2]=5.;    bflags[2] = TOP;        gflags[2] = EXTERIOR_LINE;
  	px[3]=6;     py[3]=5.;    bflags[3] = CNR3;       gflags[3] = EXTERIOR_POINT;
  	px[4]=0.;    py[4]=4.;    bflags[4] = LEFT;       gflags[4] = EXTERIOR_LINE;
  	px[5]=2.;    py[5]=4.;    bflags[5] = INTERNAL;   gflags[5] = PERIMETER_POINT;
  	px[6]=4.;    py[6]=3.5;   bflags[6] = INTERNAL;   gflags[6] = INTERIOR_LINE;
  	px[7]=6.;    py[7]=3.;    bflags[7] = RIGHT;      gflags[7] = EXTERIOR_LINE;
  	px[8]=1.8;   py[8]=3.2;   bflags[8] = INTERNAL;   gflags[8] = PERIMETER_POINT;
  	px[9]=2.8;   py[9]=3.;    bflags[9] = INTERNAL;   gflags[9] = INTERSECTION_POINT;
  	px[10]=3.5;  py[10]=2.5;  bflags[10] = INTERNAL;  gflags[10] = INTERIOR_LINE;
  	px[11]=4.7;  py[11]=2.5;  bflags[11] = INTERNAL;  gflags[11] = PERIMETER_POINT;
  	px[12]=6.;   py[12]=2.;   bflags[12] = RIGHT;     gflags[12] = EXTERIOR_LINE;
  	px[13]=0.;   py[13]=2.5;  bflags[13] = LEFT;      gflags[13] = EXTERIOR_LINE;
  	px[14]=2.;   py[14]=1.5;  bflags[14] = INTERNAL;  gflags[14] = PERIMETER_POINT;
  	px[15]=4.2;  py[15]=1.;   bflags[15] = INTERNAL;  gflags[15] = PERIMETER_POINT;
  	px[16]=6.;   py[16]=0.;   bflags[16] = CNR2;      gflags[16] = EXTERIOR_POINT;
  	px[17]=0.;   py[17]=0.;   bflags[17] = CNR1;      gflags[17] = EXTERIOR_POINT;
  	px[18]=2.5;  py[18]=0.;   bflags[18] = BOTTOM;    gflags[18] = EXTERIOR_POINT;
  	px[19]=4.3;  py[19]=0.;   bflags[19] = BOTTOM;    gflags[19] = EXTERIOR_LINE;
  	px[20]=4.8;  py[20]=3.;   bflags[20] = INTERNAL;  gflags[20] = PERIMETER_POINT;
  	px[21]=5.;   py[21]=4.5;  bflags[21] = INTERNAL;  gflags[21] = PERIMETER_POINT;
  	  	  	
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------NODE BOUNDARY FLAGS
    size_t n(0U);
    for ( auto it : bflags )
  	  vset.BFlag( n++, it );

  	//--------------------------NODE GEOMETRY FLAGS
    n = 0U;
    for ( auto it : gflags )
  	  vset.BREP_Flag( n++, it );

  	//--------------------------ELEMENTS
  	// define nodes per element
    deque< vector<size_t> > deqElements(45);
    deqElements[0]  = { 1, 5, 6 };
    deqElements[1]  = { 0, 4, 5 };
    deqElements[2]  = { 0, 5, 1 };
    deqElements[3]  = { 5, 9, 6 };
    deqElements[4]  = { 1, 6, 21, 2 };
    deqElements[5]  = { 2, 21, 3 };
    deqElements[6]  = { 3, 21, 7 };
    deqElements[7]  = { 4, 8, 5 };
    deqElements[8]  = { 5, 8, 9 };
    deqElements[9]  = { 6, 9, 10 };
    deqElements[10] = { 6, 10, 11, 20 };
    deqElements[11] = { 6, 20, 7, 21 };
    deqElements[12] = { 7, 20, 11, 12 };
    deqElements[13] = { 4, 13, 8 };
    deqElements[14] = { 8, 13, 14 };
    deqElements[15] = { 8, 14, 10 };
    deqElements[16] = { 8, 10, 9 };
    deqElements[17] = { 10, 14, 18, 15 };
    deqElements[18] = { 10, 15, 11 };
    deqElements[19] = { 11, 15, 12 };
    deqElements[20] = { 12, 15, 16 };
    deqElements[21] = { 13, 17, 14 };
    deqElements[22] = { 14, 17, 18 };
    deqElements[23] = { 15, 18, 19 };
    deqElements[24] = { 15, 19, 16 };
    deqElements[25] = { 5, 9 };
    deqElements[26] = { 8, 9 };
    deqElements[27] = { 9, 6 };
    deqElements[28] = { 6, 21 };
    deqElements[29] = { 9, 10 };
    deqElements[30] = { 11, 20 };
    deqElements[31] = { 10, 15 };
    deqElements[32] = { 14, 18 };
    deqElements[33] = { 17, 18 };
    deqElements[34] = { 18, 19 };
    deqElements[35] = { 19, 16 };
    deqElements[36] = { 16, 12 };
    deqElements[37] = { 12, 7 };
    deqElements[38] = { 7, 3 };
    deqElements[39] = { 3, 2 };
    deqElements[40] = { 2, 1 };
    deqElements[41] = { 1, 0 };
    deqElements[42] = { 0, 4 };
    deqElements[43] = { 4, 13 };
    deqElements[44] = { 13, 17 };

    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors per element
    deque<vector<int64_t> > deqElementNeighbors(45);
    deqElementNeighbors[0]  = { 3, 4, 2 };
    deqElementNeighbors[1]  = { 7, 2, LEFT };
    deqElementNeighbors[2]  = { 0, TOP, 1 };
    deqElementNeighbors[3]  = { 9, 0, 8 };
    deqElementNeighbors[4]  = { 0, 11, 5, TOP };
    deqElementNeighbors[5]  = { 6, TOP, 4 };
    deqElementNeighbors[6]  = { 11, RIGHT, 5 };
    deqElementNeighbors[7]  = { 8, 1, 13 };
    deqElementNeighbors[8]  = { 16, 3, 7 };
    deqElementNeighbors[9]  = { 16, 10, 3 };
    deqElementNeighbors[10] = { 9, 18, 12, 11 };
    deqElementNeighbors[11] = { 10, 12, 6, 4 };
    deqElementNeighbors[12] = { 11, 10, 19, RIGHT };
    deqElementNeighbors[13] = { 14, 7, LEFT };
    deqElementNeighbors[14] = { 21, 15, 13 };
    deqElementNeighbors[15] = { 17, 16, 14 };
    deqElementNeighbors[16] = { 9, 8, 15 };
    deqElementNeighbors[17] = { 15, 22, 23, 18 };
    deqElementNeighbors[18] = { 19, 10, 17 };
    deqElementNeighbors[19] = { 20, 12, 18 };
    deqElementNeighbors[20] = { 24, RIGHT, 19 };
    deqElementNeighbors[21] = { 22, 14, LEFT };
    deqElementNeighbors[22] = { BOTTOM, 17, 21 };
    deqElementNeighbors[23] = { BOTTOM, 24, 17 };
    deqElementNeighbors[24] = { BOTTOM, 20, 23 };
    // for line elements, the faces/neighbors are located with the corner nodes,
    // but opposite to the nodes with the same number
    deqElementNeighbors[25] = { 29, INTERNAL };
    deqElementNeighbors[26] = { 27, INTERNAL };
    deqElementNeighbors[27] = { 28, 26 };
    deqElementNeighbors[28] = { INTERNAL, 27 };
    deqElementNeighbors[29] = { 31, 25 };
    deqElementNeighbors[30] = { INTERNAL, INTERNAL };
    deqElementNeighbors[31] = { INTERNAL, 29 };
    deqElementNeighbors[32] = { BOTTOM, INTERNAL };
    deqElementNeighbors[33] = { 34, 44 }; // boundary edges
    deqElementNeighbors[34] = { 35, 33 };
    deqElementNeighbors[35] = { 36, 34 };
    deqElementNeighbors[36] = { 37, 35 };
    deqElementNeighbors[37] = { 38, 36 };
    deqElementNeighbors[38] = { 39, 37 };
    deqElementNeighbors[39] = { 40, 38 };
    deqElementNeighbors[40] = { 41, 39 };
    deqElementNeighbors[41] = { 42, 40 };
    deqElementNeighbors[42] = { 43, 41 };
    deqElementNeighbors[43] = { 44, 42 };
    deqElementNeighbors[44] = { 33, 43 };
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());

    // creating a matching model topology
    ModelTopology mesh_topology( "create_MeshPatchWithLineElements_VSet", true );
    // all surface elements are "MATRIX"
    mesh_topology.AddDomain( "MATRIX", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE", "ISOPARAMETRIC_LINEAR_QUADRILATERAL"},
                              vector<size_t>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24} );

    // fracture line-element regions
    mesh_topology.AddDomain( "FRAC1", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{25,29,31} );
    mesh_topology.AddDomain( "FRAC2", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{26,27,28} );
    mesh_topology.AddDomain( "FRAC3", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{30} );
    mesh_topology.AddDomain( "FRAC4", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{32} );
    // boundaries
    mesh_topology.AddDomain( "BOTTOM", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{33,34,35} );
    mesh_topology.AddDomain( "RIGHT",  set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{36,37,38} );
    mesh_topology.AddDomain( "TOP",    set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{39,40,41} );
    mesh_topology.AddDomain( "LEFT",   set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{42,43,44} );

    assert( mesh_topology.Cells() == vset.Elements() );

    // adding corresponding materials to VSet
    vector<int32_t> pmtrl(45,1); // matrix
    fill( next(pmtrl.begin(),25), next(pmtrl.begin(),31), 2 ); // fine because wrong values will be overwritten next
    fill( next(pmtrl.begin(),26), next(pmtrl.begin(),28), 3 );
    fill( next(pmtrl.begin(),33), next(pmtrl.begin(),35), 4 );
    fill( next(pmtrl.begin(),36), next(pmtrl.begin(),38), 5 );
    fill( next(pmtrl.begin(),39), next(pmtrl.begin(),41), 6 );
    fill( next(pmtrl.begin(),42), pmtrl.end(), 7 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    // adding node and element numbers for comparisons
    // element number
    PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );
    // node number
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, i ) );
    vset.AddData( "node number", node_nums );
    // element variable (tetra=3. quadrilateral-4.)
    PropertyData elmt_vars( ELEMENT, SCALAR, 2U );
    elmt_vars.Reserve( vset.Elements() );
    // tets
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    // quad
    pushBack( elmt_vars, makeScalar( ANY, 4. ) );
    // tets
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    // quad
    pushBack( elmt_vars, makeScalar( ANY, 4. ) );
    pushBack( elmt_vars, makeScalar( ANY, 4. ) );
    pushBack( elmt_vars, makeScalar( ANY, 4. ) );
    // tets
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    // quad
    pushBack( elmt_vars, makeScalar( ANY, 4. ) );
    // tets
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    pushBack( elmt_vars, makeScalar( ANY, 3. ) );
    // bars
    for ( uint32_t i{0U}; i<20U; ++i )
      pushBack( elmt_vars, makeScalar( ANY, 2. ) );
      
    vset.AddData( "element variable", elmt_vars );
      
    // vset.Out();
    
    return mesh_topology;
    
} // end create_MeshPatchWithLineElements_VSet







/// model  SPLIT22_BASIC  with box boundaries (Faces) and one through-going and one internal crossing split boundary
ModelTopology  create_BoundarySplitBoundaryPatch( VSet<2U>& vset )
 {
    //--------------------------ELEMENT TYPES
  	//add element types
    const CSMP_FEM_TYPE T(ISOPARAMETRIC_LINEAR_TRIANGLE), Q(ISOPARAMETRIC_LINEAR_QUADRILATERAL), P(ISOPARAMETRIC_LINEAR_BAR);
    deque<int8_t> vecElementTypes = { Q,Q,Q,Q,T,Q,Q,Q,Q,T,Q,Q,Q,Q,T,Q,Q,Q,Q, // elements
                                      P,P,P, P,P,P,P,P, P,P,P, P,P,P,P,P,P,  // faces
                                      P,P,P,P, P,P,P };                      // interfaces
    const int n_cells{43};
  	assert( vecElementTypes.size() == n_cells );
    const size_t     n_nodes(35); // number of nodes
  	deque<uint32_t>  npes(n_cells,4);  // default: number of nodes per element
    deque<uint32_t>  epes(n_cells,4);  // default: nbors per element
    for ( size_t i{0U}; i<vecElementTypes.size(); ++i ) {
         if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_TRIANGLE ) {
              npes[i] = 3;
              epes[i] = 3;
           }
         else if ( vecElementTypes[i] == ISOPARAMETRIC_LINEAR_BAR ) {
              npes[i] = 2;
              epes[i] = 2;
           }
      }
    const int n_faces{17}, n_interfaces{7};
  	vset.Resize( vecElementTypes, npes, epes, n_nodes, n_faces, n_interfaces );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES (35)
  	//define node coordinates
  	deque<double> px = { 0,1,3,4.5, 0,1,3,4.5, 0,1.5, 0,2,2,3,4.5, 0,1.4,2.4,2.4,3.5,4.5, 0,1.4,2.4,2.4,3.5,4.5, 0,1.5,3,4.5, 0,1.5,3,4.5 };
    assert( px.size() == n_nodes );
  	deque<double> py = { 7,7,7,7, 5.5,5.5,5.5,5.5, 4.5,4.5, 3.5,3.5,3.5,3.5,3.5, 2.5,2.3,2.2,2.2,2.1,2, 2.5,2.3,2.2,2.2,2.1,2, 1,1,1,1, 0,0,0,0 };
    assert( py.size() == n_nodes );
  	deque<double> pz(n_nodes,0.);
  	  	  	
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------NODE BOUNDARY FLAGS (35)
    const BOX_BOUNDARY B{BOTTOM}, R{RIGHT}, U{TOP}, L{LEFT}, I{INTERNAL}, N{NOT};
    //                          0  1 2  3    4 5 6 7  8 9 1011121314 151617181920 212223242526 27282930  31  32 33  34
    vector<int8_t> bflags = { CNR4,U,U,CNR3, L,N,N,R, L,I, L,I,I,N,R, L,I,I,I,I,R, L,I,I,I,I,R, L,N,I,R, CNR1,B,B,CNR2 };
    assert( bflags.size() == n_nodes );
    vset.AddBFlags( bflags.begin(), bflags.end() );
    
    
    //--------------------------TOPOTYPE NODE FLAGS (35)
    const TOPOTYPE v{MESH_VERTEX}, i{INTERSECTION_POINT}, e{EXTERIOR_POINT}, p{PERIMETER_POINT}, l{INTERIOR_LINE}, x{EXTERIOR_LINE};
    //                        0 1 2 3  4 5 6 7  8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
    vector<int8_t> gflags = { e,x,x,e, x,v,v,x, x,p, x, l, l, v, x, x, l, i, i, l, x, x, l, i, i, l, x, x, v, p, x, e, x, x, e };
    assert( gflags.size() == n_nodes );
    vset.AddBREP_Flags( gflags.begin(), gflags.end() );


  	//--------------------------ELEMENTS
  	// define nodes per element
    deque< vector<size_t> >  plist = { {0,4,5,1}, {1,5,6,2}, {2,6,7,3}, {4,8,9,5}, {5,9,6}, {8,10,11,9}, {9,12,13,6}, {6,13,14,7}, {10,15,16,11}, {11,16,17}, // 19 elements
                                       {12,18,19,13}, {13,19,20,14}, {21,27,28,22}, {28,29,23,22}, {24,29,25}, {25,29,30,26}, {27,31,32,28}, {28,32,33,29}, {29,33,34,30},
                                       {31,32}, {32,33}, {33,34}, {34,30}, {30,26}, {20,14}, {14,7}, {7,3}, {3,2}, {2,1}, {1,0}, {0,4}, {4,8}, {8,10}, {10,15}, {21,27}, {27,31}, // 17 faces
                                       {15,16,22,21}, {16,17,23,22}, {18,19,25,24}, {19,20,26,25}, {11,9,9,12}, {17,11,12,18}, {29,23,24,29} }; // 7 interfaces
    assert( plist.size() == n_cells );
    vset.AddPlist( plist.begin(), plist.end() );

     //---------------------------------NEIGHBORS
    //define neighbors per element
    deque<vector<int64_t> > pfverts = { {L,3,1,U}, {0,4,2,U}, {1,7,R,U}, {L,5,4,0}, {6,1,3}, {L,8,I,3}, {5,10,7,4}, {6,11,R,2}, {L,I,9,5}, {I,I,8}, // element neighbors
                                        {I,I,11,6}, {10,I,R,7}, {L,16,13,I}, {12,17,I,I}, {15,I,I}, {14,18,R,I}, {L,B,17,12}, {16,B,18,13}, {17,B,R,15},
                                        // face neighbors: { face-nbors, connected high-dim elmts, local face ids of high dim elmts }
                                        {20,35,16,B,1,B}, {21,19,17,B,1,B}, {22,20,18,B,1,B},
                                        {23,21,18,R,2,R}, {24,22,15,R,2,R}, {25,23,11,R,2,R}, {26,24,7,R,2,R}, {27,25,2,R,2,R},
                                        {28,26,2,U,3,U}, {29,27,1,U,3,U}, {30,28,0,U,3,U},
                                        {31,29,0,L,0,L}, {32,30,3,L,0,L}, {33,31,5,L,0,L}, {34,32,8,L,0,L}, {35,33,12,L,0,L}, {19,34,16,L,0,L},
                                        // interface neigbhors: like faces, but with extra entry for potential index of intervening element
                                        {37,L,8,12,1,3,I}, {38,36,9,13,0,3,I}, {39,37,10,14,1,1,I}, {R,38,11,15,1,3,I},
                                        {41,I,5,6,2,0,I}, {42,40,9,10,1,0,I}, {I,41,13,14,2,2,I} };
    assert( pfverts.size() == n_cells );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );

    // creating a matching model topology
    ModelTopology mesh_topology( "SPLIT22_BASIC", true );
    
    // surface elements fall in 2 domains "lower" and "upper"
    mesh_topology.AddDomain( "lower", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE", "ISOPARAMETRIC_LINEAR_QUADRILATERAL"},
                              vector<size_t>{ 12,13,14,15,16,17,18 } );
    mesh_topology.AddDomain( "upper", set<string>{"ISOPARAMETRIC_LINEAR_TRIANGLE", "ISOPARAMETRIC_LINEAR_QUADRILATERAL"},
                              vector<size_t>{ 0,1,2,3,4,5,6,7,8,9,10,11 } );
    // boundaries
    mesh_topology.AddDomain( "BOTTOM", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{19,20,21} );
    mesh_topology.AddDomain( "RIGHT",  set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{22,23,24,25,26} );
    mesh_topology.AddDomain( "TOP",    set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{27,28,29} );
    mesh_topology.AddDomain( "LEFT",   set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{30,31,32,33,34,35} );

    // split boundaries
    mesh_topology.AddDomain( "horizontal_splitboundary", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{36,37,38,39} );
    mesh_topology.AddDomain( "inclined_split_boundary", set<string>{"ISOPARAMETRIC_LINEAR_BAR"}, vector<size_t>{40,41,42} );

    assert( mesh_topology.Cells() == vset.Elements() + vset.Faces() + vset.Interfaces() );
    assert( mesh_topology.Cells() == vset.Cells() );
    
    // node manifolds
    ManifoldType SB{ ManifoldType::SPLIT_BOUNDARY }, SBE{ ManifoldType::SPLIT_BOUNDARY_END }, SBX{ ManifoldType::SPLIT_BOUNDARY_CROSSING };
    vector<pair<vector<size_t>,ManifoldType> > node_manifolds = { {{11,12},SB}, {{15,21},SBE}, {{16,22},SB},
                                                                  {{17,18,23,24},SBX}, {{19,25},SB}, {{20,26},SBE} };
    vset.AddNodeManifolds( node_manifolds.begin(),
                           node_manifolds.end() );

    // adding corresponding materials to VSet
    const size_t n_elements{19};
    vector<int32_t> pmtrl = { 1,1,1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2 };
    assert( pmtrl.size() == n_elements );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    // adding node and element numbers for comparisons
    PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t n{0U}; n<vset.Elements(); ++n ) pushBack( elmt_nums, makeScalar( ANY, n ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t n{0U}; n<vset.Vertices(); ++n ) pushBack( node_nums, makeScalar( ANY, n ) );
    vset.AddData( "node number", node_nums );
    // permeability
    PropertyData perm( ELEMENT, SCALAR, 2U );
    perm.Reserve( vset.Elements() );
    for ( size_t n{0U}; n<mesh_topology.CellsWithinDomain("upper"); ++n ) pushBack( perm, makeScalar( ANY, 1.0e-13 ) );
    for ( size_t n{0U}; n<mesh_topology.CellsWithinDomain("lower"); ++n ) pushBack( perm, makeScalar( ANY, 1.0e-12 ) );
    vset.AddData( "permeability", perm );

    vset.Out();
    
    return mesh_topology;
    
 } // end create_BoundarySplitBoundaryPatch




void create_1Hexahedron_VSet( VSet<3U>& vset, bool bSkewed )
{
  	IsoparametricLinearHexahedron iso_hexahedron;
  	
  	//elements hexahedrons
    const size_t nodes(8);  
    vset.Resize( nodes, 6, ISOPARAMETRIC_LINEAR_HEXAHEDRON, 8, 1 );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  	
  	px[0]=0;py[0]=0;pz[0]=0;
  	px[1]=1;py[1]=0;pz[1]=0;
  	px[2]=1;py[2]=1;pz[2]=0;
  	px[3]=0;py[3]=1;pz[3]=0;
  	px[4]=0;py[4]=0;pz[4]=1;
  	px[5]=1;py[5]=0;pz[5]=1;
  	px[6]=1;py[6]=1;pz[6]=1;
  	px[7]=0;py[7]=1;pz[7]=1;
  	
  	if( bSkewed )
  	 for ( unsigned int i = 0; i < 8; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
  	//--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    deque<vector<size_t> > deqElements(1);
    deqElements[0].resize(8);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    deqElements[0][4]= 5;
    deqElements[0][5]= 6;
    deqElements[0][6]= 7;
    deqElements[0][7]= 8;
  	vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(6);
  	deqElementNeighbors[0][0]= BACK_OUTSIDE;
  	deqElementNeighbors[0][1]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][2]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][3]= TOP_OUTSIDE;
  	deqElementNeighbors[0][4]= LEFT_OUTSIDE;
  	deqElementNeighbors[0][5]= FRONT_OUTSIDE;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
    vset.ResizeBFlags();
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 0, CNR1);
  	vset.BFlag( 1, CNR2);
  	vset.BFlag( 2, CNR3);
  	vset.BFlag( 3, CNR4);
  	vset.BFlag( 4, CNR5);
  	vset.BFlag( 5, CNR6);
  	vset.BFlag( 6, CNR7);
  	vset.BFlag( 7, CNR8);
 
     //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    vset.EstablishZeroBasedNumbering();
    vset.Out();
}





/**
    @test SKM 10/7/2024 - nbor connectivity and node flags are consistent with CSMP conventions
*/
void create_Hexahedra_VSet(VSet<3U>& vset, bool bSkewed )
{
    const size_t iNrOfElements(27);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	
  	//elements hexahedrons
    size_t nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid
    
    //------------------------CREATE VSET
    //this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle 	
  	vset.Resize( iso_hexahedron.Nodes(),
                 iso_hexahedron.Neighbors(),
                 iso_hexahedron.ElementType(),
                 nodes, iNrOfElements );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(auto i{0U}; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
 // 	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    for(size_t k{0U}; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(size_t i{0U}; i < iDim_i-1; i++) //x
  	{
  	 const size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
      
  	 deqElements[iElement].resize(8);
  	 
  	 deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end() );

    //---------------------------------NEIGHBORS
    //define neighbors
    const int64_t   iDim_i_(iDim_i), iDim_j_(iDim_j), iDim_k_(iDim_k), iDim_km1_2_(iDim_km1_2);
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    for( int64_t  k = 0; k < iDim_k_-1; k++) //z
  	for( int64_t  j = 0; j < iDim_j_-1; j++) //y
  	for( int64_t  i = 0; i < iDim_i_-1; i++) //x
  	{
  	 const size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(6);
  	 
     //face 0
     deqElementNeighbors[iElement][0]= k==0?BACK_OUTSIDE:(1+ iDim_km1_2_*(k-1)+(iDim_j_-1)*j+i);
     //face 1
     deqElementNeighbors[iElement][1]= j==0?BOTTOM_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*(j-1)+i);
     //face 2
     deqElementNeighbors[iElement][2]= i==iDim_i-2?RIGHT_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*j+i+1);
     //face iDim_k-1
     deqElementNeighbors[iElement][3]= j==iDim_j-2?TOP_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*(j+1)+i);
     //face 4
     deqElementNeighbors[iElement][4]= i==0?LEFT_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*j+(i-1));
     //face 5
     deqElementNeighbors[iElement][5]= k==iDim_k-2?FRONT_OUTSIDE:(1+ iDim_km1_2_*(k+1)+(iDim_j_-1)*j+i);
      	 
    }
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(auto i{0U}; i < iDim_i; i++) //x
  	{
  	  int8_t bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR4; //shouldn't this be i?
          else if(i==(iDim_i-1)) bBoundary=CNR3; //shouldn't this be i?
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i));
        vset.BFlag( iNode, bBoundary);
      }
  	}
  	
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    vset.EstablishZeroBasedNumbering();
    vset.Out();
}











void create_Square_VSet( VSet<2U>& vset, int size_sides, double dimension, bool bSkewed )
{
  if(size_sides==1)
    create_1Square_VSet( vset, dimension, bSkewed );
  else
    create_SlitRectangle_VSet( vset, size_sides, size_sides, dimension, dimension, 0, bSkewed );
}







// TODO: nbor connectivity and node flags are inconsistent with CSMP conventions
void create_SlitRectangle_VSet( VSet<2U>& vset, int x_dimension, int y_dimension,
                                double x_length, double y_length, int depth_of_slit, bool bSkewed )
{
  assert(x_dimension>0);
  assert(y_dimension>0);
  assert(depth_of_slit<x_dimension);
  assert(x_length>1.e-7);
  assert(y_length>1.e-7);
  
  IsoparametricLinearQuadrilateral iso_quadrilateral;
  	
  //node dimensions of box
  long iDim_i(x_dimension+1);//j - width
  long iDim_j(y_dimension+1);//i - length
  	
  size_t iNrOfElements( (iDim_i-1)*(iDim_j-1) );
    
  //elements
  size_t nodes(iDim_i*iDim_j); // number of nodes
  
  cout <<"\ncreate_SlitRectangle_VSet:\n";
  cout << "\n\tDim i: " << iDim_i << " Dim j: " << iDim_j << " nr of elements: " << iNrOfElements << " nodes: " << nodes;
     
  vset.Resize( iso_quadrilateral.Nodes(),
               iso_quadrilateral.Neighbors(),
               iso_quadrilateral.ElementType(), 
               nodes, iNrOfElements );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	double delta_x = x_length/static_cast<double>(x_dimension);
  	double delta_y = y_length/static_cast<double>(y_dimension);
  	
  	for(long j{0U}; j < iDim_j; j++) //y
  	for(long i{0U}; i < iDim_i; i++) //x
  	{
  	  if(bSkewed)
  	  {
  	    px[(iDim_i)*j+i]=i*delta_x + (rand()%2000)*PERTURBATION;
  	    py[(iDim_i)*j+i]=j*delta_y + (rand()%2000)*PERTURBATION;
  	    pz[(iDim_i)*j+i]=0.;
  	  }
  	  else
  	  {
  	    px[(iDim_i)*j+i]=i*delta_x;
  	    py[(iDim_i)*j+i]=j*delta_y;
  	    pz[(iDim_i)*j+i]=0.;
  	  }
  	}
  	
    //--------------------------ELEMENTS
    //define quad elements (elements 0->26), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(size_t i{0U}; i < iDim_i-1; i++) //x
  	{
  	 const size_t iElement((iDim_i-1)*j+i);
      
  	 deqElements[iElement].resize(4);
  	 
  	 deqElements[iElement][0]= 1+ (iDim_i)*j+i;
  	 deqElements[iElement][1]= 1+ (iDim_i)*j+i+1;
  	 deqElements[iElement][2]= 1+ (iDim_i)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ (iDim_i)*(j+1)+i;
  	}
  	
    //---------------------------------NEIGHBORS
    //define neighbors
    const int64_t   iDim_i_(iDim_i), iDim_j_(iDim_j);
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    for( int64_t j = 0; j < iDim_j_-1; j++ ) //y
  	for( int64_t i = 0; i < iDim_i_-1; i++ ) //x
  	{
  	 const size_t iElement((iDim_i-1)*j+i);
     
     deqElementNeighbors[iElement].resize(4);
  	 
	 	 //face 1
     deqElementNeighbors[iElement][0]= j==0?BOTTOM_OUTSIDE:(1+ (iDim_i_-1)*(j-1)+i);
     //face 2
     deqElementNeighbors[iElement][1]= i==iDim_i_-2?RIGHT_OUTSIDE:(1+ (iDim_i_-1)*j+i+1);
     //face iDim_k-1
     deqElementNeighbors[iElement][2]= j==iDim_j_-2?TOP_OUTSIDE:(1+ (iDim_i_-1)*(j+1)+i);
     //face 4
     deqElementNeighbors[iElement][3]= i==0?LEFT_OUTSIDE:(1+ (iDim_i_-1)*j+(i-1));
     
     const bool over_slit  = (j == y_dimension/2       && i+1 >= x_dimension-depth_of_slit);
     const bool under_slit = (j == (y_dimension/2 - 1) && i+1 >= x_dimension-depth_of_slit);
     
     //cout << "\n\nTo be over the slit: j("<<j<<") == " <<  y_dimension/2 << " and i("<<i<<") > " << x_dimension-depth_of_slit;
     //cout << "\n\nTo be under the slit: j("<<j<<") == " <<  y_dimension/2-1 << " and i("<<i<<") > " << x_dimension-depth_of_slit;
     //cout << "\nElement " << iElement << "-> i,j:" << i << "," << j << " under slit? " << (under_slit?"yes":"no") << " over slit? " << (over_slit?"yes":"no");
     if (over_slit)
  	  deqElementNeighbors[iElement][1]= static_cast<int32_t>(IRREGULAR_OUTSIDE);
     else if (under_slit)
  	  deqElementNeighbors[iElement][3]= static_cast<int32_t>(IRREGULAR_OUTSIDE);
  	   
    }
    
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(size_t i{0U}; i < iDim_i; i++) //x
  	{
  	  int8_t bBoundary = NOT;
  	  
      if(j==0)
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else bBoundary=BOTTOM_OUTSIDE;
      }
      else if(j==(iDim_j-1))
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else bBoundary=TOP_OUTSIDE;
      }
      else //j in the middle
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else ;//do nothing: no boundary
      }
  	 
  	 if(bBoundary!=NOT)
     {
        const size_t iNode(((iDim_i)*j+i));
        vset.BFlag( iNode, bBoundary);
     }
     }
    
    cout << "\n\tIntroduce Slit Nodes...";
    
    //go over elements, introduce slit
    for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{ 
     const bool over_slit  = (j == y_dimension/2       && i >= x_dimension-depth_of_slit);
     const bool under_slit = (j == (y_dimension/2 - 1) && i >= x_dimension-depth_of_slit);
     
     const size_t iElement((iDim_i-1)*j+i);
    
     //cout << "\nElement " << iElement << "-> i,j:" << i << "," << j << " under slit? " << (under_slit?"yes":"no") << " over slit? " << (over_slit?"yes":"no");
     
     if (over_slit)
  	 {
  	  if(i > x_dimension-depth_of_slit) // node is NOT at the end of the slit
  	  {
    	  //create new duplicate node (0)
    	  double new_px(px[(iDim_i)*j+i]), 
    	            new_py(py[(iDim_i)*j+i]), 
    	            new_pz(pz[(iDim_i)*j+i]);
    	  size_t node_number = px.size();
    	  px.push_back(new_px);
    	  py.push_back(new_py);
    	  pz.push_back(new_pz);
    	  //set new duplicate node
        deqElements[iElement][0]= 1+ node_number;
        vset.BFlag( node_number, IRREGULAR_OUTSIDE);
      }
      
      //create new duplicate node (1)
  	  double new_px = px[(iDim_i)*j+i+1]; 
  	  double new_py = py[(iDim_i)*j+i+1]; 
  	  double  new_pz = pz[(iDim_i)*j+i+1];
  	  size_t node_number = px.size();
  	  px.push_back(new_px);
  	  py.push_back(new_py);
  	  pz.push_back(new_pz);
  	  //set new duplicate node
      deqElements[iElement][1]= 1+ node_number;
      
      //set boundary
  	  vset.BFlag( node_number, IRREGULAR_OUTSIDE);
     }
  	 
  	 if (under_slit)
       {
        //create new duplicate node (2)
        double new_px(px[(iDim_i)*(j+1)+i+1]),
                 new_py(py[(iDim_i)*(j+1)+i+1]),
                 new_pz(pz[(iDim_i)*(j+1)+i+1]);
        size_t node_number = px.size();
        px.push_back(new_px);
        py.push_back(new_py);
        pz.push_back(new_pz);
        //set new duplicate node
        deqElements[iElement][2]= 1+ node_number;
        //set boundary
        vset.BFlag( node_number, IRREGULAR_OUTSIDE);

        if(i > x_dimension-depth_of_slit) // node is NOT at the end of the slit
          {
            //create new duplicate node (3)
            double new_px2 = px[(iDim_i)*(j+1)+i];
            double new_py2 = py[(iDim_i)*(j+1)+i];
            double new_pz2 = pz[(iDim_i)*(j+1)+i];
            node_number = px.size();
            px.push_back(new_px2);
            py.push_back(new_py2);
            pz.push_back(new_pz2);
            //set new duplicate node
            deqElements[iElement][3]= 1+ node_number;
            
            //set boundary
            vset.BFlag( node_number, IRREGULAR_OUTSIDE );
          }
       }
      
      }
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.AddPlist( deqElements.begin(),deqElements.end());
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	    
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

  vset.EstablishZeroBasedNumbering();
  vset.Out();
}








/**
    3D model, which is a cube of hexahedra with six pyramid elements in the middle.
    
    Creates:
    - 32 elements (26 hex + 6 pyramids)
    - 64 nodes
    
    @test boundary flags not correct yet
*/
void create_Pyramid_Hexa_VSet(VSet<3U> & vset, bool bSkewed )
{
    const int iNrOfElements(32/*26 hexahedrons + 6 pyramids*/);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	IsoparametricLinearPyramid    iso_pyramid;
  	
  	//node dimensions of box
  	const int iDim_k(4);//k - height
  	const int iDim_j(4);//j - width
  	const int iDim_i(4);//i - length
  	const int iDim_k2(iDim_k*iDim_k);//k - height
  	const int iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	const int iPyramidsPlacement(13);
  	
  	//this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle
  	
  	//elements 0->7 are hexahedrons
  	//elements 8->13 are pyramids
    unsigned int     nodes((iDim_i*iDim_j*iDim_k)+1);  //number of nodes: 64 on a 4x4x4 grid + 1 barycenter
  	deque<uint32_t>  npes(iNrOfElements);  //number of nodes per element
    deque<uint32_t>  epes(iNrOfElements);  //element type per element
    deque<int8_t>    etypes(iNrOfElements,ISOPARAMETRIC_LINEAR_HEXAHEDRON); // NB: the pyramid elements still need to be dealt with

    for( unsigned int iElement = 0u; iElement < 26U; iElement++ )
  	{
  	  npes[iElement]=iso_hexahedron.Nodes();
  	  epes[iElement]=iso_hexahedron.Neighbors();
  	}
  	for( unsigned int iElement = 26U; iElement < 32U; iElement++ )
  	{
  	  npes[iElement]=iso_pyramid.Nodes();
  	  epes[iElement]=iso_pyramid.Neighbors();
  	}
    
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int8_t> vecElementTypes(iNrOfElements);
    for( unsigned int iElement = 0u; iElement < 26U; iElement++ )
  	{
  	  vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_HEXAHEDRON;
  	}
  	for( unsigned int iElement = 26U; iElement < 32U; iElement++ )
  	{
      vecElementTypes[iElement] = ISOPARAMETRIC_LINEAR_PYRAMID;
      etypes[iElement]          = ISOPARAMETRIC_LINEAR_PYRAMID;
  	}

    vset.Resize( etypes, npes, epes, nodes, 0, 0 );
    vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );
    
    //---------------------------MATERIALS
    const int32_t material_id(5); // some plausible integer identifier
    vector<int32_t> materials( iNrOfElements, material_id );
    vset.AddPmtrl( materials.begin(), materials.end() );

    
  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	for( unsigned int k{0U}; k < iDim_k; k++ ) //z
  	for( unsigned int j{0U}; j < iDim_j; j++ ) //y
  	for( unsigned int i{0U}; i < iDim_i; i++ ) //x
  	{
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	  
  	}
  	
  	//add barycenter 
  	px[64]=3./2. + (rand()%2000)*PERTURBATION;
  	py[64]=3./2. + (rand()%2000)*PERTURBATION;
  	pz[64]=3./2. + (rand()%2000)*PERTURBATION;
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->31), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    for( unsigned int k{0U}; k < iDim_k-1; k++ ) //z
  	for( unsigned int j{0U}; j < iDim_j-1; j++ ) //y
  	for( unsigned int i{0U}; i < iDim_i-1; i++ ) //x
  	{
  	 unsigned int iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     if(i==1 && j==1 && k==1) //its the center element (6 pyramids)
  	  continue;

     if(iElement>iPyramidsPlacement)
      iElement--;
     
  	 deqElements[iElement].resize(8);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    }
  	
  	int node_center(64);
    //element 26, assign nodes per element
    vector<size_t> vecNodes(5);
    vecNodes[0]= 1+ 21;
    vecNodes[1]= 1+ 22;
    vecNodes[2]= 1+ 26;
    vecNodes[3]= 1+ 25;
    vecNodes[4]= 1U + node_center;
    deqElements[26]=vecNodes;

    //element 27
    vecNodes[0]= 1+ 37;
    vecNodes[1]= 1+ 38;
    vecNodes[2]= 1+ 22;
    vecNodes[3]= 1+ 21;
    vecNodes[4]= 1+ node_center;
    deqElements[27]=vecNodes;

    //element 28
    vecNodes[0]= 1+ 22;
    vecNodes[1]= 1+ 38;
    vecNodes[2]= 1+ 42;
    vecNodes[3]= 1+ 26;
    vecNodes[4]= 1+ node_center;
    deqElements[28]=vecNodes;

    //element 29
    vecNodes[0]= 1+ 25;
    vecNodes[1]= 1+ 26;
    vecNodes[2]= 1+ 42;
    vecNodes[3]= 1+ 41;
    vecNodes[4]= 1+ node_center;
    deqElements[29]=vecNodes;
    
    //element 30
    vecNodes[0]= 1+ 37;
    vecNodes[1]= 1+ 21;
    vecNodes[2]= 1+ 25;
    vecNodes[3]= 1+ 41;
    vecNodes[4]= 1+ node_center;
    deqElements[30]=vecNodes;

    //element 31
    vecNodes[0]= 1+ 38;
    vecNodes[1]= 1+ 37;
    vecNodes[2]= 1+ 41;
    vecNodes[3]= 1+ 42;
    vecNodes[4]= 1+ node_center;
    deqElements[31]=vecNodes;
    
    vset.AddPlist( deqElements.begin(), deqElements.end() );

    //---------------------------------NEIGHBORS
    const int   iDim_i_(iDim_i), iDim_j_(iDim_j), iDim_k_(iDim_k);
    //define neighbors
    deque<vector<int64_t> >  deqElementNeighbors(iNrOfElements);
    for( unsigned int k{0U}; k < iDim_k_-1; k++ ) //z
  	for( unsigned int j{0U}; j < iDim_j_-1; j++ ) //y
  	for( unsigned int i{0U}; i < iDim_i_-1; i++ ) //x
  	{
  	 if(i==1 && j==1 && k==1) //its the center element (6 pyramids)
  	  continue;
  	 
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     if(iElement>iPyramidsPlacement)
      iElement--;
      
     deqElementNeighbors[iElement].resize(6);
  
     //face 0
     int64_t iNeighbor(1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i);
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;
     deqElementNeighbors[iElement][0U]= (k==0) ? BACK_OUTSIDE : iNeighbor;
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;
     deqElementNeighbors[iElement][1]= (j==0) ? BOTTOM_OUTSIDE : iNeighbor;
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][2]= (i==iDim_i_-2) ? RIGHT_OUTSIDE : iNeighbor;
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][3]= (j==iDim_j_-2) ? TOP_OUTSIDE : iNeighbor;
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+(i-1);
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][4]= (i==0) ? LEFT_OUTSIDE : iNeighbor;
     //face 5
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][5]= (k==iDim_k_-2) ? FRONT_OUTSIDE : iNeighbor;
  
    }
    
    //hexa-elements with pyramid neighbors are (6): 
    
    //element 5  - face 5 -> neighbor: 26 (under)
    deqElementNeighbors[4][5]=1+26;
    
    //element 22 - face 0 -> neighbor: 31 (over)
    deqElementNeighbors[21][0]=1+31;
    
    //element 13 - face 2 -> neighbor: 30 (left)
    deqElementNeighbors[12][2]=1+30;

    //element 14 - face 4 -> neighbor: 28 (right)
    deqElementNeighbors[13][4]=1+28;
    
    //element 11 - face 3 -> neighbor: 27 (front)
    deqElementNeighbors[10][3]=1+27;
    
    //element 16 - face 1 -> neighbor: 29 (front)
    deqElementNeighbors[15][1]=1+29;
    
    //pyramid neighbors are:
    vector<int64_t> vecNeighbors(5);
    
    //element 26
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+4;
    deqElementNeighbors[26]=vecNeighbors;

    //element 27
    vecNeighbors[0]=1+31;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+26;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+10;
    deqElementNeighbors[27]=vecNeighbors;

    //element 28
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+31;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+26;
    vecNeighbors[4]=1+13;
    deqElementNeighbors[28]=vecNeighbors;

    ///element 29
    vecNeighbors[0]=1+26;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+31;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+15;
    deqElementNeighbors[29]=vecNeighbors;

    //element 30
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+26;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+31;
    vecNeighbors[4]=1+12;
    deqElementNeighbors[30]=vecNeighbors;

    //element 31
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+30;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+28;
    vecNeighbors[4]=1+21;
    deqElementNeighbors[31]=vecNeighbors;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//------------------------------------------NODE BOUNDARY FLAGS
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
    // tested: SKM 27/5/2024
  	for(size_t k{0U}; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(size_t i{0U}; i < iDim_i; i++) //x
  	{
  	  BOX_BOUNDARY bBoundary = NOT;
  	  
     // back
  	  if(k==0)
  	  {
  	    if(j==0) // along edge1
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1)) // along edge3
  	    {
          if(i==0) bBoundary=CNR4;
          else if(i==(iDim_i-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK;
  	    }
  	  }
      // front
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5; // along edge9
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8; // along edge11
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT;
  	    }
  	  }

  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT;
          else if(i==(iDim_i-1)) bBoundary=RIGHT;
          else ;//do nothing: no boundary
  	    }
  	  }

      // SKM FIX - all flags must be captured and different indexing is required
      const size_t iNode((iDim_k2*k+(iDim_j)*j+i));
      vset.BFlag( iNode, bBoundary );
// TESTING
//      cout <<" "<< iNode <<":"<< parseBoundary( bBoundary );
  	}
  	
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    vset.EstablishZeroBasedNumbering();
    
//    vset.Out();
    
} // end create_Pyramid_Hexa_VSet









void create_1Prism_VSet(VSet<3U> & vset, bool bSkewed )
{
  	IsoparametricLinearPrism iso_prism;
  	
  	//elements
    size_t nodes(iso_prism.Nodes());         //number of nodes
  	deque<uint32_t>  npes(1);  //number of nodes per element
    deque<uint32_t>  epes(1);  //elements per element
    deque<int8_t>  etypes(1,ISOPARAMETRIC_LINEAR_PRISM);

    npes[0]=iso_prism.Nodes();
  	epes[0]=iso_prism.Neighbors();
  	
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int8_t> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_PRISM;
  	
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  	
  	px[0]=0;py[0]=0;pz[0]=0;
  	px[1]=1;py[1]=0;pz[1]=0;
  	px[2]=1;py[2]=1;pz[2]=0;
  	px[3]=0;py[3]=0;pz[3]=1;
  	px[4]=1;py[4]=0;pz[4]=1;
  	px[5]=1;py[5]=1;pz[5]=1;
  	
  	if( bSkewed )
  	 for ( size_t i = 0ul; i < 6; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
  	//--------------------------ELEMENTS
    //define prism element, assign nodes per element
    deque<vector<size_t> > deqElements(1);
    deqElements[0].resize(6);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    deqElements[0][4]= 5;
    deqElements[0][5]= 6;
  	vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(6);
  	deqElementNeighbors[0][0]= BACK_OUTSIDE;
  	deqElementNeighbors[0][1]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][2]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][3]= TOP_OUTSIDE/*or LEFT_OUTSIDE*/;
  	deqElementNeighbors[0][4]= FRONT_OUTSIDE;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 1, CNR1);
  	vset.BFlag( 2, CNR2);
  	vset.BFlag( 3, CNR3);
  	vset.BFlag( 4, CNR5);
  	vset.BFlag( 5, CNR6);
  	vset.BFlag( 6, CNR7);
    
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    vset.EstablishZeroBasedNumbering();
    vset.Out();
    
} // end create_1Prism_VSet






/**
    Generates 24 hexahedra + 6 prism elements.
    The model can be distorted on demand.
    
    @note model comes with the correct box boundary flags.
    
    TODO: element neighbor connectivity does not match, nodes-per-face convention of CSMP (9/7/24)
*/
void create_Prism_Hexa_VSet( VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(30/*24 hexahedrons + 6 prisms*/);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	IsoparametricLinearPrism iso_prism;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	
  	//this is a 3D model, it is a cube of hexahedra with prisms in the middle
  	
  	//elements 0->23 are hexahedrons
  	//elements 8->29 are prisms
    const size_t   nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid
  	deque<uint32_t>  npes(iNrOfElements);  //number of nodes per element
    deque<uint32_t>  epes(iNrOfElements);  //element type per element
    deque<int8_t>  etypes(iNrOfElements,ISOPARAMETRIC_LINEAR_HEXAHEDRON);
    
    size_t iElement = 0;
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{
  	  if(i==1 && j==1) //its a prism
  	  { 
  	    npes[iElement]=iso_prism.Nodes();
  	    epes[iElement]=iso_prism.Neighbors();
  	    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    
  	    npes[iElement]=iso_prism.Nodes();
  	    epes[iElement]=iso_prism.Neighbors();
  	    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    continue; 
  	  }
  	
  	  npes[iElement]=iso_hexahedron.Nodes();
  	  epes[iElement]=iso_hexahedron.Neighbors();
  	  
  	  iElement++;
  	}
  	cout << "e:" << iElement;
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
    
    //--------------------------ELEMENT TYPES
  	//add element types
  	iElement = 0;
  
    vector<int8_t> vecElementTypes(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{
  	  if(i==1 && j==1) //its a prism
  	  {
  	    vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    
  	    vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    continue; 
  	  }
  	  
  	  vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_HEXAHEDRON;
  	  iElement++;
  	}
    vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

    
  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	for(size_t k{0UL}; k < iDim_k; k++) //z
  	for(size_t j{0UL}; j < iDim_j; j++) //y
  	for(size_t i{0UL}; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    
    iElement = 0;
    for(size_t k{0UL}; k < iDim_k-1; k++) //z
  	for(size_t j{0UL}; j < iDim_j-1; j++) //y
  	for(auto i{0UL}; i < iDim_i-1; i++) //x
  	{
  	 
     //size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     if(i==1 && j==1) //its the center element (6 pyramids)
  	 {
  	  const size_t node1 = 1+ iDim_k2*k+(iDim_j)*j+i;
  	  const size_t node2 = 1+ iDim_k2*k+(iDim_j)*j+i+1;;
  	  const size_t node3 = 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	  const size_t node4 = 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	  const size_t node5 = 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	  const size_t node6 = 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	  const size_t node7 = 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	  const size_t node8 = 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	  
  	  deqElements[iElement].resize(6);
	    deqElements[iElement][0]= node1;
  	  deqElements[iElement][1]= node2;
  	  deqElements[iElement][2]= node4;
  	  deqElements[iElement][3]= node5;
  	  deqElements[iElement][4]= node6;
  	  deqElements[iElement][5]= node8;
  	  
  	  iElement++;
      
      deqElements[iElement].resize(6);
	    deqElements[iElement][0]= node3;
  	  deqElements[iElement][1]= node4;
  	  deqElements[iElement][2]= node2;
  	  deqElements[iElement][3]= node7;
  	  deqElements[iElement][4]= node8;
  	  deqElements[iElement][5]= node6;
  	  
  	  iElement++;
      
  	  continue;
  	 }
 
  	 deqElements[iElement].resize(8);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    
     iElement++;
    
    }
    
    vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    iElement = 0;
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    deqElementNeighbors[0].resize(6);
    deqElementNeighbors[0][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[0][1] = FRONT_OUTSIDE;
    deqElementNeighbors[0][2] = 1+1;
    deqElementNeighbors[0][3] = 1+3;
    deqElementNeighbors[0][4] = LEFT_OUTSIDE;
    deqElementNeighbors[0][5] = 1+10;
    
    deqElementNeighbors[1].resize(6);
    deqElementNeighbors[1][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[1][1] = FRONT_OUTSIDE;
    deqElementNeighbors[1][2] = 1+2;
    deqElementNeighbors[1][3] = 1+4;
    deqElementNeighbors[1][4] = 1+0;
    deqElementNeighbors[1][5] = 1+11;
    
    deqElementNeighbors[2].resize(6);
    deqElementNeighbors[2][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[2][1] = FRONT_OUTSIDE;
    deqElementNeighbors[2][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[2][3] = 1+6;
    deqElementNeighbors[2][4] = 1+1;
    deqElementNeighbors[2][5] = 1+12;
    
    deqElementNeighbors[3].resize(6);
    deqElementNeighbors[3][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[3][1] = 1+0;
    deqElementNeighbors[3][2] = 1+4;
    deqElementNeighbors[3][3] = 1+7;
    deqElementNeighbors[3][4] = LEFT_OUTSIDE;
    deqElementNeighbors[3][5] = 1+13;
    
    deqElementNeighbors[6].resize(6);
    deqElementNeighbors[6][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[6][1] = 1+2;
    deqElementNeighbors[6][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[6][3] = 1+9;
    deqElementNeighbors[6][4] = 1+5;
    deqElementNeighbors[6][5] = 1+16;
    
    deqElementNeighbors[7].resize(6);
    deqElementNeighbors[7][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[7][1] = 1+3;
    deqElementNeighbors[7][2] = 1+8;
    deqElementNeighbors[7][3] = BACK_OUTSIDE;
    deqElementNeighbors[7][4] = LEFT_OUTSIDE;
    deqElementNeighbors[7][5] = 1+17;
    
    deqElementNeighbors[8].resize(6);
    deqElementNeighbors[8][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[8][1] = 1+5;
    deqElementNeighbors[8][2] = 1+9;
    deqElementNeighbors[8][3] = BACK_OUTSIDE;
    deqElementNeighbors[8][4] = 1+7;
    deqElementNeighbors[8][5] = 1+18;
    
    deqElementNeighbors[9].resize(6);
    deqElementNeighbors[9][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[9][1] = 1+6;
    deqElementNeighbors[9][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[9][3] = BACK_OUTSIDE;
    deqElementNeighbors[9][4] = 1+8;
    deqElementNeighbors[9][5] = 1+19;
    
    deqElementNeighbors[10].resize(6);
    deqElementNeighbors[10][0] = 1+0;
    deqElementNeighbors[10][1] = FRONT_OUTSIDE;
    deqElementNeighbors[10][2] = 1+11;
    deqElementNeighbors[10][3] = 1+13;
    deqElementNeighbors[10][4] = LEFT_OUTSIDE;
    deqElementNeighbors[10][5] = 1+20;
    
    deqElementNeighbors[11].resize(6);
    deqElementNeighbors[11][0] = 1+1;
    deqElementNeighbors[11][1] = FRONT_OUTSIDE;
    deqElementNeighbors[11][2] = 1+12;
    deqElementNeighbors[11][3] = 1+14;
    deqElementNeighbors[11][4] = 1+10;
    deqElementNeighbors[11][5] = 1+21;
    
    deqElementNeighbors[12].resize(6);
    deqElementNeighbors[12][0] = 1+2;
    deqElementNeighbors[12][1] = FRONT_OUTSIDE;
    deqElementNeighbors[12][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[12][3] = 1+16;
    deqElementNeighbors[12][4] = 1+11;
    deqElementNeighbors[12][5] = 1+22;
    
    deqElementNeighbors[13].resize(6);
    deqElementNeighbors[13][0] = 1+3;
    deqElementNeighbors[13][1] = 1+10;
    deqElementNeighbors[13][2] = 1+14;
    deqElementNeighbors[13][3] = 1+17;
    deqElementNeighbors[13][4] = LEFT_OUTSIDE;
    deqElementNeighbors[13][5] = 1+23;
    
    deqElementNeighbors[16].resize(6);
    deqElementNeighbors[16][0] = 1+6;
    deqElementNeighbors[16][1] = 1+12;
    deqElementNeighbors[16][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[16][3] = 1+19;
    deqElementNeighbors[16][4] = 1+15;
    deqElementNeighbors[16][5] = 1+26;
    
    deqElementNeighbors[17].resize(6);
    deqElementNeighbors[17][0] = 1+7;
    deqElementNeighbors[17][1] = 1+13;
    deqElementNeighbors[17][2] = 1+18;
    deqElementNeighbors[17][3] = BACK_OUTSIDE;
    deqElementNeighbors[17][4] = LEFT_OUTSIDE;
    deqElementNeighbors[17][5] = 1+27;
    
    deqElementNeighbors[18].resize(6);
    deqElementNeighbors[18][0] = 1+8;
    deqElementNeighbors[18][1] = 1+15;
    deqElementNeighbors[18][2] = 1+19;
    deqElementNeighbors[18][3] = BACK_OUTSIDE;
    deqElementNeighbors[18][4] = 1+17;
    deqElementNeighbors[18][5] = 1+28;
    
    deqElementNeighbors[19].resize(6);
    deqElementNeighbors[19][0] = 1+9;
    deqElementNeighbors[19][1] = 1+16;
    deqElementNeighbors[19][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[19][3] = BACK_OUTSIDE;
    deqElementNeighbors[19][4] = 1+18;
    deqElementNeighbors[19][5] = 1+29;
    
    deqElementNeighbors[20].resize(6);
    deqElementNeighbors[20][0] = 1+10;
    deqElementNeighbors[20][1] = FRONT_OUTSIDE;
    deqElementNeighbors[20][2] = 1+21;
    deqElementNeighbors[20][3] = 1+23;
    deqElementNeighbors[20][4] = LEFT_OUTSIDE;
    deqElementNeighbors[20][5] = TOP_OUTSIDE;
      
    deqElementNeighbors[21].resize(6);
    deqElementNeighbors[21][0] = 1+11;
    deqElementNeighbors[21][1] = FRONT_OUTSIDE;
    deqElementNeighbors[21][2] = 1+22;
    deqElementNeighbors[21][3] = 1+24;
    deqElementNeighbors[21][4] = 1+20;
    deqElementNeighbors[21][5] = TOP_OUTSIDE;
      
    deqElementNeighbors[22].resize(6);
    deqElementNeighbors[22][0] = 1+12;
    deqElementNeighbors[22][1] = FRONT_OUTSIDE;
    deqElementNeighbors[22][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[22][3] = 1+26;
    deqElementNeighbors[22][4] = 1+21;
    deqElementNeighbors[22][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[23].resize(6);
    deqElementNeighbors[23][0] = 1+13;
    deqElementNeighbors[23][1] = 1+20;
    deqElementNeighbors[23][2] = 1+24;
    deqElementNeighbors[23][3] = 1+27;
    deqElementNeighbors[23][4] = LEFT_OUTSIDE;
    deqElementNeighbors[23][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[26].resize(6);
    deqElementNeighbors[26][0] = 1+16;
    deqElementNeighbors[26][1] = 1+22;
    deqElementNeighbors[26][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[26][3] = 1+29;
    deqElementNeighbors[26][4] = 1+25;
    deqElementNeighbors[26][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[27].resize(6);
    deqElementNeighbors[27][0] = 1+17;
    deqElementNeighbors[27][1] = 1+23;
    deqElementNeighbors[27][2] = 1+28;
    deqElementNeighbors[27][3] = BACK_OUTSIDE;
    deqElementNeighbors[27][4] = LEFT_OUTSIDE;
    deqElementNeighbors[27][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[28].resize(6);
    deqElementNeighbors[28][0] = 1+18;
    deqElementNeighbors[28][1] = 1+25;
    deqElementNeighbors[28][2] = 1+29;
    deqElementNeighbors[28][3] = BACK_OUTSIDE;
    deqElementNeighbors[28][4] = 1+27;
    deqElementNeighbors[28][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[29].resize(6);
    deqElementNeighbors[29][0] = 1+19;
    deqElementNeighbors[29][1] = 1+26;
    deqElementNeighbors[29][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[29][3] = BACK_OUTSIDE;
    deqElementNeighbors[29][4] = 1+28;
    deqElementNeighbors[29][5] = TOP_OUTSIDE;
    
    //prism neighbors are:
    vector<int64_t> vecNeighbors(5);
    
    //element 4
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=1+1;
    vecNeighbors[2]=1+5;
    vecNeighbors[3]=1+3;
    vecNeighbors[4]=1+14;
    deqElementNeighbors[4]=vecNeighbors;
  
    //element 5
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=1+8;
    vecNeighbors[2]=1+4;
    vecNeighbors[3]=1+6;
    vecNeighbors[4]=1+15;
    deqElementNeighbors[5]=vecNeighbors;

    //element 14
    vecNeighbors[0]=1+4;
    vecNeighbors[1]=1+11;
    vecNeighbors[2]=1+15;
    vecNeighbors[3]=1+13;
    vecNeighbors[4]=1+24;
    deqElementNeighbors[14]=vecNeighbors;

    //element 15
    vecNeighbors[0]=1+5;
    vecNeighbors[1]=1+18;
    vecNeighbors[2]=1+14;
    vecNeighbors[3]=1+16;
    vecNeighbors[4]=1+25;
    deqElementNeighbors[15]=vecNeighbors;

    //element 24
    vecNeighbors[0]=1+14;
    vecNeighbors[1]=1+21;
    vecNeighbors[2]=1+25;
    vecNeighbors[3]=1+23;
    vecNeighbors[4]=TOP_OUTSIDE;
    deqElementNeighbors[24]=vecNeighbors;

    //element 25
    vecNeighbors[0]=1+15;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+24;
    vecNeighbors[3]=1+26;
    vecNeighbors[4]=TOP_OUTSIDE;
    deqElementNeighbors[25]=vecNeighbors;
 
     assert( deqElementNeighbors.size() == deqElementNeighbors.size() );
     vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
 
    
  	//------------------------------------------NODE BOUNDARY FLAGS
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
    // retested: SKM: 27/5/24
  	for(size_t k{0U}; k < iDim_k; k++) //z -> Y in CSMP
  	for(size_t j{0U}; j < iDim_j; j++) //y -> X
  	for(size_t i{0U}; i < iDim_i; i++) //x -> Z
  	{
  	  int8_t bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0) // along edge1
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1)) // along edge3
  	    {
          if(i==0) bBoundary=CNR4;
          else if(i==(iDim_i-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      // SKM_FIX -1
      const size_t iNode((iDim_k2*k+(iDim_j)*j+i));
      vset.BFlag( iNode, bBoundary );
  	}

    vset.EstablishZeroBasedNumbering();
    vset.ResizeNodes( 64 );
    
    // additional must haves
    vector<int32_t> pmtrl(vset.Elements(),1); // all the same material=1
    // fill( next(pmtrl.begin(),42), pmtrl.end(), 7 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
            // element number
    PropertyData elmt_nums( ELEMENT, SCALAR, 3U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, i ) );
    vset.AddData( "element number", elmt_nums );
    
    cout <<"\n"<<"create_Prism_Hexa_VSet: model 'Prism_Hexa': "<< endl;
    //vset.Out();

} // end create_Prism_Hexa_VSet


/* FOR WHEN THE XZ PLANE IS THE BASE PLANE OF THE HEXAHEDRON

    // SKM adjustment according to CSMP_FEM_conventions.pdf (27/5/24)
    deque<vector<int64_t> >    // first vertical plane (z=0)
                               // --------------------------
    deqElementNeighborsSKM = { {BOTTOM, 10, 1, BACK, LEFT, 3},   //  0
                               {BOTTOM, 11, 2, BACK, 0, 4},      //  1
                               {BOTTOM, 12, RIGHT, BACK, 1, 6},  //  2
                               // row 2(back)
                               {0,      13, 14, BACK, LEFT, 7},  //  3
                               // prisms (elemts 4, 5)
                               {BACK, 3, 1, 5, 14},              //  4
                               {BACK, 4, 6, 8, 15},              //  5
                               // hexahedra
                               {2,  16, RIGHT, BACK, 5, 9},      //  6
                               // top row
                               {3, 17, 8, BACK, LEFT, TOP},      //  7
                               {5, 18, 9, BACK, 7, TOP},         //  8
                               {6, 19, RIGHT, BACK, 8, TOP},     //  9  (10 elmts because one hex is plit into 2 prisms)
                               // ---------------------------
                               // second plane (middle layer)
                               // ---------------------------
                               {BOTTOM, 20, 11, 0, LEFT, 13},    // 10
                               {BOTTOM, 21, 12, 1, 10, 14},      // 11
                               {BOTTOM, 22, RIGHT, 2, 11, 16},   // 12
                               // middle row
                               {10, 23, 14, 3, LEFT, 17},        // 13
                               // prims 14, 15
                               {4, 13, 11, 15, 24},              // 14
                               {5, 14, 16, 18, 25},              // 14
                               // ----------------
                               {12, 26, RIGHT, 6, 15, 19},       // 16
                               // top row (middle layer)
                               {13, 27, 18, 7, LEFT, TOP},       // 17
                               {15, 28, 19, 8, 17, TOP},         // 18
                               {16, 29, RIGHT, 9, 18, TOP},      // 19
                               // -------------------
                               // third plane (front)
                               // -------------------
                               {BOTTOM, FRONT, 21, 10, LEFT, 23},
                               {BOTTOM, FRONT, 22, 11, 20, 24},
                               {BOTTOM, FRONT, RIGHT, 12, 21, 26},
                               // row 2
                               {20, FRONT, 24, 13, LEFT, 27},
                               // prism elements 24, 25
                               {14, 23, 21, 25, FRONT},
                               {15, 24, 26, 28, FRONT},
                               // -------------------
                               {22, FRONT, RIGHT, 16, 25, 29},
                               // top row
                               {23, FRONT, 28, 17, LEFT, TOP},
                               {25, FRONT, 29, 18, 27, TOP},
                               {26, FRONT, RIGHT, 19, 28, TOP} };
*/

// SKM attempted adjustment to CSMP_FEM_conventions.pdf (27/5/24)
// (assuming that the base plane of the hex is the xz plane)
 /*
    deque<vector<int64_t> >    // first vertical plane (z=0)
                               // --------------------------
    deqElementNeighborsSKM = { {BOTTOM, 10, 1, BACK, LEFT, 3},   //  0
                               {BOTTOM, 11, 2, BACK, 0, 4},      //  1
                               {BOTTOM, 12, RIGHT, BACK, 1, 6},  //  2
                               // row 2(back)
                               {0,      13, 14, BACK, LEFT, 7},  //  3
                               // prisms (elemts 4, 5)
                               {BACK, 3, 1, 5, 14},              //  4
                               {BACK, 4, 6, 8, 15},              //  5
                               // hexahedra
                               {2,  16, RIGHT, BACK, 5, 9},      //  6
                               // top row
                               {3, 17, 8, BACK, LEFT, TOP},      //  7
                               {5, 18, 9, BACK, 7, TOP},         //  8
                               {6, 19, RIGHT, BACK, 8, TOP},     //  9  (10 elmts because one hex is plit into 2 prisms)
                               // ---------------------------
                               // second plane (middle layer)
                               // ---------------------------
                               {BOTTOM, 20, 11, 0, LEFT, 13},    // 10
                               {BOTTOM, 21, 12, 1, 10, 14},      // 11
                               {BOTTOM, 22, RIGHT, 2, 11, 16},   // 12
                               // middle row
                               {10, 23, 14, 3, LEFT, 17},        // 13
                               // prims 14, 15
                               {4, 13, 11, 15, 24},              // 14
                               {5, 14, 16, 18, 25},              // 14
                               // ----------------
                               {12, 26, RIGHT, 6, 15, 19},       // 16
                               // top row (middle layer)
                               {13, 27, 18, 7, LEFT, TOP},       // 17
                               {15, 28, 19, 8, 17, TOP},         // 18
                               {16, 29, RIGHT, 9, 18, TOP},      // 19
                               // -------------------
                               // third plane (front)
                               // -------------------
                               {BOTTOM, FRONT, 21, 10, LEFT, 23},
                               {BOTTOM, FRONT, 22, 11, 20, 24},
                               {BOTTOM, FRONT, RIGHT, 12, 21, 26},
                               // row 2
                               {20, FRONT, 24, 13, LEFT, 27},
                               // prism elements 24, 25
                               {14, 23, 21, 25, FRONT},
                               {15, 24, 26, 28, FRONT},
                               // -------------------
                               {22, FRONT, RIGHT, 16, 25, 29},
                               // top row
                               {23, FRONT, 28, 17, LEFT, TOP},
                               {25, FRONT, 29, 18, 27, TOP},
                               {26, FRONT, RIGHT, 19, 28, TOP} };
                               
    // modify neighbor sequence because in Adriana's labelling the XY plane is the base plane of the hex
    for ( auto& eit : deqElementNeighborsSKM ) {
         // copy current entry
         vector<int64_t> swapvec = eit;
         // if it refers to hex, write it out in new order
         if ( swapvec.size() == 6 ) {
              eit[0] = swapvec[3];
              eit[1] = swapvec[0];
              // eit[2] = no change required
              eit[3] = swapvec[5];
              // eit[4] = no change required
              eit[5] = swapvec[1];
           }
      }

    assert( deqElementNeighborsSKM.size() == deqElementNeighbors.size() );
    vset.AddPfverts( deqElementNeighborsSKM.begin(), deqElementNeighborsSKM.end() );
  	
*/








void create_Prism_VSet(VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(54/*prisms*/);
    
  	IsoparametricLinearPrism iso_prism;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	
  	//this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle
  	
  	//elements 0->53 are prisms
    size_t nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid

    vset.Resize( iso_prism.Nodes(),
                 iso_prism.Neighbors(),
                 iso_prism.ElementType(), 
                 nodes, iNrOfElements );
    
  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(auto i{0U}; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
  	  if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
    
    //--------------------------ELEMENTS
    vector<int8_t> vecElementTypes(1,ISOPARAMETRIC_LINEAR_PRISM);
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

    //define prism elements (elements 0->54), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{
  	 //lower element
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
      
  	 deqElements[iElement].resize(6);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;

     //upper element
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i+27;

  	 deqElements[iElement].resize(6);
  	 deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 
    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    for(int64_t  k = 0U; k < iDim_k-1; k++) //z
  	for(int64_t  j{0U}; j < iDim_j-1; j++) //y
  	for(int64_t  i{0U}; i < iDim_i-1; i++) //x
  	{
  	 //lower element
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(5);
  
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i);
     deqElementNeighbors[iElement][0]= (k==0?BACK_OUTSIDE:iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i +27;
     deqElementNeighbors[iElement][1]= (j==0?BOTTOM_OUTSIDE:iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1 +27;
     deqElementNeighbors[iElement][2]= (i==iDim_i-2?RIGHT_OUTSIDE:iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+27;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
     //upper element
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i+27;
     
     deqElementNeighbors[iElement].resize(5);

     //face 0
     iNeighbor = 1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i +27;
     deqElementNeighbors[iElement][0]= (k==0?BACK_OUTSIDE:iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i;
     deqElementNeighbors[iElement][1]= (j==iDim_j-2?TOP_OUTSIDE:iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+(i-1);
     deqElementNeighbors[iElement][2]= (i==0?LEFT_OUTSIDE:iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i +27;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
     
    }
    
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(auto i{0U}; i < iDim_i; i++) //x
  	{
  	  int8_t bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        // SKM FIX: const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i));
        vset.BFlag( iNode, bBoundary);
      }
  	}
  	
    vset.EstablishZeroBasedNumbering();

    // SKM_FIX add ons
    // ---------------
    // adding corresponding materials to VSet
    vector<int32_t> pmtrl(vset.Elements(),1); // matrix
    // fill( next(pmtrl.begin(),42), pmtrl.end(), 7 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

//    vset.Out();

} // end







/**
      Creates cube of 27 unit hexahedra, numbered from the back XY plane to the front.
      Element number is increasing with x in the rows, and from the bottom to the top, then with the planes from the back to the front.
      
      @author SKM
      @date 31/5/2024
*/
// TODO: nbor connectivity and node flags are inconsistent with CSMP conventions
void create_RubikCube( VSet<3U>& vset )
  {
    // Rubik cube 3 x 3 x 3, starting element numbering from the origin in the back plane (XY), moving left to right, from bottom to top
    deque<double>  px{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                       0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                       0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                       0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
            
    deque<double>  py{ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,  // Bottom Layer (z = 0)
                       0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,  // Middle Layer (z = 1)
                       0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3,  // Top Layer (z = 2)
                       0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 };  // Topmost Layer (z = 3)

    deque<double>  pz{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // Bottom Layer (z = 0)
                       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // Middle Layer (z = 1)
                       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  // Top Layer (z = 2)
                       3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };  // Topmost Layer (z = 3)
        
    vector<int8_t> bflags{ // xy plane in the back
                          CNR1, EDGE1, EDGE1, CNR2,
                          EDGE4, NOT, NOT, EDGE2,
                          EDGE4, NOT, NOT, EDGE2,
                          CNR4, EDGE3, EDGE3, CNR3,
                          // in midde layer (back)
                          EDGE5, BOTTOM, BOTTOM, EDGE6,
                          LEFT, NOT, NOT, RIGHT,
                          LEFT, NOT, NOT, RIGHT,
                          EDGE8, TOP, TOP, EDGE7,
                          // middle layer (front)
                          EDGE5, BOTTOM, BOTTOM, EDGE6,
                          LEFT, NOT, NOT, RIGHT,
                          LEFT, NOT, NOT, RIGHT,
                          EDGE8, TOP, TOP, EDGE7,
                          // front (z-max)
                          CNR5, EDGE9, EDGE9, CNR6,
                          EDGE12, NOT, NOT, EDGE11,
                          EDGE12, NOT, NOT, EDGE11,
                          CNR8, EDGE11, EDGE11, CNR7 };

    assert( bflags.size() == px.size() ); // "boundary flag vector has a different size than node coordinate vector"

    vector<int8_t> vecElementTypes(1,ISOPARAMETRIC_LINEAR_HEXAHEDRON);

    deque<vector<size_t>> plist{ {0, 1, 4, 5, 16, 17, 20, 21}, {1, 2, 5, 6, 17, 18, 21, 22}, {2, 3, 6, 7, 18, 19, 22, 23},
                                 {4, 5, 8, 9, 20, 21, 24, 25}, {5, 6, 9, 10, 21, 22, 25, 26}, {6, 7, 10, 11, 22, 23, 26, 27},
                                 {8, 9, 12, 13, 24, 25, 28, 29}, {9, 10, 13, 14, 25, 26, 29, 30}, {10, 11, 14, 15, 26, 27, 30, 31},
                                 {16, 17, 20, 21, 32, 33, 36, 37}, {17, 18, 21, 22, 33, 34, 37, 38}, {18, 19, 22, 23, 34, 35, 38, 39},
                                 {20, 21, 24, 25, 36, 37, 40, 41}, {21, 22, 25, 26, 37, 38, 41, 42}, {22, 23, 26, 27, 38, 39, 42, 43},
                                 {24, 25, 28, 29, 40, 41, 44, 45}, {25, 26, 29, 30, 41, 42, 45, 46}, {26, 27, 30, 31, 42, 43, 46, 47},
                                 {32, 33, 36, 37, 48, 49, 52, 53}, {33, 34, 37, 38, 49, 50, 53, 54}, {34, 35, 38, 39, 50, 51, 54, 55},
                                 {36, 37, 40, 41, 52, 53, 56, 57}, {37, 38, 41, 42, 53, 54, 57, 58}, {38, 39, 42, 43, 54, 55, 58, 59},
                                 {40, 41, 44, 45, 56, 57, 60, 61}, {41, 42, 45, 46, 57, 58, 61, 62}, {42, 43, 46, 47, 58, 59, 62, 63} };

    const auto BA{BACK}, BO{BOTTOM}, L{LEFT}, R{RIGHT}, T{TOP}, F{FRONT};
    deque<vector<int64_t>> pfverts = { {BO,9,1,BA,L,3}, {BO,10,2,BA,0,4}, {BO,11,R,BA,1,5}, // backplane
                                       {0,12,4,BA,L,6}, {1,13,5,BA,3,7}, {2,14,R,BA,4,8},
                                       {3,15,7,BA,L,T}, {4,16,8,BA,6,T}, {5,17,R,BA,7,T},
                                       {BO,18,10,0,L,12}, {BO,19,11,1,9,13}, {BO,20,R,2,10,14}, // middle plane
                                       {9,21,13,3,L,15}, {10,22,14,4,12,16}, {11,23,R,5,13,17},
                                       {12,24,16,6,L,T}, {13,25,17,7,15,T}, {14,26,R,8,16,T},
                                       {9,F,19,0,L,21}, {10,F,20,1,18,22}, {11,F,R,2,19,23}, // front plane
                                       {12,F,22,3,L,24}, {13,F,23,4,21,25}, {14,F,R,5,22,26},
                                       {15,F,25,6,L,T}, {16,F,26,7,24,T}, {17,F,R,8,25,T} };

    const size_t n_elements{plist.size()}, n_nodes{px.size()};
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	
    //------------------------CREATE VSET
    //this is a 3D model, it is a cube of unit-cell hexahedra
  	vset.Resize( iso_hexahedron.Nodes(),
                 iso_hexahedron.Neighbors(),
                 iso_hexahedron.ElementType(),
                 n_nodes, n_elements );

 	  vset.AddXYZ( px, py, pz );
    vset.AddBFlags( bflags.begin(), bflags.end() );

  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );
  	vset.AddPlist( plist.begin(), plist.end() );
    vset.AddPfverts( pfverts.begin(), pfverts.end());

    // all elements are of the same material labelled 1
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    vset.Out();

} // end create_RubikCube





/**
       Decomposition of a hexahedron into 6 tetrahedra.
       Only 6 elements!
       
       @test is consistent with CSMP UG, SKM 3/6/2024
*/
void create_Tetra_VSet( VSet<3U>& vset )
 {
    const size_t iNrOfElements{6}, nodes{8};
  	IsoparametricLinearTetrahedron  iso_tet;

    vset.Resize( iso_tet.Nodes(),
                 iso_tet.Neighbors(),
                 iso_tet.ElementType(),
                 nodes, iNrOfElements );

  	//-----------------------NODES
  	//define nodes
    const size_t n_nodes{8};
  	deque<double> px(n_nodes);
  	deque<double> py(n_nodes);
  	deque<double> pz(n_nodes);

    // node coordinates
    px[0] =-1.0;
    px[1] = 1.0;
    px[2] = 1.0;
    px[3] =-1.0;
    px[4] =-1.0;
    px[5] = 1.0;
    px[6] = 1.0;
    px[7] =-1.0;
    // nodal y-coordinates
    py[0] =-1.0;
    py[1] =-1.0;
    py[2] = 1.0;
    py[3] = 1.0;
    py[4] =-1.0;
    py[5] =-1.0;
    py[6] = 1.0;
    py[7] = 1.0;
    // nodal z-coordinates
    pz[0] = -1.0;
    pz[1] = -1.0;
    pz[2] = -1.0;
    pz[3] = -1.0;
    pz[4] =  1.0;
    pz[5] =  1.0;
    pz[6] =  1.0;
    pz[7] =  1.0;

   	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
    
    // SKM 3/6/24 - made compliant with CSMP UG boundary flags (fig.5)
    vset.BFlag(0,CNR1);
    vset.BFlag(1,CNR2);
    vset.BFlag(2,CNR3);
    vset.BFlag(3,CNR4);
    vset.BFlag(4,CNR5);
    vset.BFlag(5,CNR6);
    vset.BFlag(6,CNR7);
    vset.BFlag(7,CNR8);


    // -------------------------PELMT
    vector<int8_t> vecElementTypes(1,ISOPARAMETRIC_LINEAR_TETRAHEDRON);
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );
   

    //--------------------------ELEMENTS ('plist')
    //define tetrahedral elements (1..6), assign nodes per element (see hexa decomposition)
    deque<vector<size_t> >  deqElements{ {0,1,3,4}, {4,5,1,3}, {4,5,3,7}, {1,2,3,6}, {3,1,6,5}, {5,6,3,7} };

  	vset.AddPlist( deqElements.begin(),deqElements.end());


    //--------------------------ELEMENT NEIGHBORS          0                   1                  2
    deque<vector<int64_t> >  deqElementNeighbors{ {1,LEFT,BOTTOM,BACK}, {4,0,2,BOTTOM}, {5,LEFT,FRONT,1},
                                                 //        3                 4              5
                                                  {TOP,4,RIGHT,BACK}, {RIGHT,5,1,3}, {TOP,2,FRONT,4} };

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());

    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    cout <<"\n"<<"testCreateTetra_VSet: model 'Tetra': "<< endl;
    vset.Out();

 } // end testCreateTetra_VSet







void create_Pyramid_VSet( VSet<3U>& vset, bool bSkewed )
{
  	IsoparametricLinearPyramid iso_pyramid;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
    
    const size_t iNrOfCells((iDim_i-1)*(iDim_j-1)*(iDim_k-1));
    const size_t iNrOfElements(iNrOfCells*6U/*pyramids*/);
     	
  	//this is a 3D model, it is a cube of 27 sets of six pyramids glued by the apex
  	
  	//elements 0->161 are pyramids
    size_t nodes(iDim_i*iDim_j*iDim_k+iNrOfCells);  //number of nodes: 64 on a 4x4x4 grid + one barycenter for each cell

    vset.Resize( iso_pyramid.Nodes(),
                 iso_pyramid.Neighbors(),
                 iso_pyramid.ElementType(),
                 nodes, iNrOfElements );
    
  	//-----------------------NODES
  	//define nodes
  	deque<double> px(nodes);
  	deque<double> py(nodes);
  	deque<double> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(size_t i{0U}; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//add barycenters, one barycenter per cell
  	for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
  	  //const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
  	  
  	  if(bSkewed)
  	  {
    	  px[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=i + 1./2. + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=j + 1./2. + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=k + 1./2. + (rand()%2000)*PERTURBATION;
      }
      else
      {
        px[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=i + 1./2.;
  	    py[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=j + 1./2.;
  	    pz[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=k + 1./2.;
      }
    }  	
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
    //--------------------------ELEMENTS
    //define pyramid elements (elements 0->54), assign nodes per element
    deque<vector<size_t> > deqElements(iNrOfElements);
    for(size_t k{0U}; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(size_t i{0U}; i < iDim_i-1; i++) //x
  	{
  	 //we have 6 pyramids per cell
  	 
  	 //pyramid 0
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 1
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*1;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;
  	 
  	 //pyramid 2
  	 cout << "Element:(i="<<i<<",j="<<j<<",k="<<k<<")" << iDim_km1_2*k+(iDim_j-1)*j+i + 27*2 << endl;
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*2;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
	   cout << "Coord [" << iElement << "][0]: " << deqElements[iElement][0] << endl;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
	   cout << "Coord [" << iElement << "][1]: " << deqElements[iElement][1] << endl;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
	   cout << "Coord [" << iElement << "][2]: " << deqElements[iElement][2] << endl;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   cout << "Coord [" << iElement << "][3]: " << deqElements[iElement][3] << endl;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;
	   cout << "Coord [" << iElement << "][4]: " << deqElements[iElement][4] << endl;
  	 
  	 //pyramid 3
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*3;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 4
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*4;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 5
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*5;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j{0U}; j < iDim_j-1; j++) //y
  	for(auto i{0U}; i < iDim_i-1; i++) //x
  	{
  	 //pyramid 0
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1);
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][4]= (k==0?BACK_OUTSIDE:iNeighbor);
     
     //pyramid 1
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i +iNrOfCells*3;
     deqElementNeighbors[iElement][4]= (j==0?BOTTOM_OUTSIDE:iNeighbor);
     
      //pyramid 2
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1 +iNrOfCells*4;
     deqElementNeighbors[iElement][4]= (i==iDim_i-2?RIGHT_OUTSIDE:iNeighbor);
      
      //pyramid 3
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i +iNrOfCells*1;
     deqElementNeighbors[iElement][4]= (j==iDim_j-2?TOP_OUTSIDE:iNeighbor);
     
     //pyramid 4
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i-1 +iNrOfCells*2;
     deqElementNeighbors[iElement][4]= (i==0?LEFT_OUTSIDE:iNeighbor);
     
     //pyramid 5
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
    }
    
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j{0U}; j < iDim_j; j++) //y
  	for(auto i{0U}; i < iDim_i; i++) //x
  	{
  	  int8_t bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        // SKM FIX const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i));
        vset.BFlag( iNode, bBoundary);
      }
  	}
 
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    vset.EstablishZeroBasedNumbering();
    
    cout <<"\n"<<"create_Pyramid_VSet: model 'Pyra_Hexa': "<< endl;
    vset.Out();
    
} // end create_Pyramid_VSet



/**
      From ANSYS meshed FracBox with boundary surfaces and several tetrahedra
      that span model corners so that their rectification can be explored.
      
      @author SKM
      @date 2/2/22
*/
void create_FracBox( ModelTopology& topology, VSet<3U>& vset )
 {
    set<string> femTypes_matrix{"ISOSPARAMETRIC_LINEAR_TETRAHEDRON"},
                femTypes_surfaces{"ISOSPARAMETRIC_LINEAR_TRIANGLE"},
                femTypes_lines{"ISOSPARAMETRIC_LINEAR_BAR"};

    // MODEL TOPOLOGY
    topology.ModelName("FracBox");
    
    vector<size_t> elmtIdx_matrix(1363), elmtIdx_lines(1871-1782); // elmtIdx_surfaces(1782-1362)
    
    vector<size_t> elmtIdx_fracture{ 1516, 1519, 1525, 1526, 1527, 1528, 1537, 1539, 1545, 1546,
                                     1547, 1549, 1550, 1551, 1552, 1562, 1563, 1565, 1566, 1567,
                                     1595, 1596, 1597, 1630, 1631, 1765, 1766, 1767, 1768, 1769,
                                     1770, 1771, 1772, 1773, 1774, 1775, 1778, 1782 };
    assert( elmtIdx_fracture.size() == 38 );
    
    vector<size_t> elmtIdx_boundary1{ 1363, 1364, 1367, 1368, 1372, 1373, 1404, 1406, 1407, 1409,
                                      1439, 1440, 1441, 1442, 1443, 1444, 1445, 1446, 1447, 1448,
                                      1449, 1451, 1452, 1495, 1498, 1505, 1506, 1507, 1508, 1515,
                                      1517, 1518, 1520, 1521, 1523, 1524, 1568, 1581, 1583, 1587,
                                      1632, 1636, 1643, 1655, 1710, 1728, 1729, 1730, 1731, 1762,
                                      1779, 1780, 1781 };
    assert( elmtIdx_boundary1.size() == 53 );

    vector<size_t> elmtIdx_boundary2{ 1371, 1420, 1421, 1424, 1453, 1454, 1455, 1456, 1457, 1458,
                                      1459, 1460, 1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468,
                                      1469, 1474, 1475, 1478, 1480, 1486, 1487, 1488, 1489, 1535,
                                      1536, 1538, 1540, 1541, 1544, 1561, 1564, 1572, 1574, 1598,
                                      1614, 1617, 1618, 1619, 1688, 1726, 1732, 1733, 1734, 1735,
                                      1747, 1776, 1777 };
    assert( elmtIdx_boundary2.size() == 53 );

    vector<size_t> elmtIdx_boundary3{ 1365, 1366, 1390, 1392, 1394, 1395, 1397, 1398, 1399, 1402,
                                      1403, 1405, 1408, 1410, 1422, 1423, 1425, 1476, 1477, 1479,
                                      1496, 1497, 1499, 1522, 1530, 1531, 1532, 1533, 1534, 1542,
                                      1543, 1553, 1554, 1555, 1556, 1557, 1558, 1569, 1570, 1571,
                                      1616, 1621, 1622, 1623, 1624, 1625, 1629, 1644, 1651, 1652,
                                      1653, 1669, 1672, 1673, 1736, 1737, 1738, 1740, 1750, 1751,
                                      1752, 1753 };
    assert( elmtIdx_boundary3.size() == 62 );

    vector<size_t> elmtIdx_boundary4{ 1374, 1375, 1379, 1385, 1432, 1434, 1435, 1436, 1450, 1576,
                                      1577, 1578, 1579, 1580, 1585, 1586, 1588, 1589, 1590, 1591,
                                      1592, 1593, 1594, 1599, 1600, 1601, 1602, 1603, 1604, 1605,
                                      1606, 1607, 1608, 1609, 1610, 1611, 1612, 1658, 1659, 1660,
                                      1661, 1662, 1663, 1664, 1665, 1666, 1667, 1668, 1676, 1678,
                                      1679, 1681, 1682, 1683, 1684, 1685, 1686, 1687, 1690, 1696,
                                      1697, 1698, 1700, 1701, 1703, 1704, 1705, 1706, 1708, 1709,
                                      1711, 1712, 1714, 1715, 1719, 1720, 1721, 1722, 1723, 1724,
                                      1725 };
    assert( elmtIdx_boundary4.size() == 81 );

    vector<size_t> elmtIdx_boundary5{ 1376, 1377, 1378, 1380, 1381, 1382, 1383, 1396, 1400, 1401,
                                      1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1426,
                                      1427, 1428, 1429, 1430, 1431, 1509, 1510, 1511, 1512, 1513,
                                      1514, 1613, 1615, 1620, 1626, 1627, 1628, 1633, 1634, 1635,
                                      1637, 1638, 1639, 1640, 1641, 1642, 1645, 1646, 1647, 1648,
                                      1649, 1650, 1654, 1656, 1657, 1670, 1671, 1674, 1675, 1677,
                                      1680, 1689 };
    assert( elmtIdx_boundary5.size() == 62 );

    vector<size_t> elmtIdx_boundary6{ 1369, 1370, 1384, 1386, 1387, 1388, 1389, 1391, 1393, 1433,
                                      1437, 1438, 1470, 1471, 1472, 1473, 1481, 1482, 1483, 1484,
                                      1485, 1490, 1491, 1492, 1493, 1494, 1500, 1501, 1502, 1503,
                                      1504, 1529, 1548, 1559, 1560, 1573, 1575, 1582, 1584, 1691,
                                      1692, 1693, 1694, 1695, 1699, 1702, 1707, 1713, 1716, 1717,
                                      1718, 1727, 1739, 1741, 1742, 1743, 1744, 1745, 1746, 1748,
                                      1749, 1754, 1755, 1756, 1757, 1758, 1759, 1760, 1761, 1763,
                                      1764 };
    assert( elmtIdx_boundary6.size() == 71 );

    iota( elmtIdx_matrix.begin(), elmtIdx_matrix.end(), 0 );
//    iota( elmtIdx_surfaces.begin(), elmtIdx_surfaces.end(), 1363 );
    iota( elmtIdx_lines.begin(), elmtIdx_lines.end(), 1783 );

    topology.AddDomain( "MATRIX", femTypes_matrix, elmtIdx_matrix );
//    topology.AddDomain( "SURFACES", femTypes_surfaces, elmtIdx_surfaces ); // don't include because it would duplicate elements
    topology.AddDomain( "LINES", femTypes_lines, elmtIdx_lines );
    
    // specific regions
    topology.AddDomain( "FRACTURE", femTypes_surfaces, elmtIdx_fracture );
    // boundaries
    topology.AddDomain( "BOUNDARY1", femTypes_surfaces, elmtIdx_boundary1 );
    topology.AddDomain( "BOUNDARY2", femTypes_surfaces, elmtIdx_boundary2 );
    topology.AddDomain( "BOUNDARY3", femTypes_surfaces, elmtIdx_boundary3 );
    topology.AddDomain( "BOUNDARY4", femTypes_surfaces, elmtIdx_boundary4 );
    topology.AddDomain( "BOUNDARY5", femTypes_surfaces, elmtIdx_boundary5 );
    topology.AddDomain( "BOUNDARY6", femTypes_surfaces, elmtIdx_boundary6 );

   // VSET (types 25, 20, and 17)
   deque<int8_t> etypes{25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,
    25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17};

deque<uint32_t> npes{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

deque<uint32_t> epes{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

vset.Resize( etypes, npes, epes, 340, 0, 0 );

deque<double> px{10,-20,-14.7863,-16.4106,-17.9974,-20,-17.9978,-20,-18.1302,-20,-20,-18.6752,5.17197,4.51914,6.94351,7.99967,10,10,8.08487,10,8.12336,10,8.35642,7.99596,8.87834,-0.25,-9.75,-4.99898,-7.375,-4.96419,-3.6807,-2.625,-8.46843,-5,-5,-7.375,-20,-16.8652,-15.5225,-12.1065,-16.4565,-13.7137,-9.47408,-11.3335,10,6.875,6.17034,4.29862,6.17526,2.125,2.3088,-0.0711133,0.540269,10,-20,-20,10,6.05116,6.28089,-0.612867,1.55337,-20,-17.0708,-8.32689,-11.8521,-16.9777,-12.5418,-11.2241,10,6.27476,5.52513,1.12491,3.03348,6.875,2.12558,0.784153,-14.9974,-11.8192,-20,-17.0451,-12.2767,-16.4675,-10.4037,-8.44247,-8.96755,-20,-9.84178,-13.0033,-20,-16.5016,-14.8307,-16.5465,-16.7806,-20,-11.8309,-11.8677,-9.59169,-15.9227,-14.7642,-20,-20,-2.10792,-4.45709,-7.33115,-2.69461,-4.68386,-7.15291,-1.2522,-6.20454,-2.09285,-5,-0.411006,10,10,6.6344,10,10,-6.85097,-0.770139,-8.90791,-5,-10.0846,4.36597,1.69417,0.809702,-0.159916,-2.84229,1.77378,-0.25,2.16231,-2.625,-5,-3.03468,10,3.97979,7.02475,6.81378,10,-0.848291,0.218963,5.20927,10,-5.83225,-5,-1.09728,-7.57378,-14.3302,-12.1086,-14.5,-20,-16.5727,-18,-9.75,-15.7978,-11.3191,-5,-7.64289,-9.4228,-6.92254,-20,-20,4.5,7.69126,6.1603,8,4.45011,4.5,2.20374,-1.88823,-2.625,-5.43208,10,10,-2.14004,-2.625,-7.33022,-9.85179,-0.590538,-20,-20,-14.5043,-16.6879,-20,-20,-15.2367,-16.7205,-16.8994,-12.125,-12.5808,4.60757,6.66077,10,10,10,6.84099,6.413,0.703625,10,10,-0.828068,2.62263,0.371657,-0.25,-4.87766,-1.3456,-3.70971,4.80448,1.97791,2.73514,-0.464512,-1.45254,1.27204,-0.25,2.125,-1.65999,-5,-2.623,-2.625,10,10,10,4.9568,4.5,6.73536,6.875,10,2.08693,-0.158708,4.56346,5.68536,2.0892,10,4.56306,7.89394,7.45308,10,2.125,-0.451669,3.33755,7.00725,10,6.6955,10,-5,-11.2328,-9.75,-7.55794,-8.95647,-7.96507,-2.51988,-7.375,-4.95885,-14.5,-10.8266,-12.3076,-9.7589,-4.51998,-9.25436,-9.75,-12.125,-6.40167,-7.375,-20,-20,-20,-16.1714,-14.2059,-16.5024,-16.8301,-20,-12.125,-9.56144,-14.3925,-15.4132,-9.88683,-20,-13.9147,-17.8449,-16.875,-20,-12.2554,-9.52856,-16.0557,-20,-16.5515,-20,-20,-9.75,-7.2959,-8.12145,-5.01272,-4.63098,-4.69563,-0.19917,-9.71139,-9.75,-20,-14.6837,-14.3558,-20,-9.36466,-4.27531,5.43602,4.5,10,0.48896,-5.03467,10,4.48054,10,10,-14.4896,-20,-20,-5.16275,-4.24729,-0.25,-9.75,5.32655,-0.25,10,-3.99403,-13.717,-20,-12.8931,5.39208,-12.839,-5,-17.3625,3.4475,-13.9554,-11.7223,-5.05554,-2.55807,2.73119,-10.5618,10,-7.24118,-13.8747,-20};

deque<double> py{-1.87282,1.84099,-0.543371,-0.0434405,-0.0539297,2.95632,1.12129,0.667044,1.10917,-2.87765,-0.351769,-0.528406,-1.20451,0.427888,0.0528217,-0.00683389,-4.41135,-1.21829,-1.19107,-1.18562,-1.1888,1.53778,1.17556,1.37627,-0.00221258,0,0,-4.06863,-2.375,4.31846,1.85897,-2.375,2.47574,0,0,2.375,-10,-7.17215,-10,-7.24084,-10,-5.38478,-10,-6.7568,-10,-7.125,-10,-4.97826,-10,-7.125,-10,-10,-5.93378,1.26697,-0.640183,1.85132,-10,-6.83173,-10,-10,-6.45402,-10,-5.22909,-10,-6.38795,-10,-10,-10,10,6.31108,10,10,7.72788,7.125,7.11829,10,6.49753,7.12648,10,7.60762,7.06215,10,6.58985,10,10,10,10,10,10,10,3.92305,6.46722,7.41736,10,5.10262,7.38096,10,10,10,10,6.11666,10,10,7.65357,7.40272,5.06645,10,10,5.74477,7.49481,10,10,4.46009,10,6.52006,10,10,10,6.07033,6.16267,5,10,4.37136,1.99866,10,4.40185,7.6106,7.39584,0,2.67735,2.375,4.75,1.7575,10,4.33581,7.13863,7.20893,10,5.68543,10,10,5.86868,-2.27714,-6,-2.77168,1.87974,-1.31412,2.96744,0,-1.76874,2.75014,0,0,-3.9061,-1.72558,-4.75,-6.16717,-10,-7.41423,-6.4921,-10,0,-3.33952,3.07777,0,-2.04182,-4.75,-2.45106,-7.04514,-7.125,-10,-5.93606,-10,3.77235,2.375,2.3768,4.84044,5.70601,1.27067,2.23261,4.74636,2.33177,6.59565,5.46378,3.7507,2.59148,-2.35294,2.375,1.8009,4.13493,2.64828,5.31329,-0.315971,1.86796,2.44894,-2.60071,1.42289,-6.95138,-1.87256,-10,-10,-0.443097,-4.75,0.648293,-10,-4.87212,-4.28062,-1.58979,-10,-4.7508,-5.77195,-5.25329,0,-2.375,-2.48827,-4.75,-2.373,-7.125,-5.68223,-10,0,-3.39241,0,-2.6592,2.375,-1.00749,2.40322,0.141005,4.41814,2.07647,5.63151,-10,-4.06247,-7.63056,-6.6694,-10,-7.125,-10,-10,-10,1.83935,-2.61527,-5.1836,-10,-2.86978,-4.75,-2.24867,-10,-6.12456,-10,-7.125,-3.72035,-4.75,-1.4352,-10,-5.05327,-5.99003,-5.6688,0,-2.375,-1.8746,-2.375,-10,-5.60224,0,-4.58171,0.166213,-2.35376,2.58388,0.324203,2.375,-0.19321,4.9172,1.13185,1.99355,-10,-4.36432,-7.02899,-7.125,-10,-6.97182,-10,-10,-1.75295,-2.75211,-6.93868,10,4.75,6.79602,5.82313,4.32349,5.99233,10,5.11305,4.54,4.75,-5.71122,-10,0.213791,-5,-5.27233,-10,-10,0,-5.52499,-4.83094,-10,-3.64552,-0.0849293,4.04177,5.58605,-0.399898,5.9379,4.24812,0.315394,-10,-4.75,-4.75,10,4.75,5,10,10,5,10,10,-10,0,-10,-10,-10,10,10,10,10,-2.3379,-7.71165,-10,5.32574,-3.68022};

deque<double> pz{5.69268,13.547,10,11.9315,10,7.98355,8.51612,10.1406,11.0563,9.34521,6.61555,7.58385,10,10,12.3241,10,10.1984,8.86021,8.61557,11.8973,11.0234,9.85633,11.2472,8.7109,9.88055,7.5,7.5,7.44738,5.125,7.42513,10,5.125,10,12.25,2.75,5.125,16.4608,14.4801,15.562,14.5643,20,16.5907,16.564,20,15.8067,14.625,14.3347,16.9931,20,14.625,20,16.3883,20,6.33383,3.54018,4.96101,0,0,0,0,0,0,0,0,0,0,3.68278,0,0,0,0,0,0,14.625,14.9884,20,0,0,16.0208,13.9694,14.7653,20,20,12.2455,20,12.0589,2.69376,0,3.93654,0,10,5.02131,10.0638,7.91784,10,5.18663,7.28382,15.6801,4.64705,0,0,20,14.2041,14.3614,14.3905,16.9995,16.5375,16.9377,20,20,20,12.8347,15.5599,20,20,16.4096,12.7815,0,0,0,0,0,2.81911,5.01396,3.03495,6.27952,5.57185,5.33445,2.75,0,5.125,2.75,0,4.49355,10,9.71428,5.05267,8.94477,10,7.85568,11.7523,0,0,0,0,0,0,0,2.75,0,0,0,2.75,2.35902,4.57699,2.75,0,3.04152,6.20205,0,3.30938,2.75,0,0,0,0,2.75,5.20807,0,5.125,4.19262,0,3.56776,20,14.625,14.0102,16.9339,17.2775,2.22756,0,2.65025,5.38098,7.19323,14.7351,13.6729,10,4.52859,5.125,10,12.7574,10,7.55141,3.16138,0,5.16905,5.63988,10,4.2189,0,2.62276,0,20,17,20,20,20,12.5448,10,11.8182,7.486,10,10,12.25,14.625,10,12.25,13.4263,14.625,13.575,20,20,20,12.25,14.182,14.625,15.6137,13.7924,17.2303,17.099,20,20,11.622,10,9.85394,6.67557,7.14501,5.125,7.63796,3.9687,4.7083,12.7612,10,7.01469,20,20,17,20,20,20,12.3549,14.625,17.384,12.25,10,11.7859,7.69632,10,10,12.25,14.625,10,14.625,20,16.5879,20,20,12.4663,14.6741,14.7109,16.1341,14.625,17.0183,16.9421,20,20,12.8946,10,5.3807,9.875,8.2488,5.27616,7.35753,3.01606,12.6562,10,6.36547,20,7.5,5.42884,10,12.3153,10,4.29663,12.2479,12.6778,2.75,11.8517,7.70507,17.2267,20,12.5121,16.2721,8.88154,17,20,12.9717,8.77369,3.3478,7.47676,3.56032,11.1819,7.37507,11.7322,2.95062,17.274,0,2.75,2.75,6.00492,2.75,20,17.5359,11.1464,20,20,20,20,7.5,11.1902,17.0062,0,15.3037,9.20819,0,17.0473,0,10.0204,12.3255,20,5.17809};


vset.AddXYZ( px, py, pz );

//'plist' nodes that make up the elements
deque<vector<size_t> >  plist( 1872 );
	plist[0] = { 194, 193, 309, 192 };
	plist[1] = { 210, 216, 305, 214 };
	plist[2] = { 210, 217, 305, 216 };
	plist[3] = { 339, 153, 186, 277 };
	plist[4] = { 148, 178, 181, 150 };
	plist[5] = { 90, 2, 311, 188 };
	plist[6] = { 268, 312, 185, 8 };
	plist[7] = { 256, 158, 257, 27 };
	plist[8] = { 291, 103, 289, 290 };
	plist[9] = { 288, 332, 292, 96 };
	plist[10] = { 283, 8, 7, 1 };
	plist[11] = { 4, 9, 284, 283 };
	plist[12] = { 91, 100, 313, 88 };
	plist[13] = { 335, 154, 317, 152 };
	plist[14] = { 282, 62, 330, 65 };
	plist[15] = { 334, 230, 325, 75 };
	plist[16] = { 337, 249, 215, 250 };
	plist[17] = { 220, 303, 221, 229 };
	plist[18] = { 45, 223, 206, 47 };
	plist[19] = { 45, 223, 218, 206 };
	plist[20] = { 208, 232, 302, 211 };
	plist[21] = { 54, 151, 149, 178 };
	plist[22] = { 148, 54, 178, 151 };
	plist[23] = { 138, 290, 30, 293 };
	plist[24] = { 138, 29, 30, 291 };
	plist[25] = { 257, 27, 260, 256 };
	plist[26] = { 5, 8, 185, 312 };
	plist[27] = { 14, 190, 224, 22 };
	plist[28] = { 20, 22, 24, 15 };
	plist[29] = { 19, 22, 21, 24 };
	plist[30] = { 283, 284, 4, 3 };
	plist[31] = { 135, 139, 140, 318 };
	plist[32] = { 186, 339, 54, 149 };
	plist[33] = { 339, 62, 153, 277 };
	plist[34] = { 274, 298, 244, 273 };
	plist[35] = { 338, 82, 272, 324 };
	plist[36] = { 338, 82, 274, 272 };
	plist[37] = { 263, 40, 299, 262 };
	plist[38] = { 296, 267, 263, 37 };
	plist[39] = { 263, 275, 296, 37 };
	plist[40] = { 263, 275, 37, 36 };
	plist[41] = { 337, 39, 254, 42 };
	plist[42] = { 337, 281, 254, 257 };
	plist[43] = { 337, 257, 254, 300 };
	plist[44] = { 235, 234, 242, 197 };
	plist[45] = { 218, 233, 206, 16 };
	plist[46] = { 336, 231, 233, 235 };
	plist[47] = { 182, 91, 313, 88 };
	plist[48] = { 5, 181, 185, 6 };
	plist[49] = { 5, 55, 182, 181 };
	plist[50] = { 182, 55, 313, 181 };
	plist[51] = { 335, 146, 152, 147 };
	plist[52] = { 335, 64, 153, 146 };
	plist[53] = { 335, 153, 154, 146 };
	plist[54] = { 334, 107, 177, 75 };
	plist[55] = { 334, 73, 228, 74 };
	plist[56] = { 333, 131, 120, 118 };
	plist[57] = { 333, 126, 124, 292 };
	plist[58] = { 333, 71, 124, 118 };
	plist[59] = { 288, 29, 291, 289 };
	plist[60] = { 332, 111, 138, 139 };
	plist[61] = { 126, 288, 332, 292 };
	plist[62] = { 322, 80, 79, 184 };
	plist[63] = { 92, 184, 322, 79 };
	plist[64] = { 322, 80, 94, 83 };
	plist[65] = { 331, 106, 80, 83 };
	plist[66] = { 331, 82, 324, 80 };
	plist[67] = { 282, 159, 62, 65 };
	plist[68] = { 330, 67, 64, 66 };
	plist[69] = { 47, 329, 49, 45 };
	plist[70] = { 329, 52, 47, 49 };
	plist[71] = { 254, 276, 257, 297 };
	plist[72] = { 254, 278, 276, 297 };
	plist[73] = { 328, 37, 275, 38 };
	plist[74] = { 328, 296, 278, 279 };
	plist[75] = { 328, 279, 278, 297 };
	plist[76] = { 112, 115, 114, 320 };
	plist[77] = { 189, 224, 310, 73 };
	plist[78] = { 20, 19, 223, 16 };
	plist[79] = { 166, 162, 307, 57 };
	plist[80] = { 233, 218, 231, 336 };
	plist[81] = { 122, 136, 309, 69 };
	plist[82] = { 180, 150, 313, 76 };
	plist[83] = { 252, 267, 296, 37 };
	plist[84] = { 205, 52, 202, 201 };
	plist[85] = { 109, 230, 177, 75 };
	plist[86] = { 274, 108, 176, 82 };
	plist[87] = { 248, 244, 245, 43 };
	plist[88] = { 142, 152, 317, 28 };
	plist[89] = { 144, 155, 316, 31 };
	plist[90] = { 315, 170, 168, 143 };
	plist[91] = { 218, 233, 45, 206 };
	plist[92] = { 196, 212, 293, 30 };
	plist[93] = { 224, 222, 303, 226 };
	plist[94] = { 195, 12, 308, 232 };
	plist[95] = { 308, 190, 13, 15 };
	plist[96] = { 226, 222, 303, 213 };
	plist[97] = { 218, 336, 233, 16 };
	plist[98] = { 224, 112, 310, 73 };
	plist[99] = { 223, 225, 303, 221 };
	plist[100] = { 163, 122, 309, 69 };
	plist[101] = { 307, 192, 195, 0 };
	plist[102] = { 167, 128, 316, 31 };
	plist[103] = { 163, 161, 129, 165 };
	plist[104] = { 194, 161, 308, 123 };
	plist[105] = { 195, 161, 308, 194 };
	plist[106] = { 163, 161, 165, 164 };
	plist[107] = { 164, 161, 162, 192 };
	plist[108] = { 234, 166, 307, 57 };
	plist[109] = { 210, 211, 209, 214 };
	plist[110] = { 208, 233, 46, 302 };
	plist[111] = { 195, 232, 308, 167 };
	plist[112] = { 18, 12, 308, 195 };
	plist[113] = { 239, 197, 235, 172 };
	plist[114] = { 233, 235, 302, 234 };
	plist[115] = { 135, 189, 310, 73 };
	plist[116] = { 191, 194, 309, 53 };
	plist[117] = { 123, 125, 319, 130 };
	plist[118] = { 135, 127, 139, 318 };
	plist[119] = { 13, 308, 196, 134 };
	plist[120] = { 190, 134, 194, 308 };
	plist[121] = { 135, 136, 134, 127 };
	plist[122] = { 116, 135, 310, 73 };
	plist[123] = { 136, 133, 309, 69 };
	plist[124] = { 140, 137, 135, 318 };
	plist[125] = { 178, 181, 313, 55 };
	plist[126] = { 187, 152, 295, 35 };
	plist[127] = { 186, 148, 311, 154 };
	plist[128] = { 181, 148, 311, 186 };
	plist[129] = { 146, 148, 150, 151 };
	plist[130] = { 150, 313, 100, 179 };
	plist[131] = { 268, 184, 312, 79 };
	plist[132] = { 253, 258, 300, 260 };
	plist[133] = { 267, 266, 298, 259 };
	plist[134] = { 311, 284, 2, 4 };
	plist[135] = { 259, 266, 298, 270 };
	plist[136] = { 283, 284, 267, 296 };
	plist[137] = { 265, 40, 41, 326 };
	plist[138] = { 267, 269, 298, 268 };
	plist[139] = { 265, 269, 298, 267 };
	plist[140] = { 278, 252, 296, 37 };
	plist[141] = { 339, 284, 186, 11 };
	plist[142] = { 154, 255, 317, 28 };
	plist[143] = { 276, 297, 278, 277 };
	plist[144] = { 2, 252, 284, 276 };
	plist[145] = { 284, 276, 186, 311 };
	plist[146] = { 254, 37, 252, 278 };
	plist[147] = { 277, 279, 297, 278 };
	plist[148] = { 160, 277, 297, 282 };
	plist[149] = { 91, 180, 313, 76 };
	plist[150] = { 6, 4, 185, 8 };
	plist[151] = { 289, 94, 287, 32 };
	plist[152] = { 92, 184, 90, 322 };
	plist[153] = { 181, 90, 311, 187 };
	plist[154] = { 6, 185, 311, 181 };
	plist[155] = { 91, 89, 76, 100 };
	plist[156] = { 91, 100, 76, 313 };
	plist[157] = { 92, 85, 312, 79 };
	plist[158] = { 92, 93, 322, 98 };
	plist[159] = { 158, 155, 317, 28 };
	plist[160] = { 170, 143, 156, 155 };
	plist[161] = { 210, 249, 305, 217 };
	plist[162] = { 217, 249, 301, 250 };
	plist[163] = { 169, 209, 316, 31 };
	plist[164] = { 236, 237, 302, 211 };
	plist[165] = { 260, 257, 300, 253 };
	plist[166] = { 255, 281, 257, 297 };
	plist[167] = { 289, 290, 294, 32 };
	plist[168] = { 332, 138, 126, 139 };
	plist[169] = { 103, 102, 321, 104 };
	plist[170] = { 126, 131, 319, 130 };
	plist[171] = { 117, 86, 288, 292 };
	plist[172] = { 131, 120, 119, 117 };
	plist[173] = { 30, 138, 293, 196 };
	plist[174] = { 135, 127, 134, 139 };
	plist[175] = { 288, 287, 295, 35 };
	plist[176] = { 95, 96, 322, 94 };
	plist[177] = { 109, 173, 177, 230 };
	plist[178] = { 216, 212, 305, 214 };
	plist[179] = { 227, 213, 47, 202 };
	plist[180] = { 173, 227, 314, 174 };
	plist[181] = { 203, 227, 314, 173 };
	plist[182] = { 213, 227, 303, 226 };
	plist[183] = { 229, 303, 221, 201 };
	plist[184] = { 274, 246, 314, 203 };
	plist[185] = { 174, 290, 293, 30 };
	plist[186] = { 109, 105, 321, 104 };
	plist[187] = { 314, 176, 274, 105 };
	plist[188] = { 173, 105, 314, 203 };
	plist[189] = { 104, 105, 321, 103 };
	plist[190] = { 106, 110, 108, 321 };
	plist[191] = { 132, 128, 319, 130 };
	plist[192] = { 248, 246, 245, 244 };
	plist[193] = { 175, 258, 294, 32 };
	plist[194] = { 176, 298, 272, 270 };
	plist[195] = { 246, 271, 314, 261 };
	plist[196] = { 274, 246, 271, 314 };
	plist[197] = { 270, 271, 298, 259 };
	plist[198] = { 261, 215, 300, 260 };
	plist[199] = { 248, 251, 301, 250 };
	plist[200] = { 251, 227, 216, 314 };
	plist[201] = { 246, 251, 314, 203 };
	plist[202] = { 250, 251, 301, 217 };
	plist[203] = { 204, 301, 205, 51 };
	plist[204] = { 145, 131, 295, 35 };
	plist[205] = { 320, 115, 114, 113 };
	plist[206] = { 115, 113, 325, 114 };
	plist[207] = { 74, 226, 293, 189 };
	plist[208] = { 228, 224, 303, 226 };
	plist[209] = { 334, 177, 230, 75 };
	plist[210] = { 228, 189, 73, 224 };
	plist[211] = { 189, 228, 74, 226 };
	plist[212] = { 189, 228, 226, 224 };
	plist[213] = { 228, 189, 74, 73 };
	plist[214] = { 19, 17, 24, 21 };
	plist[215] = { 60, 165, 316, 166 };
	plist[216] = { 165, 57, 166, 60 };
	plist[217] = { 162, 57, 166, 165 };
	plist[218] = { 211, 206, 208, 232 };
	plist[219] = { 329, 45, 48, 46 };
	plist[220] = { 208, 45, 49, 46 };
	plist[221] = { 208, 206, 233, 232 };
	plist[222] = { 206, 208, 49, 305 };
	plist[223] = { 208, 233, 206, 45 };
	plist[224] = { 72, 127, 319, 122 };
	plist[225] = { 72, 318, 70, 124 };
	plist[226] = { 72, 136, 122, 69 };
	plist[227] = { 72, 136, 69, 70 };
	plist[228] = { 72, 127, 122, 136 };
	plist[229] = { 72, 136, 70, 318 };
	plist[230] = { 77, 147, 295, 180 };
	plist[231] = { 147, 76, 180, 77 };
	plist[232] = { 150, 76, 180, 147 };
	plist[233] = { 39, 259, 300, 252 };
	plist[234] = { 41, 267, 298, 259 };
	plist[235] = { 41, 252, 37, 267 };
	plist[236] = { 252, 41, 39, 259 };
	plist[237] = { 252, 41, 259, 267 };
	plist[238] = { 41, 252, 39, 37 };
	plist[239] = { 64, 280, 317, 153 };
	plist[240] = { 282, 153, 277, 280 };
	plist[241] = { 282, 64, 153, 280 };
	plist[242] = { 282, 64, 62, 153 };
	plist[243] = { 282, 64, 280, 66 };
	plist[244] = { 282, 280, 297, 66 };
	plist[245] = { 282, 330, 64, 66 };
	plist[246] = { 80, 94, 294, 184 };
	plist[247] = { 92, 93, 98, 91 };
	plist[248] = { 323, 97, 81, 272 };
	plist[249] = { 90, 185, 184, 312 };
	plist[250] = { 322, 80, 184, 94 };
	plist[251] = { 322, 95, 90, 92 };
	plist[252] = { 92, 312, 184, 79 };
	plist[253] = { 64, 156, 317, 157 };
	plist[254] = { 157, 63, 64, 156 };
	plist[255] = { 170, 63, 156, 315 };
	plist[256] = { 63, 157, 64, 67 };
	plist[257] = { 249, 51, 49, 217 };
	plist[258] = { 249, 305, 208, 210 };
	plist[259] = { 249, 51, 217, 301 };
	plist[260] = { 60, 236, 316, 199 };
	plist[261] = { 302, 237, 236, 238 };
	plist[262] = { 57, 200, 238, 60 };
	plist[263] = { 199, 238, 60, 236 };
	plist[264] = { 238, 237, 236, 199 };
	plist[265] = { 238, 199, 60, 200 };
	plist[266] = { 326, 38, 39, 42 };
	plist[267] = { 254, 257, 281, 297 };
	plist[268] = { 80, 103, 294, 83 };
	plist[269] = { 321, 102, 103, 106 };
	plist[270] = { 83, 106, 80, 103 };
	plist[271] = { 106, 102, 103, 83 };
	plist[272] = { 322, 331, 97, 80 };
	plist[273] = { 72, 118, 319, 124 };
	plist[274] = { 124, 70, 72, 71 };
	plist[275] = { 124, 71, 72, 118 };
	plist[276] = { 131, 288, 292, 117 };
	plist[277] = { 140, 111, 138, 293 };
	plist[278] = { 334, 74, 111, 140 };
	plist[279] = { 140, 293, 189, 74 };
	plist[280] = { 140, 139, 138, 111 };
	plist[281] = { 334, 74, 140, 73 };
	plist[282] = { 77, 95, 295, 86 };
	plist[283] = { 322, 96, 95, 98 };
	plist[284] = { 76, 87, 98, 77 };
	plist[285] = { 86, 98, 77, 95 };
	plist[286] = { 98, 96, 95, 86 };
	plist[287] = { 98, 86, 77, 87 };
	plist[288] = { 49, 213, 305, 202 };
	plist[289] = { 227, 47, 201, 202 };
	plist[290] = { 47, 202, 52, 201 };
	plist[291] = { 202, 47, 49, 213 };
	plist[292] = { 227, 47, 303, 201 };
	plist[293] = { 47, 202, 49, 52 };
	plist[294] = { 74, 104, 293, 177 };
	plist[295] = { 107, 109, 321, 104 };
	plist[296] = { 107, 177, 75, 109 };
	plist[297] = { 177, 107, 74, 104 };
	plist[298] = { 177, 107, 104, 109 };
	plist[299] = { 132, 131, 319, 118 };
	plist[300] = { 80, 270, 294, 176 };
	plist[301] = { 176, 271, 274, 298 };
	plist[302] = { 272, 176, 82, 274 };
	plist[303] = { 176, 272, 80, 270 };
	plist[304] = { 176, 298, 274, 272 };
	plist[305] = { 272, 176, 80, 82 };
	plist[306] = { 39, 250, 300, 245 };
	plist[307] = { 42, 248, 301, 250 };
	plist[308] = { 42, 245, 43, 248 };
	plist[309] = { 245, 42, 39, 250 };
	plist[310] = { 245, 42, 250, 248 };
	plist[311] = { 42, 245, 39, 43 };
	plist[312] = { 47, 48, 50, 329 };
	plist[313] = { 263, 40, 262, 36 };
	plist[314] = { 145, 35, 295, 152 };
	plist[315] = { 28, 35, 327, 34 };
	plist[316] = { 152, 34, 28, 142 };
	plist[317] = { 34, 152, 35, 145 };
	plist[318] = { 34, 152, 145, 142 };
	plist[319] = { 152, 34, 35, 28 };
	plist[320] = { 119, 147, 295, 77 };
	plist[321] = { 335, 142, 152, 317 };
	plist[322] = { 147, 145, 295, 152 };
	plist[323] = { 335, 152, 145, 147 };
	plist[324] = { 335, 152, 142, 145 };
	plist[325] = { 335, 317, 156, 142 };
	plist[326] = { 142, 28, 317, 155 };
	plist[327] = { 31, 28, 327, 34 };
	plist[328] = { 155, 34, 31, 144 };
	plist[329] = { 34, 155, 28, 142 };
	plist[330] = { 34, 155, 142, 144 };
	plist[331] = { 155, 34, 28, 31 };
	plist[332] = { 335, 317, 64, 156 };
	plist[333] = { 168, 155, 316, 144 };
	plist[334] = { 156, 142, 317, 155 };
	plist[335] = { 143, 170, 168, 155 };
	plist[336] = { 155, 143, 142, 144 };
	plist[337] = { 155, 143, 144, 168 };
	plist[338] = { 143, 155, 142, 156 };
	plist[339] = { 213, 49, 305, 206 };
	plist[340] = { 223, 47, 303, 213 };
	plist[341] = { 47, 48, 329, 45 };
	plist[342] = { 304, 48, 221, 47 };
	plist[343] = { 47, 206, 49, 213 };
	plist[344] = { 47, 206, 213, 223 };
	plist[345] = { 206, 47, 49, 45 };
	plist[346] = { 308, 23, 15, 18 };
	plist[347] = { 24, 53, 23, 21 };
	plist[348] = { 20, 241, 18, 16 };
	plist[349] = { 24, 18, 53, 17 };
	plist[350] = { 207, 214, 305, 212 };
	plist[351] = { 196, 25, 308, 207 };
	plist[352] = { 30, 214, 327, 25 };
	plist[353] = { 196, 214, 212, 30 };
	plist[354] = { 214, 196, 25, 30 };
	plist[355] = { 207, 214, 212, 196 };
	plist[356] = { 214, 207, 25, 196 };
	plist[357] = { 134, 13, 190, 189 };
	plist[358] = { 224, 189, 310, 190 };
	plist[359] = { 196, 293, 189, 134 };
	plist[360] = { 222, 189, 226, 224 };
	plist[361] = { 196, 189, 226, 222 };
	plist[362] = { 14, 190, 15, 13 };
	plist[363] = { 196, 189, 13, 134 };
	plist[364] = { 223, 213, 303, 222 };
	plist[365] = { 241, 206, 16, 223 };
	plist[366] = { 207, 213, 305, 206 };
	plist[367] = { 206, 207, 222, 12 };
	plist[368] = { 206, 222, 213, 223 };
	plist[369] = { 223, 206, 222, 12 };
	plist[370] = { 222, 206, 213, 207 };
	plist[371] = { 14, 225, 224, 303 };
	plist[372] = { 310, 224, 190, 22 };
	plist[373] = { 223, 16, 241, 20 };
	plist[374] = { 13, 14, 222, 12 };
	plist[375] = { 14, 225, 240, 224 };
	plist[376] = { 240, 14, 224, 22 };
	plist[377] = { 222, 14, 223, 12 };
	plist[378] = { 13, 207, 308, 12 };
	plist[379] = { 226, 212, 293, 196 };
	plist[380] = { 213, 207, 305, 212 };
	plist[381] = { 222, 212, 213, 226 };
	plist[382] = { 212, 222, 207, 196 };
	plist[383] = { 212, 222, 196, 226 };
	plist[384] = { 222, 212, 207, 213 };
	plist[385] = { 44, 304, 48, 219 };
	plist[386] = { 45, 44, 218, 304 };
	plist[387] = { 304, 47, 223, 45 };
	plist[388] = { 325, 73, 115, 114 };
	plist[389] = { 116, 112, 73, 310 };
	plist[390] = { 112, 115, 73, 114 };
	plist[391] = { 304, 45, 223, 218 };
	plist[392] = { 47, 223, 303, 221 };
	plist[393] = { 329, 52, 50, 47 };
	plist[394] = { 304, 223, 225, 218 };
	plist[395] = { 304, 47, 221, 223 };
	plist[396] = { 229, 114, 320, 112 };
	plist[397] = { 224, 228, 303, 229 };
	plist[398] = { 112, 228, 73, 224 };
	plist[399] = { 228, 112, 114, 229 };
	plist[400] = { 228, 112, 229, 224 };
	plist[401] = { 112, 228, 114, 73 };
	plist[402] = { 14, 303, 222, 223 };
	plist[403] = { 22, 240, 21, 310 };
	plist[404] = { 14, 19, 240, 225 };
	plist[405] = { 14, 225, 303, 223 };
	plist[406] = { 304, 223, 221, 225 };
	plist[407] = { 220, 224, 225, 303 };
	plist[408] = { 228, 229, 230, 201 };
	plist[409] = { 304, 225, 221, 220 };
	plist[410] = { 223, 218, 16, 19 };
	plist[411] = { 304, 48, 47, 45 };
	plist[412] = { 218, 225, 19, 223 };
	plist[413] = { 220, 112, 224, 229 };
	plist[414] = { 220, 112, 229, 320 };
	plist[415] = { 220, 303, 225, 221 };
	plist[416] = { 129, 72, 319, 122 };
	plist[417] = { 129, 69, 122, 163 };
	plist[418] = { 72, 69, 122, 129 };
	plist[419] = { 53, 194, 309, 192 };
	plist[420] = { 194, 18, 53, 23 };
	plist[421] = { 194, 18, 195, 0 };
	plist[422] = { 194, 0, 195, 192 };
	plist[423] = { 123, 130, 319, 128 };
	plist[424] = { 167, 25, 308, 123 };
	plist[425] = { 31, 130, 327, 25 };
	plist[426] = { 128, 25, 31, 167 };
	plist[427] = { 25, 128, 130, 123 };
	plist[428] = { 25, 128, 123, 167 };
	plist[429] = { 128, 25, 130, 31 };
	plist[430] = { 195, 167, 308, 161 };
	plist[431] = { 162, 166, 307, 195 };
	plist[432] = { 165, 167, 316, 166 };
	plist[433] = { 161, 166, 165, 162 };
	plist[434] = { 166, 161, 167, 195 };
	plist[435] = { 166, 161, 195, 162 };
	plist[436] = { 161, 166, 167, 165 };
	plist[437] = { 194, 122, 309, 163 };
	plist[438] = { 123, 129, 319, 122 };
	plist[439] = { 161, 122, 123, 194 };
	plist[440] = { 122, 161, 129, 163 };
	plist[441] = { 122, 161, 163, 194 };
	plist[442] = { 161, 122, 129, 123 };
	plist[443] = { 195, 192, 307, 162 };
	plist[444] = { 194, 193, 163, 309 };
	plist[445] = { 161, 192, 194, 195 };
	plist[446] = { 194, 164, 193, 192 };
	plist[447] = { 192, 161, 162, 195 };
	plist[448] = { 194, 161, 164, 192 };
	plist[449] = { 167, 123, 308, 161 };
	plist[450] = { 165, 128, 316, 167 };
	plist[451] = { 129, 123, 319, 128 };
	plist[452] = { 161, 128, 129, 165 };
	plist[453] = { 128, 161, 123, 167 };
	plist[454] = { 128, 161, 167, 165 };
	plist[455] = { 161, 128, 123, 129 };
	plist[456] = { 162, 171, 307, 57 };
	plist[457] = { 163, 309, 193, 141 };
	plist[458] = { 162, 192, 307, 198 };
	plist[459] = { 194, 164, 163, 193 };
	plist[460] = { 192, 164, 193, 198 };
	plist[461] = { 192, 164, 198, 162 };
	plist[462] = { 194, 161, 163, 164 };
	plist[463] = { 162, 198, 307, 171 };
	plist[464] = { 163, 309, 141, 69 };
	plist[465] = { 236, 60, 316, 166 };
	plist[466] = { 234, 238, 302, 236 };
	plist[467] = { 58, 200, 238, 57 };
	plist[468] = { 166, 238, 57, 234 };
	plist[469] = { 238, 166, 60, 236 };
	plist[470] = { 238, 166, 236, 234 };
	plist[471] = { 166, 238, 60, 57 };
	plist[472] = { 194, 18, 0, 53 };
	plist[473] = { 241, 16, 242, 18 };
	plist[474] = { 18, 0, 53, 17 };
	plist[475] = { 194, 0, 192, 53 };
	plist[476] = { 18, 242, 0, 17 };
	plist[477] = { 167, 31, 316, 209 };
	plist[478] = { 207, 25, 308, 167 };
	plist[479] = { 214, 31, 327, 25 };
	plist[480] = { 209, 25, 214, 207 };
	plist[481] = { 25, 209, 31, 167 };
	plist[482] = { 25, 209, 167, 207 };
	plist[483] = { 209, 25, 31, 214 };
	plist[484] = { 206, 12, 241, 232 };
	plist[485] = { 233, 206, 16, 241 };
	plist[486] = { 211, 207, 206, 232 };
	plist[487] = { 234, 302, 233, 232 };
	plist[488] = { 12, 232, 207, 308 };
	plist[489] = { 206, 232, 241, 233 };
	plist[490] = { 211, 206, 305, 208 };
	plist[491] = { 236, 232, 211, 302 };
	plist[492] = { 195, 166, 307, 234 };
	plist[493] = { 167, 236, 316, 166 };
	plist[494] = { 236, 195, 234, 166 };
	plist[495] = { 236, 195, 232, 234 };
	plist[496] = { 236, 167, 232, 195 };
	plist[497] = { 236, 167, 195, 166 };
	plist[498] = { 302, 234, 236, 232 };
	plist[499] = { 16, 233, 241, 242 };
	plist[500] = { 195, 234, 307, 242 };
	plist[501] = { 241, 234, 233, 232 };
	plist[502] = { 241, 242, 233, 234 };
	plist[503] = { 241, 234, 232, 195 };
	plist[504] = { 241, 242, 234, 195 };
	plist[505] = { 207, 167, 308, 232 };
	plist[506] = { 214, 211, 209, 207 };
	plist[507] = { 236, 167, 316, 209 };
	plist[508] = { 211, 167, 232, 236 };
	plist[509] = { 167, 211, 209, 236 };
	plist[510] = { 207, 167, 232, 211 };
	plist[511] = { 167, 207, 209, 211 };
	plist[512] = { 231, 218, 45, 44 };
	plist[513] = { 218, 206, 223, 16 };
	plist[514] = { 46, 44, 45, 48 };
	plist[515] = { 172, 56, 57, 171 };
	plist[516] = { 172, 56, 58, 57 };
	plist[517] = { 57, 307, 234, 197 };
	plist[518] = { 233, 218, 45, 231 };
	plist[519] = { 231, 233, 235, 302 };
	plist[520] = { 231, 46, 45, 233 };
	plist[521] = { 231, 233, 302, 46 };
	plist[522] = { 46, 231, 45, 44 };
	plist[523] = { 57, 197, 172, 171 };
	plist[524] = { 239, 238, 302, 234 };
	plist[525] = { 239, 57, 197, 172 };
	plist[526] = { 239, 238, 234, 57 };
	plist[527] = { 239, 57, 172, 58 };
	plist[528] = { 239, 238, 57, 58 };
	plist[529] = { 234, 197, 307, 242 };
	plist[530] = { 242, 235, 233, 234 };
	plist[531] = { 16, 233, 242, 336 };
	plist[532] = { 242, 235, 336, 233 };
	plist[533] = { 239, 234, 302, 235 };
	plist[534] = { 57, 307, 197, 171 };
	plist[535] = { 239, 197, 234, 235 };
	plist[536] = { 239, 57, 234, 197 };
	plist[537] = { 196, 293, 226, 189 };
	plist[538] = { 293, 134, 138, 140 };
	plist[539] = { 334, 73, 140, 115 };
	plist[540] = { 189, 140, 73, 135 };
	plist[541] = { 140, 111, 293, 74 };
	plist[542] = { 135, 139, 134, 140 };
	plist[543] = { 189, 140, 74, 73 };
	plist[544] = { 194, 308, 18, 23 };
	plist[545] = { 23, 53, 194, 191 };
	plist[546] = { 53, 23, 21, 191 };
	plist[547] = { 138, 126, 29, 291 };
	plist[548] = { 123, 25, 308, 196 };
	plist[549] = { 130, 30, 327, 25 };
	plist[550] = { 125, 25, 130, 123 };
	plist[551] = { 25, 125, 30, 196 };
	plist[552] = { 25, 125, 196, 123 };
	plist[553] = { 125, 25, 30, 130 };
	plist[554] = { 194, 123, 308, 134 };
	plist[555] = { 136, 122, 309, 194 };
	plist[556] = { 127, 123, 319, 122 };
	plist[557] = { 134, 122, 127, 136 };
	plist[558] = { 122, 134, 123, 194 };
	plist[559] = { 122, 134, 194, 136 };
	plist[560] = { 134, 122, 123, 127 };
	plist[561] = { 134, 140, 139, 138 };
	plist[562] = { 190, 189, 310, 135 };
	plist[563] = { 293, 134, 140, 189 };
	plist[564] = { 196, 189, 222, 13 };
	plist[565] = { 134, 189, 135, 140 };
	plist[566] = { 189, 134, 135, 190 };
	plist[567] = { 293, 196, 138, 134 };
	plist[568] = { 135, 136, 127, 318 };
	plist[569] = { 194, 191, 309, 136 };
	plist[570] = { 190, 135, 310, 191 };
	plist[571] = { 194, 134, 135, 136 };
	plist[572] = { 194, 135, 190, 191 };
	plist[573] = { 194, 135, 191, 136 };
	plist[574] = { 194, 134, 190, 135 };
	plist[575] = { 123, 196, 308, 134 };
	plist[576] = { 127, 125, 319, 123 };
	plist[577] = { 138, 30, 125, 196 };
	plist[578] = { 138, 126, 139, 127 };
	plist[579] = { 125, 134, 196, 123 };
	plist[580] = { 125, 134, 138, 196 };
	plist[581] = { 125, 127, 134, 123 };
	plist[582] = { 140, 115, 73, 116 };
	plist[583] = { 116, 115, 73, 112 };
	plist[584] = { 141, 69, 309, 133 };
	plist[585] = { 68, 133, 69, 141 };
	plist[586] = { 116, 137, 135, 140 };
	plist[587] = { 334, 73, 115, 325 };
	plist[588] = { 140, 116, 73, 135 };
	plist[589] = { 334, 228, 230, 74 };
	plist[590] = { 137, 133, 136, 318 };
	plist[591] = { 133, 70, 69, 136 };
	plist[592] = { 318, 133, 136, 70 };
	plist[593] = { 70, 133, 69, 68 };
	plist[594] = { 135, 136, 318, 137 };
	plist[595] = { 191, 137, 136, 135 };
	plist[596] = { 133, 191, 136, 309 };
	plist[597] = { 310, 135, 116, 137 };
	plist[598] = { 310, 135, 137, 191 };
	plist[599] = { 133, 137, 136, 191 };
	plist[600] = { 335, 154, 152, 146 };
	plist[601] = { 62, 149, 153, 151 };
	plist[602] = { 148, 54, 181, 178 };
	plist[603] = { 181, 10, 55, 54 };
	plist[604] = { 154, 28, 317, 152 };
	plist[605] = { 187, 26, 311, 154 };
	plist[606] = { 35, 28, 327, 26 };
	plist[607] = { 152, 26, 35, 187 };
	plist[608] = { 26, 152, 28, 154 };
	plist[609] = { 26, 152, 154, 187 };
	plist[610] = { 152, 26, 28, 35 };
	plist[611] = { 181, 187, 311, 148 };
	plist[612] = { 150, 180, 313, 181 };
	plist[613] = { 147, 187, 295, 180 };
	plist[614] = { 148, 180, 147, 150 };
	plist[615] = { 180, 148, 187, 181 };
	plist[616] = { 180, 148, 181, 150 };
	plist[617] = { 148, 180, 187, 147 };
	plist[618] = { 153, 151, 148, 146 };
	plist[619] = { 339, 153, 149, 186 };
	plist[620] = { 335, 64, 317, 153 };
	plist[621] = { 148, 153, 154, 186 };
	plist[622] = { 186, 151, 148, 153 };
	plist[623] = { 148, 153, 146, 154 };
	plist[624] = { 148, 186, 181, 54 };
	plist[625] = { 181, 178, 313, 150 };
	plist[626] = { 151, 186, 149, 153 };
	plist[627] = { 148, 186, 54, 151 };
	plist[628] = { 146, 148, 147, 150 };
	plist[629] = { 148, 178, 150, 151 };
	plist[630] = { 187, 154, 311, 148 };
	plist[631] = { 147, 152, 295, 187 };
	plist[632] = { 335, 153, 317, 154 };
	plist[633] = { 148, 152, 146, 147 };
	plist[634] = { 152, 148, 154, 187 };
	plist[635] = { 152, 148, 187, 147 };
	plist[636] = { 148, 152, 154, 146 };
	plist[637] = { 339, 159, 62, 277 };
	plist[638] = { 179, 150, 313, 178 };
	plist[639] = { 54, 55, 181, 178 };
	plist[640] = { 151, 178, 150, 179 };
	plist[641] = { 151, 179, 149, 178 };
	plist[642] = { 150, 313, 76, 100 };
	plist[643] = { 272, 80, 97, 324 };
	plist[644] = { 270, 80, 294, 184 };
	plist[645] = { 184, 266, 270, 298 };
	plist[646] = { 184, 272, 79, 268 };
	plist[647] = { 272, 184, 80, 270 };
	plist[648] = { 184, 298, 270, 272 };
	plist[649] = { 184, 272, 80, 79 };
	plist[650] = { 311, 11, 4, 6 };
	plist[651] = { 7, 8, 4, 6 };
	plist[652] = { 4, 11, 9, 7 };
	plist[653] = { 188, 32, 294, 258 };
	plist[654] = { 253, 26, 311, 188 };
	plist[655] = { 260, 32, 327, 26 };
	plist[656] = { 253, 32, 258, 260 };
	plist[657] = { 32, 253, 26, 260 };
	plist[658] = { 188, 32, 258, 253 };
	plist[659] = { 32, 188, 26, 253 };
	plist[660] = { 311, 284, 4, 11 };
	plist[661] = { 267, 252, 296, 284 };
	plist[662] = { 259, 253, 300, 252 };
	plist[663] = { 266, 252, 259, 267 };
	plist[664] = { 253, 252, 266, 2 };
	plist[665] = { 3, 284, 4, 2 };
	plist[666] = { 266, 252, 253, 259 };
	plist[667] = { 184, 298, 272, 268 };
	plist[668] = { 268, 185, 266, 3 };
	plist[669] = { 188, 270, 294, 184 };
	plist[670] = { 90, 184, 185, 266 };
	plist[671] = { 184, 266, 298, 268 };
	plist[672] = { 268, 185, 184, 266 };
	plist[673] = { 266, 184, 270, 188 };
	plist[674] = { 267, 268, 298, 266 };
	plist[675] = { 2, 3, 266, 185 };
	plist[676] = { 269, 3, 267, 283 };
	plist[677] = { 283, 3, 267, 284 };
	plist[678] = { 268, 312, 184, 185 };
	plist[679] = { 90, 2, 185, 311 };
	plist[680] = { 259, 258, 300, 253 };
	plist[681] = { 270, 188, 294, 258 };
	plist[682] = { 266, 258, 270, 259 };
	plist[683] = { 258, 266, 188, 253 };
	plist[684] = { 258, 266, 253, 259 };
	plist[685] = { 266, 258, 188, 270 };
	plist[686] = { 263, 36, 37, 40 };
	plist[687] = { 78, 323, 81, 286 };
	plist[688] = { 79, 78, 183, 323 };
	plist[689] = { 323, 78, 97, 79 };
	plist[690] = { 296, 267, 283, 263 };
	plist[691] = { 41, 265, 298, 267 };
	plist[692] = { 265, 40, 263, 37 };
	plist[693] = { 265, 37, 263, 267 };
	plist[694] = { 265, 37, 267, 41 };
	plist[695] = { 268, 79, 312, 183 };
	plist[696] = { 273, 272, 298, 268 };
	plist[697] = { 272, 324, 82, 80 };
	plist[698] = { 79, 323, 183, 272 };
	plist[699] = { 272, 183, 79, 268 };
	plist[700] = { 323, 268, 272, 273 };
	plist[701] = { 323, 97, 272, 79 };
	plist[702] = { 273, 265, 298, 244 };
	plist[703] = { 299, 265, 269, 264 };
	plist[704] = { 323, 268, 183, 272 };
	plist[705] = { 265, 273, 298, 269 };
	plist[706] = { 323, 269, 273, 264 };
	plist[707] = { 284, 283, 9, 296 };
	plist[708] = { 1, 8, 7, 312 };
	plist[709] = { 269, 268, 3, 1 };
	plist[710] = { 269, 268, 267, 3 };
	plist[711] = { 299, 265, 263, 269 };
	plist[712] = { 269, 263, 267, 265 };
	plist[713] = { 269, 263, 283, 267 };
	plist[714] = { 273, 268, 298, 269 };
	plist[715] = { 323, 269, 268, 273 };
	plist[716] = { 328, 254, 278, 37 };
	plist[717] = { 276, 300, 253, 252 };
	plist[718] = { 252, 254, 278, 276 };
	plist[719] = { 41, 39, 326, 38 };
	plist[720] = { 328, 254, 37, 38 };
	plist[721] = { 300, 276, 257, 254 };
	plist[722] = { 254, 300, 252, 39 };
	plist[723] = { 254, 39, 37, 38 };
	plist[724] = { 11, 311, 181, 6 };
	plist[725] = { 10, 181, 5, 11 };
	plist[726] = { 339, 9, 284, 11 };
	plist[727] = { 285, 284, 277, 339 };
	plist[728] = { 154, 26, 311, 253 };
	plist[729] = { 28, 260, 327, 26 };
	plist[730] = { 255, 26, 28, 154 };
	plist[731] = { 26, 255, 260, 253 };
	plist[732] = { 26, 255, 253, 154 };
	plist[733] = { 255, 26, 260, 28 };
	plist[734] = { 186, 154, 311, 276 };
	plist[735] = { 280, 154, 317, 153 };
	plist[736] = { 280, 186, 154, 153 };
	plist[737] = { 280, 186, 276, 154 };
	plist[738] = { 280, 277, 276, 186 };
	plist[739] = { 280, 277, 186, 153 };
	plist[740] = { 284, 252, 296, 278 };
	plist[741] = { 254, 300, 276, 252 };
	plist[742] = { 276, 2, 253, 311 };
	plist[743] = { 297, 276, 280, 277 };
	plist[744] = { 252, 276, 278, 284 };
	plist[745] = { 276, 300, 257, 253 };
	plist[746] = { 254, 39, 252, 37 };
	plist[747] = { 279, 285, 277, 160 };
	plist[748] = { 284, 285, 278, 296 };
	plist[749] = { 285, 278, 277, 284 };
	plist[750] = { 284, 276, 278, 277 };
	plist[751] = { 284, 276, 277, 186 };
	plist[752] = { 154, 253, 311, 276 };
	plist[753] = { 280, 255, 317, 154 };
	plist[754] = { 257, 260, 255, 253 };
	plist[755] = { 255, 281, 297, 280 };
	plist[756] = { 255, 276, 253, 154 };
	plist[757] = { 255, 276, 257, 253 };
	plist[758] = { 255, 280, 276, 154 };
	plist[759] = { 282, 159, 65, 160 };
	plist[760] = { 65, 160, 159, 61 };
	plist[761] = { 282, 277, 62, 159 };
	plist[762] = { 282, 280, 277, 297 };
	plist[763] = { 282, 277, 159, 160 };
	plist[764] = { 62, 153, 277, 282 };
	plist[765] = { 282, 64, 330, 62 };
	plist[766] = { 328, 296, 279, 275 };
	plist[767] = { 38, 40, 37, 36 };
	plist[768] = { 275, 37, 36, 38 };
	plist[769] = { 41, 38, 326, 40 };
	plist[770] = { 279, 285, 278, 277 };
	plist[771] = { 160, 285, 277, 159 };
	plist[772] = { 339, 62, 149, 153 };
	plist[773] = { 160, 277, 279, 297 };
	plist[774] = { 328, 37, 296, 275 };
	plist[775] = { 328, 37, 278, 296 };
	plist[776] = { 296, 278, 279, 285 };
	plist[777] = { 95, 77, 295, 180 };
	plist[778] = { 92, 98, 322, 95 };
	plist[779] = { 89, 87, 98, 76 };
	plist[780] = { 180, 98, 76, 91 };
	plist[781] = { 98, 180, 77, 95 };
	plist[782] = { 98, 180, 95, 91 };
	plist[783] = { 180, 98, 77, 76 };
	plist[784] = { 5, 7, 8, 312 };
	plist[785] = { 5, 181, 182, 185 };
	plist[786] = { 185, 5, 312, 182 };
	plist[787] = { 181, 5, 11, 6 };
	plist[788] = { 187, 35, 295, 287 };
	plist[789] = { 188, 26, 311, 187 };
	plist[790] = { 32, 35, 327, 26 };
	plist[791] = { 287, 26, 32, 188 };
	plist[792] = { 26, 287, 35, 187 };
	plist[793] = { 26, 287, 187, 188 };
	plist[794] = { 287, 26, 35, 32 };
	plist[795] = { 90, 184, 266, 188 };
	plist[796] = { 322, 90, 94, 184 };
	plist[797] = { 94, 188, 294, 184 };
	plist[798] = { 92, 312, 90, 184 };
	plist[799] = { 90, 266, 2, 188 };
	plist[800] = { 90, 312, 92, 182 };
	plist[801] = { 90, 184, 188, 94 };
	plist[802] = { 322, 95, 94, 90 };
	plist[803] = { 181, 180, 313, 91 };
	plist[804] = { 187, 95, 295, 180 };
	plist[805] = { 187, 91, 181, 90 };
	plist[806] = { 187, 180, 181, 91 };
	plist[807] = { 187, 91, 90, 95 };
	plist[808] = { 187, 180, 91, 95 };
	plist[809] = { 92, 91, 95, 90 };
	plist[810] = { 181, 182, 185, 90 };
	plist[811] = { 90, 182, 92, 91 };
	plist[812] = { 90, 312, 182, 185 };
	plist[813] = { 181, 182, 90, 91 };
	plist[814] = { 188, 187, 311, 90 };
	plist[815] = { 32, 94, 287, 188 };
	plist[816] = { 95, 187, 295, 287 };
	plist[817] = { 187, 94, 90, 188 };
	plist[818] = { 94, 187, 287, 188 };
	plist[819] = { 95, 187, 287, 94 };
	plist[820] = { 187, 95, 90, 94 };
	plist[821] = { 89, 88, 99, 100 };
	plist[822] = { 91, 98, 76, 89 };
	plist[823] = { 85, 183, 79, 78 };
	plist[824] = { 85, 312, 79, 183 };
	plist[825] = { 323, 78, 81, 97 };
	plist[826] = { 88, 93, 91, 98 };
	plist[827] = { 91, 98, 89, 88 };
	plist[828] = { 91, 89, 100, 88 };
	plist[829] = { 85, 92, 93, 322 };
	plist[830] = { 322, 80, 97, 79 };
	plist[831] = { 322, 79, 97, 85 };
	plist[832] = { 97, 85, 79, 78 };
	plist[833] = { 92, 98, 95, 91 };
	plist[834] = { 182, 93, 92, 91 };
	plist[835] = { 85, 182, 92, 312 };
	plist[836] = { 88, 182, 91, 93 };
	plist[837] = { 182, 181, 313, 91 };
	plist[838] = { 85, 93, 92, 182 };
	plist[839] = { 168, 60, 316, 199 };
	plist[840] = { 59, 199, 60, 168 };
	plist[841] = { 199, 59, 60, 200 };
	plist[842] = { 169, 31, 316, 155 };
	plist[843] = { 256, 169, 158, 27 };
	plist[844] = { 28, 31, 327, 27 };
	plist[845] = { 155, 27, 28, 158 };
	plist[846] = { 27, 155, 31, 169 };
	plist[847] = { 27, 155, 169, 158 };
	plist[848] = { 155, 27, 31, 28 };
	plist[849] = { 156, 158, 317, 157 };
	plist[850] = { 170, 157, 156, 63 };
	plist[851] = { 281, 170, 158, 306 };
	plist[852] = { 170, 157, 158, 156 };
	plist[853] = { 199, 315, 168, 59 };
	plist[854] = { 169, 168, 316, 199 };
	plist[855] = { 199, 170, 168, 315 };
	plist[856] = { 170, 199, 168, 169 };
	plist[857] = { 158, 169, 306, 170 };
	plist[858] = { 156, 155, 317, 158 };
	plist[859] = { 168, 169, 316, 155 };
	plist[860] = { 170, 315, 156, 143 };
	plist[861] = { 155, 170, 169, 158 };
	plist[862] = { 155, 170, 158, 156 };
	plist[863] = { 170, 155, 169, 168 };
	plist[864] = { 337, 39, 42, 250 };
	plist[865] = { 337, 250, 301, 249 };
	plist[866] = { 337, 300, 254, 39 };
	plist[867] = { 215, 257, 337, 300 };
	plist[868] = { 256, 158, 306, 257 };
	plist[869] = { 214, 260, 327, 27 };
	plist[870] = { 256, 214, 215, 210 };
	plist[871] = { 214, 256, 27, 210 };
	plist[872] = { 260, 214, 215, 256 };
	plist[873] = { 214, 260, 27, 256 };
	plist[874] = { 249, 49, 305, 217 };
	plist[875] = { 249, 49, 208, 305 };
	plist[876] = { 256, 249, 306, 237 };
	plist[877] = { 249, 51, 208, 49 };
	plist[878] = { 337, 250, 42, 301 };
	plist[879] = { 215, 257, 256, 337 };
	plist[880] = { 337, 215, 300, 250 };
	plist[881] = { 256, 237, 169, 210 };
	plist[882] = { 215, 260, 257, 300 };
	plist[883] = { 249, 215, 250, 217 };
	plist[884] = { 215, 249, 256, 210 };
	plist[885] = { 210, 249, 217, 215 };
	plist[886] = { 337, 249, 256, 215 };
	plist[887] = { 329, 49, 208, 51 };
	plist[888] = { 211, 207, 305, 206 };
	plist[889] = { 329, 49, 46, 208 };
	plist[890] = { 208, 233, 302, 232 };
	plist[891] = { 211, 214, 305, 207 };
	plist[892] = { 256, 210, 169, 27 };
	plist[893] = { 31, 214, 327, 27 };
	plist[894] = { 209, 27, 31, 169 };
	plist[895] = { 27, 209, 214, 210 };
	plist[896] = { 27, 209, 210, 169 };
	plist[897] = { 209, 27, 214, 31 };
	plist[898] = { 237, 170, 169, 199 };
	plist[899] = { 236, 169, 316, 199 };
	plist[900] = { 237, 306, 169, 170 };
	plist[901] = { 237, 199, 169, 236 };
	plist[902] = { 208, 45, 206, 49 };
	plist[903] = { 210, 211, 305, 208 };
	plist[904] = { 237, 208, 210, 249 };
	plist[905] = { 302, 211, 237, 208 };
	plist[906] = { 237, 208, 211, 210 };
	plist[907] = { 256, 237, 306, 169 };
	plist[908] = { 236, 209, 316, 169 };
	plist[909] = { 211, 210, 305, 214 };
	plist[910] = { 237, 209, 211, 236 };
	plist[911] = { 209, 237, 210, 169 };
	plist[912] = { 209, 237, 169, 236 };
	plist[913] = { 237, 209, 210, 211 };
	plist[914] = { 280, 64, 317, 157 };
	plist[915] = { 66, 157, 64, 280 };
	plist[916] = { 157, 66, 64, 67 };
	plist[917] = { 158, 28, 317, 255 };
	plist[918] = { 257, 337, 306, 256 };
	plist[919] = { 260, 28, 327, 27 };
	plist[920] = { 257, 306, 281, 158 };
	plist[921] = { 27, 255, 28, 158 };
	plist[922] = { 257, 27, 255, 260 };
	plist[923] = { 255, 27, 28, 260 };
	plist[924] = { 337, 256, 249, 306 };
	plist[925] = { 257, 337, 281, 306 };
	plist[926] = { 254, 328, 278, 297 };
	plist[927] = { 256, 249, 237, 210 };
	plist[928] = { 337, 300, 39, 250 };
	plist[929] = { 281, 66, 280, 157 };
	plist[930] = { 158, 280, 317, 157 };
	plist[931] = { 281, 157, 158, 170 };
	plist[932] = { 281, 297, 280, 66 };
	plist[933] = { 281, 157, 280, 158 };
	plist[934] = { 256, 169, 306, 158 };
	plist[935] = { 215, 260, 256, 257 };
	plist[936] = { 280, 158, 317, 255 };
	plist[937] = { 255, 297, 276, 280 };
	plist[938] = { 257, 158, 255, 27 };
	plist[939] = { 257, 281, 255, 158 };
	plist[940] = { 281, 255, 158, 280 };
	plist[941] = { 104, 74, 293, 111 };
	plist[942] = { 111, 102, 104, 107 };
	plist[943] = { 334, 107, 74, 177 };
	plist[944] = { 107, 111, 74, 104 };
	plist[945] = { 138, 104, 290, 293 };
	plist[946] = { 291, 288, 126, 29 };
	plist[947] = { 32, 30, 327, 29 };
	plist[948] = { 291, 32, 290, 289 };
	plist[949] = { 32, 291, 29, 289 };
	plist[950] = { 30, 32, 290, 291 };
	plist[951] = { 32, 30, 29, 291 };
	plist[952] = { 291, 332, 83, 102 };
	plist[953] = { 103, 289, 294, 83 };
	plist[954] = { 289, 83, 332, 96 };
	plist[955] = { 291, 83, 103, 102 };
	plist[956] = { 107, 102, 104, 321 };
	plist[957] = { 138, 290, 291, 30 };
	plist[958] = { 102, 111, 104, 291 };
	plist[959] = { 291, 83, 289, 103 };
	plist[960] = { 103, 290, 294, 289 };
	plist[961] = { 138, 111, 104, 293 };
	plist[962] = { 289, 83, 291, 332 };
	plist[963] = { 103, 291, 104, 290 };
	plist[964] = { 103, 291, 102, 104 };
	plist[965] = { 331, 106, 84, 82 };
	plist[966] = { 334, 74, 107, 111 };
	plist[967] = { 331, 80, 324, 97 };
	plist[968] = { 119, 77, 295, 86 };
	plist[969] = { 86, 117, 119, 121 };
	plist[970] = { 121, 86, 77, 119 };
	plist[971] = { 86, 121, 77, 87 };
	plist[972] = { 288, 35, 295, 131 };
	plist[973] = { 126, 288, 291, 332 };
	plist[974] = { 130, 35, 327, 29 };
	plist[975] = { 131, 29, 130, 126 };
	plist[976] = { 29, 131, 35, 288 };
	plist[977] = { 29, 131, 288, 126 };
	plist[978] = { 131, 29, 35, 130 };
	plist[979] = { 118, 126, 319, 124 };
	plist[980] = { 124, 292, 126, 139 };
	plist[981] = { 131, 292, 126, 333 };
	plist[982] = { 288, 119, 295, 86 };
	plist[983] = { 117, 86, 119, 288 };
	plist[984] = { 131, 288, 126, 292 };
	plist[985] = { 118, 131, 319, 126 };
	plist[986] = { 119, 288, 295, 131 };
	plist[987] = { 131, 292, 333, 117 };
	plist[988] = { 333, 131, 117, 120 };
	plist[989] = { 117, 131, 288, 119 };
	plist[990] = { 127, 72, 319, 124 };
	plist[991] = { 72, 318, 124, 127 };
	plist[992] = { 126, 130, 319, 125 };
	plist[993] = { 30, 130, 327, 29 };
	plist[994] = { 332, 138, 291, 126 };
	plist[995] = { 29, 125, 130, 126 };
	plist[996] = { 138, 29, 125, 30 };
	plist[997] = { 125, 29, 130, 30 };
	plist[998] = { 332, 291, 111, 102 };
	plist[999] = { 138, 104, 291, 290 };
	plist[1000] = { 332, 111, 291, 138 };
	plist[1001] = { 72, 136, 318, 127 };
	plist[1002] = { 126, 127, 319, 124 };
	plist[1003] = { 333, 118, 124, 126 };
	plist[1004] = { 124, 139, 127, 318 };
	plist[1005] = { 139, 124, 127, 126 };
	plist[1006] = { 126, 332, 139, 292 };
	plist[1007] = { 138, 111, 291, 104 };
	plist[1008] = { 127, 126, 319, 125 };
	plist[1009] = { 138, 127, 139, 134 };
	plist[1010] = { 138, 126, 125, 29 };
	plist[1011] = { 138, 126, 127, 125 };
	plist[1012] = { 138, 127, 134, 125 };
	plist[1013] = { 94, 80, 294, 83 };
	plist[1014] = { 322, 79, 85, 92 };
	plist[1015] = { 94, 32, 294, 188 };
	plist[1016] = { 289, 332, 291, 288 };
	plist[1017] = { 35, 32, 327, 29 };
	plist[1018] = { 287, 29, 35, 288 };
	plist[1019] = { 29, 287, 32, 289 };
	plist[1020] = { 29, 287, 289, 288 };
	plist[1021] = { 287, 29, 32, 35 };
	plist[1022] = { 292, 86, 288, 96 };
	plist[1023] = { 95, 288, 295, 86 };
	plist[1024] = { 96, 86, 288, 95 };
	plist[1025] = { 322, 331, 80, 83 };
	plist[1026] = { 289, 94, 294, 83 };
	plist[1027] = { 322, 94, 96, 83 };
	plist[1028] = { 96, 83, 94, 289 };
	plist[1029] = { 289, 332, 288, 96 };
	plist[1030] = { 95, 287, 295, 288 };
	plist[1031] = { 94, 289, 294, 32 };
	plist[1032] = { 96, 287, 94, 95 };
	plist[1033] = { 287, 96, 289, 288 };
	plist[1034] = { 287, 96, 288, 95 };
	plist[1035] = { 96, 287, 289, 94 };
	plist[1036] = { 226, 74, 293, 177 };
	plist[1037] = { 227, 177, 230, 228 };
	plist[1038] = { 228, 325, 114, 73 };
	plist[1039] = { 74, 230, 177, 228 };
	plist[1040] = { 228, 177, 74, 226 };
	plist[1041] = { 227, 228, 303, 226 };
	plist[1042] = { 174, 30, 293, 212 };
	plist[1043] = { 216, 33, 314, 174 };
	plist[1044] = { 260, 30, 327, 32 };
	plist[1045] = { 212, 33, 214, 216 };
	plist[1046] = { 33, 212, 30, 174 };
	plist[1047] = { 33, 212, 174, 216 };
	plist[1048] = { 212, 33, 30, 214 };
	plist[1049] = { 251, 202, 216, 227 };
	plist[1050] = { 227, 251, 205, 203 };
	plist[1051] = { 213, 216, 305, 202 };
	plist[1052] = { 227, 213, 303, 47 };
	plist[1053] = { 251, 227, 314, 203 };
	plist[1054] = { 227, 205, 201, 203 };
	plist[1055] = { 227, 202, 216, 213 };
	plist[1056] = { 228, 229, 201, 303 };
	plist[1057] = { 174, 226, 293, 177 };
	plist[1058] = { 227, 177, 174, 173 };
	plist[1059] = { 227, 177, 228, 226 };
	plist[1060] = { 177, 227, 230, 173 };
	plist[1061] = { 227, 177, 226, 174 };
	plist[1062] = { 228, 201, 230, 227 };
	plist[1063] = { 201, 173, 227, 203 };
	plist[1064] = { 230, 173, 227, 201 };
	plist[1065] = { 216, 174, 314, 227 };
	plist[1066] = { 213, 212, 305, 216 };
	plist[1067] = { 226, 174, 293, 212 };
	plist[1068] = { 227, 212, 226, 213 };
	plist[1069] = { 212, 227, 174, 216 };
	plist[1070] = { 212, 227, 216, 213 };
	plist[1071] = { 227, 212, 174, 226 };
	plist[1072] = { 228, 201, 227, 303 };
	plist[1073] = { 334, 228, 325, 230 };
	plist[1074] = { 325, 230, 228, 114 };
	plist[1075] = { 201, 47, 303, 221 };
	plist[1076] = { 221, 50, 47, 48 };
	plist[1077] = { 221, 52, 47, 50 };
	plist[1078] = { 220, 224, 303, 229 };
	plist[1079] = { 220, 112, 225, 224 };
	plist[1080] = { 334, 73, 325, 228 };
	plist[1081] = { 221, 201, 47, 52 };
	plist[1082] = { 45, 44, 304, 48 };
	plist[1083] = { 329, 51, 52, 49 };
	plist[1084] = { 103, 80, 294, 176 };
	plist[1085] = { 108, 106, 321, 103 };
	plist[1086] = { 176, 106, 82, 108 };
	plist[1087] = { 106, 176, 80, 103 };
	plist[1088] = { 106, 176, 103, 108 };
	plist[1089] = { 176, 106, 80, 82 };
	plist[1090] = { 175, 32, 294, 290 };
	plist[1091] = { 174, 33, 314, 175 };
	plist[1092] = { 214, 260, 33, 30 };
	plist[1093] = { 290, 33, 30, 174 };
	plist[1094] = { 33, 290, 32, 175 };
	plist[1095] = { 33, 290, 175, 174 };
	plist[1096] = { 290, 33, 32, 30 };
	plist[1097] = { 173, 174, 314, 105 };
	plist[1098] = { 104, 174, 293, 177 };
	plist[1099] = { 105, 177, 104, 109 };
	plist[1100] = { 177, 105, 174, 173 };
	plist[1101] = { 177, 105, 173, 109 };
	plist[1102] = { 105, 177, 174, 104 };
	plist[1103] = { 108, 103, 321, 105 };
	plist[1104] = { 175, 103, 294, 176 };
	plist[1105] = { 314, 176, 271, 274 };
	plist[1106] = { 176, 105, 103, 108 };
	plist[1107] = { 105, 274, 203, 108 };
	plist[1108] = { 105, 176, 103, 175 };
	plist[1109] = { 109, 321, 105, 108 };
	plist[1110] = { 108, 173, 105, 109 };
	plist[1111] = { 174, 175, 314, 105 };
	plist[1112] = { 104, 290, 293, 174 };
	plist[1113] = { 103, 175, 294, 290 };
	plist[1114] = { 105, 290, 103, 104 };
	plist[1115] = { 290, 105, 175, 174 };
	plist[1116] = { 290, 105, 174, 104 };
	plist[1117] = { 105, 290, 175, 103 };
	plist[1118] = { 84, 110, 108, 106 };
	plist[1119] = { 331, 82, 84, 324 };
	plist[1120] = { 106, 84, 82, 108 };
	plist[1121] = { 109, 107, 321, 101 };
	plist[1122] = { 109, 75, 107, 101 };
	plist[1123] = { 321, 109, 110, 108 };
	plist[1124] = { 108, 203, 105, 173 };
	plist[1125] = { 109, 101, 321, 110 };
	plist[1126] = { 144, 60, 316, 168 };
	plist[1127] = { 144, 31, 316, 128 };
	plist[1128] = { 145, 132, 34, 142 };
	plist[1129] = { 130, 31, 327, 34 };
	plist[1130] = { 128, 34, 130, 132 };
	plist[1131] = { 34, 128, 31, 144 };
	plist[1132] = { 34, 128, 144, 132 };
	plist[1133] = { 128, 34, 31, 130 };
	plist[1134] = { 144, 165, 316, 60 };
	plist[1135] = { 129, 128, 319, 132 };
	plist[1136] = { 165, 144, 316, 128 };
	plist[1137] = { 144, 129, 128, 165 };
	plist[1138] = { 144, 132, 128, 129 };
	plist[1139] = { 259, 39, 300, 245 };
	plist[1140] = { 244, 41, 298, 259 };
	plist[1141] = { 39, 326, 43, 41 };
	plist[1142] = { 245, 41, 43, 244 };
	plist[1143] = { 41, 245, 39, 259 };
	plist[1144] = { 41, 245, 259, 244 };
	plist[1145] = { 245, 41, 39, 43 };
	plist[1146] = { 261, 260, 300, 258 };
	plist[1147] = { 175, 33, 314, 261 };
	plist[1148] = { 30, 260, 33, 32 };
	plist[1149] = { 258, 33, 32, 175 };
	plist[1150] = { 33, 258, 260, 261 };
	plist[1151] = { 33, 258, 261, 175 };
	plist[1152] = { 258, 33, 260, 32 };
	plist[1153] = { 314, 175, 271, 176 };
	plist[1154] = { 270, 175, 294, 176 };
	plist[1155] = { 176, 271, 298, 270 };
	plist[1156] = { 314, 175, 176, 105 };
	plist[1157] = { 105, 274, 108, 176 };
	plist[1158] = { 271, 176, 175, 270 };
	plist[1159] = { 244, 259, 298, 271 };
	plist[1160] = { 261, 259, 300, 245 };
	plist[1161] = { 271, 245, 261, 246 };
	plist[1162] = { 245, 271, 259, 244 };
	plist[1163] = { 245, 271, 244, 246 };
	plist[1164] = { 271, 245, 259, 261 };
	plist[1165] = { 274, 244, 298, 271 };
	plist[1166] = { 274, 314, 105, 203 };
	plist[1167] = { 244, 246, 271, 274 };
	plist[1168] = { 175, 261, 314, 271 };
	plist[1169] = { 270, 258, 294, 175 };
	plist[1170] = { 259, 261, 300, 258 };
	plist[1171] = { 271, 258, 259, 270 };
	plist[1172] = { 258, 271, 261, 175 };
	plist[1173] = { 258, 271, 175, 270 };
	plist[1174] = { 271, 258, 261, 259 };
	plist[1175] = { 265, 43, 41, 244 };
	plist[1176] = { 338, 272, 323, 81 };
	plist[1177] = { 41, 244, 298, 265 };
	plist[1178] = { 41, 38, 40, 37 };
	plist[1179] = { 265, 326, 41, 43 };
	plist[1180] = { 263, 265, 299, 40 };
	plist[1181] = { 41, 39, 38, 37 };
	plist[1182] = { 273, 274, 298, 272 };
	plist[1183] = { 272, 97, 81, 324 };
	plist[1184] = { 338, 272, 81, 324 };
	plist[1185] = { 272, 80, 79, 97 };
	plist[1186] = { 273, 274, 272, 338 };
	plist[1187] = { 273, 265, 264, 269 };
	plist[1188] = { 323, 269, 183, 268 };
	plist[1189] = { 217, 49, 305, 202 };
	plist[1190] = { 205, 51, 301, 217 };
	plist[1191] = { 329, 51, 50, 52 };
	plist[1192] = { 202, 51, 52, 205 };
	plist[1193] = { 51, 202, 49, 217 };
	plist[1194] = { 51, 202, 217, 205 };
	plist[1195] = { 202, 51, 49, 52 };
	plist[1196] = { 210, 216, 214, 215 };
	plist[1197] = { 261, 33, 314, 216 };
	plist[1198] = { 260, 214, 327, 30 };
	plist[1199] = { 215, 33, 260, 261 };
	plist[1200] = { 33, 215, 214, 216 };
	plist[1201] = { 33, 215, 216, 261 };
	plist[1202] = { 215, 33, 214, 260 };
	plist[1203] = { 246, 261, 314, 251 };
	plist[1204] = { 250, 261, 300, 245 };
	plist[1205] = { 251, 245, 250, 248 };
	plist[1206] = { 245, 251, 261, 246 };
	plist[1207] = { 245, 251, 246, 248 };
	plist[1208] = { 251, 245, 261, 250 };
	plist[1209] = { 205, 217, 301, 251 };
	plist[1210] = { 216, 217, 305, 202 };
	plist[1211] = { 227, 251, 202, 205 };
	plist[1212] = { 202, 251, 217, 205 };
	plist[1213] = { 227, 205, 202, 201 };
	plist[1214] = { 251, 202, 217, 216 };
	plist[1215] = { 248, 301, 251, 205 };
	plist[1216] = { 248, 205, 251, 246 };
	plist[1217] = { 261, 216, 314, 251 };
	plist[1218] = { 250, 215, 300, 261 };
	plist[1219] = { 210, 217, 216, 215 };
	plist[1220] = { 251, 215, 217, 250 };
	plist[1221] = { 215, 251, 216, 261 };
	plist[1222] = { 215, 251, 261, 250 };
	plist[1223] = { 251, 215, 216, 217 };
	plist[1224] = { 51, 204, 52, 205 };
	plist[1225] = { 51, 204, 50, 52 };
	plist[1226] = { 247, 42, 43, 248 };
	plist[1227] = { 243, 42, 248, 301 };
	plist[1228] = { 247, 42, 326, 43 };
	plist[1229] = { 301, 248, 243, 205 };
	plist[1230] = { 204, 243, 205, 301 };
	plist[1231] = { 243, 247, 248, 42 };
	plist[1232] = { 118, 129, 319, 132 };
	plist[1233] = { 132, 130, 319, 131 };
	plist[1234] = { 35, 130, 327, 34 };
	plist[1235] = { 131, 34, 35, 145 };
	plist[1236] = { 34, 131, 130, 132 };
	plist[1237] = { 131, 132, 34, 145 };
	plist[1238] = { 131, 34, 130, 35 };
	plist[1239] = { 119, 131, 295, 145 };
	plist[1240] = { 119, 145, 295, 147 };
	plist[1241] = { 132, 131, 120, 145 };
	plist[1242] = { 131, 132, 120, 118 };
	plist[1243] = { 333, 126, 131, 118 };
	plist[1244] = { 131, 119, 120, 145 };
	plist[1245] = { 118, 72, 319, 129 };
	plist[1246] = { 19, 22, 240, 21 };
	plist[1247] = { 23, 21, 190, 22 };
	plist[1248] = { 21, 23, 24, 22 };
	plist[1249] = { 24, 18, 23, 53 };
	plist[1250] = { 24, 20, 18, 17 };
	plist[1251] = { 24, 18, 15, 23 };
	plist[1252] = { 20, 22, 19, 24 };
	plist[1253] = { 24, 20, 15, 18 };
	plist[1254] = { 22, 15, 190, 23 };
	plist[1255] = { 20, 17, 19, 16 };
	plist[1256] = { 19, 14, 20, 223 };
	plist[1257] = { 24, 53, 21, 17 };
	plist[1258] = { 15, 24, 23, 22 };
	plist[1259] = { 19, 22, 20, 14 };
	plist[1260] = { 14, 190, 22, 15 };
	plist[1261] = { 12, 20, 241, 15 };
	plist[1262] = { 14, 15, 12, 13 };
	plist[1263] = { 22, 14, 15, 20 };
	plist[1264] = { 14, 15, 20, 12 };
	plist[1265] = { 4, 6, 11, 7 };
	plist[1266] = { 4, 9, 283, 7 };
	plist[1267] = { 10, 9, 339, 11 };
	plist[1268] = { 7, 6, 11, 5 };
	plist[1269] = { 9, 284, 11, 4 };
	plist[1270] = { 3, 268, 8, 1 };
	plist[1271] = { 3, 4, 185, 2 };
	plist[1272] = { 283, 4, 8, 3 };
	plist[1273] = { 3, 4, 8, 185 };
	plist[1274] = { 20, 24, 19, 17 };
	plist[1275] = { 11, 5, 10, 7 };
	plist[1276] = { 11, 7, 10, 9 };
	plist[1277] = { 190, 21, 310, 22 };
	plist[1278] = { 196, 13, 222, 207 };
	plist[1279] = { 207, 13, 222, 12 };
	plist[1280] = { 12, 15, 241, 18 };
	plist[1281] = { 13, 12, 308, 15 };
	plist[1282] = { 253, 2, 266, 188 };
	plist[1283] = { 90, 266, 185, 2 };
	plist[1284] = { 2, 185, 311, 4 };
	plist[1285] = { 240, 22, 224, 310 };
	plist[1286] = { 190, 308, 194, 23 };
	plist[1287] = { 12, 18, 241, 195 };
	plist[1288] = { 194, 308, 195, 18 };
	plist[1289] = { 13, 196, 308, 207 };
	plist[1290] = { 308, 190, 15, 23 };
	plist[1291] = { 14, 222, 224, 189 };
	plist[1292] = { 12, 206, 207, 232 };
	plist[1293] = { 206, 223, 241, 12 };
	plist[1294] = { 14, 19, 225, 223 };
	plist[1295] = { 14, 303, 224, 222 };
	plist[1296] = { 14, 189, 190, 13 };
	plist[1297] = { 14, 189, 224, 190 };
	plist[1298] = { 14, 222, 189, 13 };
	plist[1299] = { 19, 22, 14, 240 };
	plist[1300] = { 241, 12, 223, 20 };
	plist[1301] = { 14, 20, 223, 12 };
	plist[1302] = { 15, 12, 308, 18 };
	plist[1303] = { 18, 242, 195, 0 };
	plist[1304] = { 16, 17, 242, 18 };
	plist[1305] = { 242, 195, 241, 18 };
	plist[1306] = { 191, 310, 190, 21 };
	plist[1307] = { 21, 190, 191, 23 };
	plist[1308] = { 194, 191, 190, 23 };
	plist[1309] = { 10, 181, 55, 5 };
	plist[1310] = { 5, 7, 6, 8 };
	plist[1311] = { 284, 311, 186, 11 };
	plist[1312] = { 11, 311, 186, 181 };
	plist[1313] = { 10, 186, 181, 11 };
	plist[1314] = { 2, 253, 311, 188 };
	plist[1315] = { 252, 2, 253, 276 };
	plist[1316] = { 3, 266, 267, 252 };
	plist[1317] = { 181, 185, 311, 90 };
	plist[1318] = { 267, 268, 266, 3 };
	plist[1319] = { 3, 252, 284, 2 };
	plist[1320] = { 3, 252, 267, 284 };
	plist[1321] = { 3, 266, 252, 2 };
	plist[1322] = { 283, 8, 1, 3 };
	plist[1323] = { 268, 185, 3, 8 };
	plist[1324] = { 283, 4, 7, 8 };
	plist[1325] = { 284, 285, 296, 9 };
	plist[1326] = { 285, 284, 339, 9 };
	plist[1327] = { 10, 11, 339, 186 };
	plist[1328] = { 339, 284, 277, 186 };
	plist[1329] = { 4, 185, 311, 6 };
	plist[1330] = { 39, 326, 42, 43 };
	plist[1331] = { 334, 177, 74, 230 };
	plist[1332] = { 54, 151, 186, 149 };
	plist[1333] = { 265, 40, 37, 41 };
	plist[1334] = { 255, 297, 257, 276 };
	plist[1335] = { 208, 233, 45, 46 };
	plist[1336] = { 20, 18, 17, 16 };
	plist[1337] = { 20, 241, 15, 18 };
	plist[1338] = { 5, 8, 6, 185 };
	plist[1339] = { 190, 134, 308, 13 };
	plist[1340] = { 12, 195, 241, 232 };
	plist[1341] = { 284, 276, 311, 2 };
	plist[1342] = { 329, 49, 45, 46 };
	plist[1343] = { 42, 38, 39, 254 };
	plist[1344] = { 165, 164, 161, 162 };
	plist[1345] = { 331, 106, 82, 80 };
	plist[1346] = { 142, 132, 34, 144 };
	plist[1347] = { 62, 151, 153, 146 };
	plist[1348] = { 62, 146, 153, 64 };
	plist[1349] = { 339, 159, 149, 62 };
	plist[1350] = { 339, 277, 285, 159 };
	plist[1351] = { 246, 205, 251, 203 };
	plist[1352] = { 338, 272, 273, 323 };
	plist[1353] = { 114, 230, 228, 229 };
	plist[1354] = { 181, 10, 54, 186 };
	plist[1355] = { 310, 112, 224, 240 };
	plist[1356] = { 112, 225, 224, 240 };
	plist[1357] = { 242, 307, 195, 0 };
	plist[1358] = { 183, 312, 268, 1 };
	plist[1359] = { 269, 183, 268, 1 };
	plist[1360] = { 8, 1, 268, 312 };
	plist[1361] = { 269, 3, 283, 1 };
	plist[1362] = { 186, 339, 10, 54 };
 // triangular elements
	plist[1363] = { 149, 54, 339 };
	plist[1364] = { 339, 10, 9 };
	plist[1365] = { 338, 324, 81 };
	plist[1366] = { 82, 338, 274 };
	plist[1367] = { 263, 296, 283 };
	plist[1368] = { 36, 275, 263 };
	plist[1369] = { 337, 281, 254 };
	plist[1370] = { 337, 42, 301 };
	plist[1371] = { 231, 336, 235 };
	plist[1372] = { 182, 313, 88 };
	plist[1373] = { 55, 182, 5 };
	plist[1374] = { 335, 147, 145 };
	plist[1375] = { 64, 146, 335 };
	plist[1376] = { 325, 75, 334 };
	plist[1377] = { 334, 107, 111 };
	plist[1378] = { 333, 124, 292 };
	plist[1379] = { 333, 118, 120 };
	plist[1380] = { 292, 332, 96 };
	plist[1381] = { 332, 139, 111 };
	plist[1382] = { 331, 84, 324 };
	plist[1383] = { 83, 106, 331 };
	plist[1384] = { 330, 66, 67 };
	plist[1385] = { 64, 67, 330 };
	plist[1386] = { 329, 50, 48 };
	plist[1387] = { 329, 51, 50 };
	plist[1388] = { 297, 279, 328 };
	plist[1389] = { 275, 38, 328 };
	plist[1390] = { 326, 265, 40 };
	plist[1391] = { 326, 38, 40 };
	plist[1392] = { 326, 43, 265 };
	plist[1393] = { 329, 48, 46 };
	plist[1394] = { 50, 221, 52 };
	plist[1395] = { 50, 48, 221 };
	plist[1396] = { 115, 325, 334 };
	plist[1397] = { 325, 75, 230 };
	plist[1398] = { 325, 230, 114 };
	plist[1399] = { 338, 81, 323 };
	plist[1400] = { 324, 81, 97 };
	plist[1401] = { 331, 324, 97 };
	plist[1402] = { 324, 82, 84 };
	plist[1403] = { 324, 338, 82 };
	plist[1404] = { 323, 183, 78 };
	plist[1405] = { 323, 81, 286 };
	plist[1406] = { 323, 78, 286 };
	plist[1407] = { 323, 269, 183 };
	plist[1408] = { 323, 264, 273 };
	plist[1409] = { 323, 264, 269 };
	plist[1410] = { 338, 323, 273 };
	plist[1411] = { 322, 93, 98 };
	plist[1412] = { 322, 98, 96 };
	plist[1413] = { 85, 322, 97 };
	plist[1414] = { 83, 322, 96 };
	plist[1415] = { 321, 101, 110 };
	plist[1416] = { 321, 110, 106 };
	plist[1417] = { 321, 106, 102 };
	plist[1418] = { 321, 102, 107 };
	plist[1419] = { 321, 107, 101 };
	plist[1420] = { 320, 115, 112 };
	plist[1421] = { 320, 113, 115 };
	plist[1422] = { 320, 113, 114 };
	plist[1423] = { 320, 114, 229 };
	plist[1424] = { 320, 112, 220 };
	plist[1425] = { 320, 229, 220 };
	plist[1426] = { 318, 137, 140 };
	plist[1427] = { 318, 140, 139 };
	plist[1428] = { 318, 133, 137 };
	plist[1429] = { 318, 70, 133 };
	plist[1430] = { 318, 139, 124 };
	plist[1431] = { 318, 124, 70 };
	plist[1432] = { 315, 143, 168 };
	plist[1433] = { 315, 63, 170 };
	plist[1434] = { 315, 168, 59 };
	plist[1435] = { 315, 63, 156 };
	plist[1436] = { 315, 156, 143 };
	plist[1437] = { 315, 170, 199 };
	plist[1438] = { 315, 199, 59 };
	plist[1439] = { 313, 55, 178 };
	plist[1440] = { 313, 179, 100 };
	plist[1441] = { 313, 100, 88 };
	plist[1442] = { 313, 178, 179 };
	plist[1443] = { 313, 182, 55 };
	plist[1444] = { 312, 85, 183 };
	plist[1445] = { 312, 5, 182 };
	plist[1446] = { 312, 7, 5 };
	plist[1447] = { 312, 1, 7 };
	plist[1448] = { 312, 183, 1 };
	plist[1449] = { 312, 182, 85 };
	plist[1450] = { 149, 62, 159 };
	plist[1451] = { 159, 149, 339 };
	plist[1452] = { 159, 339, 285 };
	plist[1453] = { 310, 112, 116 };
	plist[1454] = { 310, 21, 240 };
	plist[1455] = { 310, 240, 112 };
	plist[1456] = { 310, 191, 21 };
	plist[1457] = { 310, 137, 191 };
	plist[1458] = { 309, 193, 192 };
	plist[1459] = { 309, 53, 191 };
	plist[1460] = { 309, 192, 53 };
	plist[1461] = { 309, 141, 193 };
	plist[1462] = { 309, 133, 141 };
	plist[1463] = { 309, 191, 133 };
	plist[1464] = { 307, 197, 242 };
	plist[1465] = { 307, 0, 192 };
	plist[1466] = { 307, 192, 198 };
	plist[1467] = { 307, 198, 171 };
	plist[1468] = { 307, 242, 0 };
	plist[1469] = { 307, 171, 197 };
	plist[1470] = { 306, 237, 170 };
	plist[1471] = { 306, 337, 249 };
	plist[1472] = { 306, 249, 237 };
	plist[1473] = { 306, 170, 281 };
	plist[1474] = { 304, 218, 44 };
	plist[1475] = { 304, 44, 219 };
	plist[1476] = { 304, 48, 219 };
	plist[1477] = { 304, 221, 48 };
	plist[1478] = { 304, 225, 218 };
	plist[1479] = { 304, 220, 221 };
	plist[1480] = { 304, 220, 225 };
	plist[1481] = { 302, 239, 238 };
	plist[1482] = { 302, 238, 237 };
	plist[1483] = { 231, 302, 46 };
	plist[1484] = { 302, 235, 239 };
	plist[1485] = { 208, 302, 237 };
	plist[1486] = { 218, 336, 231 };
	plist[1487] = { 242, 336, 16 };
	plist[1488] = { 16, 218, 19 };
	plist[1489] = { 218, 16, 336 };
	plist[1490] = { 301, 51, 249 };
	plist[1491] = { 301, 204, 51 };
	plist[1492] = { 337, 301, 249 };
	plist[1493] = { 301, 243, 204 };
	plist[1494] = { 301, 42, 243 };
	plist[1495] = { 263, 299, 262 };
	plist[1496] = { 299, 262, 40 };
	plist[1497] = { 299, 40, 265 };
	plist[1498] = { 269, 299, 263 };
	plist[1499] = { 299, 265, 264 };
	plist[1500] = { 297, 160, 279 };
	plist[1501] = { 254, 328, 38 };
	plist[1502] = { 297, 66, 282 };
	plist[1503] = { 254, 297, 328 };
	plist[1504] = { 297, 281, 66 };
	plist[1505] = { 296, 9, 283 };
	plist[1506] = { 296, 263, 275 };
	plist[1507] = { 296, 285, 9 };
	plist[1508] = { 296, 279, 285 };
	plist[1509] = { 292, 86, 117 };
	plist[1510] = { 292, 124, 139 };
	plist[1511] = { 83, 96, 332 };
	plist[1512] = { 292, 139, 332 };
	plist[1513] = { 333, 292, 117 };
	plist[1514] = { 292, 96, 86 };
	plist[1515] = { 339, 9, 285 };
	plist[1516] = { 284, 4, 2 };
	plist[1517] = { 283, 9, 7 };
	plist[1518] = { 283, 7, 1 };
	plist[1519] = { 276, 284, 2 };
	plist[1520] = { 296, 275, 279 };
	plist[1521] = { 269, 283, 1 };
	plist[1522] = { 264, 265, 273 };
	plist[1523] = { 263, 283, 269 };
	plist[1524] = { 269, 264, 299 };
	plist[1525] = { 256, 257, 260 };
	plist[1526] = { 253, 260, 257 };
	plist[1527] = { 253, 276, 2 };
	plist[1528] = { 253, 257, 276 };
	plist[1529] = { 306, 281, 337 };
	plist[1530] = { 244, 248, 246 };
	plist[1531] = { 244, 273, 265 };
	plist[1532] = { 244, 246, 274 };
	plist[1533] = { 274, 338, 273 };
	plist[1534] = { 243, 248, 247 };
	plist[1535] = { 242, 17, 0 };
	plist[1536] = { 242, 16, 17 };
	plist[1537] = { 241, 12, 15 };
	plist[1538] = { 240, 21, 19 };
	plist[1539] = { 232, 12, 241 };
	plist[1540] = { 235, 242, 197 };
	plist[1541] = { 225, 240, 19 };
	plist[1542] = { 201, 229, 230 };
	plist[1543] = { 220, 229, 221 };
	plist[1544] = { 218, 225, 19 };
	plist[1545] = { 214, 256, 260 };
	plist[1546] = { 210, 214, 211 };
	plist[1547] = { 210, 256, 214 };
	plist[1548] = { 208, 237, 249 };
	plist[1549] = { 207, 13, 12 };
	plist[1550] = { 207, 12, 232 };
	plist[1551] = { 207, 232, 211 };
	plist[1552] = { 207, 211, 214 };
	plist[1553] = { 205, 246, 248 };
	plist[1554] = { 205, 248, 243 };
	plist[1555] = { 204, 205, 243 };
	plist[1556] = { 203, 246, 205 };
	plist[1557] = { 201, 221, 229 };
	plist[1558] = { 201, 203, 205 };
	plist[1559] = { 199, 237, 238 };
	plist[1560] = { 199, 238, 200 };
	plist[1561] = { 235, 336, 242 };
	plist[1562] = { 196, 207, 214 };
	plist[1563] = { 196, 13, 207 };
	plist[1564] = { 192, 193, 198 };
	plist[1565] = { 190, 15, 13 };
	plist[1566] = { 90, 188, 2 };
	plist[1567] = { 188, 253, 2 };
	plist[1568] = { 183, 269, 1 };
	plist[1569] = { 203, 274, 246 };
	plist[1570] = { 173, 203, 201 };
	plist[1571] = { 173, 201, 230 };
	plist[1572] = { 172, 235, 197 };
	plist[1573] = { 172, 239, 235 };
	plist[1574] = { 171, 172, 197 };
	plist[1575] = { 170, 237, 199 };
	plist[1576] = { 164, 193, 198 };
	plist[1577] = { 163, 193, 164 };
	plist[1578] = { 162, 165, 164 };
	plist[1579] = { 162, 198, 171 };
	plist[1580] = { 162, 164, 198 };
	plist[1581] = { 160, 285, 279 };
	plist[1582] = { 297, 282, 160 };
	plist[1583] = { 159, 285, 160 };
	plist[1584] = { 157, 281, 170 };
	plist[1585] = { 150, 151, 179 };
	plist[1586] = { 149, 179, 151 };
	plist[1587] = { 149, 179, 178 };
	plist[1588] = { 146, 150, 147 };
	plist[1589] = { 146, 151, 150 };
	plist[1590] = { 143, 144, 168 };
	plist[1591] = { 142, 144, 143 };
	plist[1592] = { 335, 145, 142 };
	plist[1593] = { 335, 146, 147 };
	plist[1594] = { 142, 143, 156 };
	plist[1595] = { 134, 13, 196 };
	plist[1596] = { 134, 190, 13 };
	plist[1597] = { 134, 196, 138 };
	plist[1598] = { 133, 191, 137 };
	plist[1599] = { 132, 144, 142 };
	plist[1600] = { 132, 142, 145 };
	plist[1601] = { 165, 163, 164 };
	plist[1602] = { 165, 129, 163 };
	plist[1603] = { 129, 165, 144 };
	plist[1604] = { 129, 144, 132 };
	plist[1605] = { 120, 132, 145 };
	plist[1606] = { 119, 145, 147 };
	plist[1607] = { 119, 120, 145 };
	plist[1608] = { 118, 129, 132 };
	plist[1609] = { 118, 132, 120 };
	plist[1610] = { 117, 120, 119 };
	plist[1611] = { 117, 119, 121 };
	plist[1612] = { 333, 120, 117 };
	plist[1613] = { 116, 140, 137 };
	plist[1614] = { 310, 116, 137 };
	plist[1615] = { 115, 140, 116 };
	plist[1616] = { 114, 230, 229 };
	plist[1617] = { 112, 115, 116 };
	plist[1618] = { 112, 225, 220 };
	plist[1619] = { 112, 240, 225 };
	plist[1620] = { 111, 139, 140 };
	plist[1621] = { 109, 173, 230 };
	plist[1622] = { 108, 274, 203 };
	plist[1623] = { 108, 173, 109 };
	plist[1624] = { 108, 203, 173 };
	plist[1625] = { 108, 109, 110 };
	plist[1626] = { 334, 111, 140 };
	plist[1627] = { 332, 111, 102 };
	plist[1628] = { 102, 111, 107 };
	plist[1629] = { 101, 110, 109 };
	plist[1630] = { 90, 2, 185 };
	plist[1631] = { 90, 94, 188 };
	plist[1632] = { 88, 100, 99 };
	plist[1633] = { 88, 99, 89 };
	plist[1634] = { 88, 98, 93 };
	plist[1635] = { 88, 89, 98 };
	plist[1636] = { 88, 93, 182 };
	plist[1637] = { 87, 98, 89 };
	plist[1638] = { 86, 96, 98 };
	plist[1639] = { 86, 98, 87 };
	plist[1640] = { 86, 121, 117 };
	plist[1641] = { 86, 87, 121 };
	plist[1642] = { 85, 93, 322 };
	plist[1643] = { 85, 182, 93 };
	plist[1644] = { 84, 108, 110 };
	plist[1645] = { 84, 106, 110 };
	plist[1646] = { 84, 331, 106 };
	plist[1647] = { 83, 102, 106 };
	plist[1648] = { 83, 332, 102 };
	plist[1649] = { 322, 83, 331 };
	plist[1650] = { 322, 331, 97 };
	plist[1651] = { 82, 274, 108 };
	plist[1652] = { 82, 108, 84 };
	plist[1653] = { 274, 273, 244 };
	plist[1654] = { 78, 81, 286 };
	plist[1655] = { 78, 183, 85 };
	plist[1656] = { 78, 85, 97 };
	plist[1657] = { 78, 97, 81 };
	plist[1658] = { 77, 87, 121 };
	plist[1659] = { 77, 121, 119 };
	plist[1660] = { 77, 119, 147 };
	plist[1661] = { 150, 100, 76 };
	plist[1662] = { 76, 77, 147 };
	plist[1663] = { 76, 147, 150 };
	plist[1664] = { 76, 87, 77 };
	plist[1665] = { 150, 179, 100 };
	plist[1666] = { 100, 89, 76 };
	plist[1667] = { 76, 89, 87 };
	plist[1668] = { 100, 99, 89 };
	plist[1669] = { 75, 109, 230 };
	plist[1670] = { 75, 107, 334 };
	plist[1671] = { 75, 101, 107 };
	plist[1672] = { 75, 101, 109 };
	plist[1673] = { 325, 114, 113 };
	plist[1674] = { 325, 115, 113 };
	plist[1675] = { 115, 334, 140 };
	plist[1676] = { 72, 129, 118 };
	plist[1677] = { 71, 124, 333 };
	plist[1678] = { 71, 118, 333 };
	plist[1679] = { 71, 72, 118 };
	plist[1680] = { 70, 124, 71 };
	plist[1681] = { 70, 72, 71 };
	plist[1682] = { 141, 163, 69 };
	plist[1683] = { 69, 72, 70 };
	plist[1684] = { 69, 163, 129 };
	plist[1685] = { 69, 129, 72 };
	plist[1686] = { 141, 193, 163 };
	plist[1687] = { 68, 141, 69 };
	plist[1688] = { 68, 141, 133 };
	plist[1689] = { 68, 133, 70 };
	plist[1690] = { 68, 69, 70 };
	plist[1691] = { 66, 157, 67 };
	plist[1692] = { 66, 281, 157 };
	plist[1693] = { 65, 160, 282 };
	plist[1694] = { 282, 66, 330 };
	plist[1695] = { 282, 330, 65 };
	plist[1696] = { 156, 335, 142 };
	plist[1697] = { 62, 64, 330 };
	plist[1698] = { 156, 64, 335 };
	plist[1699] = { 63, 157, 170 };
	plist[1700] = { 63, 67, 64 };
	plist[1701] = { 63, 64, 156 };
	plist[1702] = { 63, 67, 157 };
	plist[1703] = { 62, 146, 64 };
	plist[1704] = { 62, 149, 151 };
	plist[1705] = { 62, 151, 146 };
	plist[1706] = { 62, 330, 65 };
	plist[1707] = { 61, 160, 65 };
	plist[1708] = { 65, 159, 62 };
	plist[1709] = { 65, 61, 159 };
	plist[1710] = { 61, 159, 160 };
	plist[1711] = { 60, 144, 165 };
	plist[1712] = { 60, 168, 144 };
	plist[1713] = { 59, 199, 200 };
	plist[1714] = { 59, 168, 60 };
	plist[1715] = { 59, 60, 200 };
	plist[1716] = { 58, 239, 172 };
	plist[1717] = { 58, 200, 238 };
	plist[1718] = { 58, 238, 239 };
	plist[1719] = { 57, 60, 165 };
	plist[1720] = { 57, 165, 162 };
	plist[1721] = { 57, 200, 60 };
	plist[1722] = { 57, 162, 171 };
	plist[1723] = { 57, 58, 200 };
	plist[1724] = { 56, 57, 171 };
	plist[1725] = { 56, 58, 57 };
	plist[1726] = { 56, 172, 171 };
	plist[1727] = { 56, 58, 172 };
	plist[1728] = { 54, 149, 178 };
	plist[1729] = { 54, 178, 55 };
	plist[1730] = { 54, 10, 339 };
	plist[1731] = { 54, 55, 10 };
	plist[1732] = { 53, 21, 191 };
	plist[1733] = { 53, 192, 0 };
	plist[1734] = { 53, 0, 17 };
	plist[1735] = { 53, 17, 21 };
	plist[1736] = { 52, 201, 205 };
	plist[1737] = { 52, 205, 204 };
	plist[1738] = { 52, 221, 201 };
	plist[1739] = { 51, 208, 249 };
	plist[1740] = { 50, 52, 204 };
	plist[1741] = { 50, 51, 204 };
	plist[1742] = { 231, 235, 302 };
	plist[1743] = { 208, 51, 329 };
	plist[1744] = { 208, 329, 46 };
	plist[1745] = { 208, 46, 302 };
	plist[1746] = { 44, 48, 219 };
	plist[1747] = { 44, 218, 231 };
	plist[1748] = { 44, 231, 46 };
	plist[1749] = { 44, 46, 48 };
	plist[1750] = { 43, 248, 244 };
	plist[1751] = { 43, 244, 265 };
	plist[1752] = { 43, 247, 248 };
	plist[1753] = { 326, 247, 43 };
	plist[1754] = { 42, 337, 254 };
	plist[1755] = { 42, 247, 243 };
	plist[1756] = { 42, 326, 247 };
	plist[1757] = { 254, 281, 297 };
	plist[1758] = { 38, 326, 42 };
	plist[1759] = { 328, 279, 275 };
	plist[1760] = { 38, 42, 254 };
	plist[1761] = { 36, 262, 40 };
	plist[1762] = { 263, 262, 36 };
	plist[1763] = { 275, 36, 38 };
	plist[1764] = { 36, 40, 38 };
	plist[1765] = { 32, 94, 289 };
	plist[1766] = { 32, 260, 253 };
	plist[1767] = { 32, 253, 188 };
	plist[1768] = { 32, 188, 94 };
	plist[1769] = { 32, 289, 291 };
	plist[1770] = { 30, 138, 196 };
	plist[1771] = { 30, 196, 214 };
	plist[1772] = { 30, 291, 138 };
	plist[1773] = { 30, 32, 291 };
	plist[1774] = { 30, 260, 32 };
	plist[1775] = { 30, 214, 260 };
	plist[1776] = { 21, 17, 19 };
	plist[1777] = { 19, 17, 16 };
	plist[1778] = { 15, 12, 13 };
	plist[1779] = { 10, 5, 7 };
	plist[1780] = { 10, 7, 9 };
	plist[1781] = { 10, 55, 5 };
	plist[1782] = { 4, 185, 2 };
 // line elements
	plist[1783] = { 333, 71 };
	plist[1784] = { 330, 67 };
	plist[1785] = { 40, 262 };
	plist[1786] = { 36, 262 };
	plist[1787] = { 36, 275 };
	plist[1788] = { 160, 279 };
	plist[1789] = { 141, 193 };
	plist[1790] = { 68, 141 };
	plist[1791] = { 68, 133 };
	plist[1792] = { 68, 70 };
	plist[1793] = { 243, 247 };
	plist[1794] = { 172, 235 };
	plist[1795] = { 149, 179 };
	plist[1796] = { 117, 121 };
	plist[1797] = { 115, 116 };
	plist[1798] = { 113, 115 };
	plist[1799] = { 88, 93 };
	plist[1800] = { 88, 99 };
	plist[1801] = { 85, 93 };
	plist[1802] = { 325, 75 };
	plist[1803] = { 324, 84 };
	plist[1804] = { 324, 81 };
	plist[1805] = { 84, 110 };
	plist[1806] = { 81, 286 };
	plist[1807] = { 78, 85 };
	plist[1808] = { 78, 286 };
	plist[1809] = { 333, 117 };
	plist[1810] = { 61, 160 };
	plist[1811] = { 61, 65 };
	plist[1812] = { 61, 159 };
	plist[1813] = { 59, 200 };
	plist[1814] = { 326, 247 };
	plist[1815] = { 193, 198 };
	plist[1816] = { 315, 59 };
	plist[1817] = { 315, 63 };
	plist[1818] = { 101, 110 };
	plist[1819] = { 99, 100 };
	plist[1820] = { 87, 121 };
	plist[1821] = { 87, 89 };
	plist[1822] = { 325, 113 };
	plist[1823] = { 70, 71 };
	plist[1824] = { 63, 67 };
	plist[1825] = { 231, 235 };
	plist[1826] = { 50, 204 };
	plist[1827] = { 323, 264 };
	plist[1828] = { 323, 286 };
	plist[1829] = { 133, 137 };
	plist[1830] = { 100, 179 };
	plist[1831] = { 50, 48 };
	plist[1832] = { 304, 219 };
	plist[1833] = { 304, 220 };
	plist[1834] = { 299, 262 };
	plist[1835] = { 299, 264 };
	plist[1836] = { 89, 99 };
	plist[1837] = { 330, 65 };
	plist[1838] = { 58, 200 };
	plist[1839] = { 48, 219 };
	plist[1840] = { 44, 219 };
	plist[1841] = { 44, 231 };
	plist[1842] = { 204, 243 };
	plist[1843] = { 149, 159 };
	plist[1844] = { 116, 137 };
	plist[1845] = { 320, 220 };
	plist[1846] = { 320, 113 };
	plist[1847] = { 171, 198 };
	plist[1848] = { 326, 40 };
	plist[1849] = { 275, 279 };
	plist[1850] = { 75, 101 };
	plist[1851] = { 56, 172 };
	plist[1852] = { 56, 58 };
	plist[1853] = { 56, 171 };
	plist[1854] = { 210, 256 };
	plist[1855] = { 210, 211 };
	plist[1856] = { 190, 15 };
	plist[1857] = { 257, 276 };
	plist[1858] = { 90, 94 };
	plist[1859] = { 90, 185 };
	plist[1860] = { 232, 241 };
	plist[1861] = { 211, 232 };
	plist[1862] = { 94, 289 };
	plist[1863] = { 241, 15 };
	plist[1864] = { 289, 291 };
	plist[1865] = { 284, 4 };
	plist[1866] = { 134, 138 };
	plist[1867] = { 134, 190 };
	plist[1868] = { 138, 291 };
	plist[1869] = { 276, 284 };
	plist[1870] = { 256, 257 };
	plist[1871] = { 185, 4 };

vset.AddPlist( plist.begin(), plist.end());

//'pfverts' neighbors of the faces of each element
deque<vector<int64_t> >  pfverts( 1872 );
	pfverts[0] = { -1, 419, 446, 444 };
	pfverts[1] = { 178, 909, 1196, 2 };
	pfverts[2] = { 1210, 1, 1219, 161 };
	pfverts[3] = { 739, 1328, 33, 619 };
	pfverts[4] = { 625, 616, 629, 602 };
	pfverts[5] = { 1314, 814, 799, 679 };
	pfverts[6] = { 26, 1323, 1360, 678 };
	pfverts[7] = { 938, 25, 843, 868 };
	pfverts[8] = { 960, 948, 963, 959 };
	pfverts[9] = { -1, 1022, 1029, 61 };
	pfverts[10] = { 708, -1, 1322, 1324 };
	pfverts[11] = { 707, 30, 1266, 1269 };
	pfverts[12] = { -1, 47, 828, 156 };
	pfverts[13] = { 604, 321, 600, 632 };
	pfverts[14] = { -1, -1, 67, 765 };
	pfverts[15] = { -1, -1, 209, 1073 };
	pfverts[16] = { 883, 880, 865, 886 };
	pfverts[17] = { 183, -1, 1078, 415 };
	pfverts[18] = { 344, 345, 387, 19 };
	pfverts[19] = { 513, 91, 18, 391 };
	pfverts[20] = { 491, 905, 218, 890 };
	pfverts[21] = { 641, -1, 22, 1332 };
	pfverts[22] = { 21, 629, 627, 602 };
	pfverts[23] = { 185, 173, 945, 957 };
	pfverts[24] = { 951, 957, 547, 996 };
	pfverts[25] = { 873, 935, 7, 922 };
	pfverts[26] = { 6, 786, 784, 1338 };
	pfverts[27] = { 372, 376, 1260, 1297 };
	pfverts[28] = { 1258, 1253, 1263, 1252 };
	pfverts[29] = { 1248, 214, 1252, 1246 };
	pfverts[30] = { 665, 1272, 677, 11 };
	pfverts[31] = { -1, 124, 118, 542 };
	pfverts[32] = { -1, 1332, 619, 1362 };
	pfverts[33] = { 764, 3, 637, 772 };
	pfverts[34] = { 702, -1, 1182, 1165 };
	pfverts[35] = { 697, 1184, -1, 36 };
	pfverts[36] = { 302, 1186, 35, -1 };
	pfverts[37] = { -1, -1, 313, 1180 };
	pfverts[38] = { 693, 39, 83, 690 };
	pfverts[39] = { 774, 38, 40, -1 };
	pfverts[40] = { 768, 686, -1, 39 };
	pfverts[41] = { 1343, -1, 864, 866 };
	pfverts[42] = { 267, 43, 925, -1 };
	pfverts[43] = { 721, 866, 867, 42 };
	pfverts[44] = { 529, -1, 535, 530 };
	pfverts[45] = { 485, 513, 97, 91 };
	pfverts[46] = { 519, 532, -1, 80 };
	pfverts[47] = { 12, -1, 836, 837 };
	pfverts[48] = { 154, 1338, 787, 785 };
	pfverts[49] = { 50, 785, 1309, -1 };
	pfverts[50] = { 125, 837, 49, -1 };
	pfverts[51] = { 633, 323, -1, 600 };
	pfverts[52] = { 1348, 53, -1, 620 };
	pfverts[53] = { 623, 600, 52, 632 };
	pfverts[54] = { 296, 209, -1, 943 };
	pfverts[55] = { 213, 589, 281, 1080 };
	pfverts[56] = { 1242, -1, 1243, 988 };
	pfverts[57] = { 980, -1, 981, 1003 };
	pfverts[58] = { 275, 1003, -1, -1 };
	pfverts[59] = { 949, 1016, 1020, 946 };
	pfverts[60] = { 280, 168, -1, 1000 };
	pfverts[61] = { 9, 1006, 984, 973 };
	pfverts[62] = { 649, 63, 250, 830 };
	pfverts[63] = { 62, 1014, 252, 152 };
	pfverts[64] = { 1013, 1027, 1025, 250 };
	pfverts[65] = { 270, 1025, -1, 1345 };
	pfverts[66] = { 697, 967, 1345, 1119 };
	pfverts[67] = { -1, 14, 759, 761 };
	pfverts[68] = { 916, 245, -1, -1 };
	pfverts[69] = { 1342, 345, 341, 70 };
	pfverts[70] = { 293, 69, 1083, 393 };
	pfverts[71] = { 1334, 267, 72, 721 };
	pfverts[72] = { 143, 71, 926, 718 };
	pfverts[73] = { 768, -1, 720, 774 };
	pfverts[74] = { 776, 75, 766, 775 };
	pfverts[75] = { 147, 926, -1, 74 };
	pfverts[76] = { 205, 396, -1, 390 };
	pfverts[77] = { 98, 115, 210, 358 };
	pfverts[78] = { 410, 373, 1255, 1256 };
	pfverts[79] = { 456, 108, 217, 431 };
	pfverts[80] = { -1, 46, 97, 518 };
	pfverts[81] = { 123, 100, 226, 555 };
	pfverts[82] = { 642, 149, 232, 612 };
	pfverts[83] = { 38, 140, 235, 661 };
	pfverts[84] = { 290, 1213, -1, 1192 };
	pfverts[85] = { 209, 296, -1, 177 };
	pfverts[86] = { 1086, 302, -1, 1157 };
	pfverts[87] = { 1142, 308, -1, 192 };
	pfverts[88] = { 604, 326, 316, 321 };
	pfverts[89] = { 842, 1127, 328, 333 };
	pfverts[90] = { 335, -1, 860, 855 };
	pfverts[91] = { 223, 19, 45, 518 };
	pfverts[92] = { 1042, 173, 353, 379 };
	pfverts[93] = { 96, 208, 360, 1295 };
	pfverts[94] = { 488, 111, 1340, 112 };
	pfverts[95] = { 362, 1281, 1290, 1339 };
	pfverts[96] = { 364, 182, 381, 93 };
	pfverts[97] = { 531, 45, -1, 80 };
	pfverts[98] = { 389, 77, 398, 1355 };
	pfverts[99] = { 415, 392, 406, 405 };
	pfverts[100] = { 81, 464, 417, 437 };
	pfverts[101] = { 422, 1357, -1, 443 };
	pfverts[102] = { 1127, 477, 426, 450 };
	pfverts[103] = { 452, -1, 106, 440 };
	pfverts[104] = { 449, 554, 439, 105 };
	pfverts[105] = { 104, 1288, 445, 430 };
	pfverts[106] = { 1344, -1, 462, 103 };
	pfverts[107] = { 447, 461, 448, 1344 };
	pfverts[108] = { 79, 517, 468, 492 };
	pfverts[109] = { 506, 895, 909, 913 };
	pfverts[110] = { 521, -1, 890, 1335 };
	pfverts[111] = { 505, 430, 496, 94 };
	pfverts[112] = { 94, 1288, 1287, 1302 };
	pfverts[113] = { -1, -1, 525, 535 };
	pfverts[114] = { 533, 487, 530, 519 };
	pfverts[115] = { 77, 122, 540, 562 };
	pfverts[116] = { 419, -1, 545, 569 };
	pfverts[117] = { 992, 423, 550, 576 };
	pfverts[118] = { 1004, 31, 568, 174 };
	pfverts[119] = { 575, 363, 1339, 1289 };
	pfverts[120] = { 554, 1286, 1339, 574 };
	pfverts[121] = { 557, 174, 568, 571 };
	pfverts[122] = { 115, 389, 588, 597 };
	pfverts[123] = { 584, 81, 591, 596 };
	pfverts[124] = { 594, 31, -1, 586 };
	pfverts[125] = { 50, -1, 639, 625 };
	pfverts[126] = { 314, 788, 607, 631 };
	pfverts[127] = { 630, 734, 621, 128 };
	pfverts[128] = { 127, 1312, 624, 611 };
	pfverts[129] = { 629, -1, 618, 628 };
	pfverts[130] = { -1, -1, 638, 642 };
	pfverts[131] = { 252, 695, 646, 678 };
	pfverts[132] = { 1146, 165, 656, 680 };
	pfverts[133] = { 135, 234, 663, 674 };
	pfverts[134] = { 665, 1284, 660, 1341 };
	pfverts[135] = { 645, 197, 682, 133 };
	pfverts[136] = { 661, 690, 707, 677 };
	pfverts[137] = { 769, 1179, -1, 1333 };
	pfverts[138] = { 714, 674, 710, 139 };
	pfverts[139] = { 138, 691, 712, 705 };
	pfverts[140] = { 83, 775, 146, 740 };
	pfverts[141] = { 1311, 1327, 726, 1328 };
	pfverts[142] = { 917, 604, 730, 753 };
	pfverts[143] = { 147, 750, 743, 72 };
	pfverts[144] = { 744, 1341, 1315, 1319 };
	pfverts[145] = { 734, 1311, 1341, 751 };
	pfverts[146] = { 140, 718, 716, 746 };
	pfverts[147] = { 75, 143, 770, 773 };
	pfverts[148] = { 762, -1, 763, 773 };
	pfverts[149] = { 82, 156, 780, 803 };
	pfverts[150] = { 1273, 1338, 651, 1329 };
	pfverts[151] = { 815, 1019, 1031, 1035 };
	pfverts[152] = { 796, 251, 63, 798 };
	pfverts[153] = { 814, 611, 805, 1317 };
	pfverts[154] = { 1317, 724, 48, 1329 };
	pfverts[155] = { -1, 156, 828, 822 };
	pfverts[156] = { 642, 149, 12, 155 };
	pfverts[157] = { 824, 252, 1014, 835 };
	pfverts[158] = { -1, 778, 247, 829 };
	pfverts[159] = { 326, 917, 845, 858 };
	pfverts[160] = { 338, 862, 335, 860 };
	pfverts[161] = { 874, 2, 885, 258 };
	pfverts[162] = { 865, 202, 883, 259 };
	pfverts[163] = { 477, 842, 894, 908 };
	pfverts[164] = { 905, 491, 910, 261 };
	pfverts[165] = { 745, 132, 754, 882 };
	pfverts[166] = { 267, 1334, 755, 939 };
	pfverts[167] = { 1090, 1031, 948, 960 };
	pfverts[168] = { 578, 1006, 60, 994 };
	pfverts[169] = { 956, 189, 964, 269 };
	pfverts[170] = { 1233, 992, 975, 985 };
	pfverts[171] = { 1022, 276, -1, 983 };
	pfverts[172] = { -1, 989, 988, 1244 };
	pfverts[173] = { 567, 92, 577, 23 };
	pfverts[174] = { 1009, 542, 118, 121 };
	pfverts[175] = { 788, 972, 1018, 1030 };
	pfverts[176] = { 1027, 802, 1032, 283 };
	pfverts[177] = { 1060, 85, -1, 1101 };
	pfverts[178] = { 350, 1, 1045, 1066 };
	pfverts[179] = { 291, 289, 1055, 1052 };
	pfverts[180] = { 1065, 1097, 1058, 181 };
	pfverts[181] = { 180, 188, 1063, 1053 };
	pfverts[182] = { 1041, 96, 1068, 1052 };
	pfverts[183] = { 1075, -1, 1056, 17 };
	pfverts[184] = { 201, 1166, -1, 196 };
	pfverts[185] = { 23, 1042, 1093, 1112 };
	pfverts[186] = { 189, 295, 1099, 1109 };
	pfverts[187] = { 1157, 1166, 1156, 1105 };
	pfverts[188] = { 1166, 181, 1124, 1097 };
	pfverts[189] = { 1103, 169, 1114, 186 };
	pfverts[190] = { 1123, 1085, -1, 1118 };
	pfverts[191] = { 423, 1233, 1130, 1135 };
	pfverts[192] = { 1163, 87, -1, 1207 };
	pfverts[193] = { 653, 1090, 1149, 1169 };
	pfverts[194] = { 648, 303, 1155, 304 };
	pfverts[195] = { 1168, 1203, 1161, 196 };
	pfverts[196] = { 195, 1105, 184, 1167 };
	pfverts[197] = { 1159, 135, 1171, 1155 };
	pfverts[198] = { 882, 1146, 1199, 1218 };
	pfverts[199] = { 202, 307, 1205, 1215 };
	pfverts[200] = { 1065, 1217, 1053, 1049 };
	pfverts[201] = { 1053, 184, 1351, 1203 };
	pfverts[202] = { 1209, 162, 1220, 199 };
	pfverts[203] = { 1190, 1224, -1, 1230 };
	pfverts[204] = { 972, 314, 1235, 1239 };
	pfverts[205] = { 206, -1, -1, 76 };
	pfverts[206] = { -1, 388, 205, -1 };
	pfverts[207] = { 537, 279, 211, 1036 };
	pfverts[208] = { 93, 1041, 212, 397 };
	pfverts[209] = { 85, 15, 54, 1331 };
	pfverts[210] = { 77, 398, 212, 213 };
	pfverts[211] = { 1040, 207, 212, 213 };
	pfverts[212] = { 208, 360, 210, 211 };
	pfverts[213] = { 543, 55, 210, 211 };
	pfverts[214] = { 1257, 29, -1, 1274 };
	pfverts[215] = { 432, 465, 216, 1134 };
	pfverts[216] = { 471, 215, -1, 217 };
	pfverts[217] = { 216, 433, -1, 79 };
	pfverts[218] = { 221, 20, 486, 490 };
	pfverts[219] = { 514, -1, 1342, 341 };
	pfverts[220] = { 1342, 889, 1335, 902 };
	pfverts[221] = { 489, 890, 218, 223 };
	pfverts[222] = { 875, 339, 490, 902 };
	pfverts[223] = { 91, 902, 1335, 221 };
	pfverts[224] = { 556, 416, 228, 990 };
	pfverts[225] = { -1, 274, 991, 229 };
	pfverts[226] = { 81, 418, 227, 228 };
	pfverts[227] = { 591, -1, 229, 226 };
	pfverts[228] = { 557, 226, 1001, 224 };
	pfverts[229] = { 592, 225, 1001, 227 };
	pfverts[230] = { 613, 777, 231, 320 };
	pfverts[231] = { 783, 230, -1, 232 };
	pfverts[232] = { 231, 614, -1, 82 };
	pfverts[233] = { 662, 722, 236, 1139 };
	pfverts[234] = { 133, 1140, 237, 691 };
	pfverts[235] = { 83, 694, 237, 238 };
	pfverts[236] = { 1143, 233, 237, 238 };
	pfverts[237] = { 234, 663, 235, 236 };
	pfverts[238] = { 746, 1181, 235, 236 };
	pfverts[239] = { 735, 620, 241, 914 };
	pfverts[240] = { 739, 762, 241, 764 };
	pfverts[241] = { 239, 240, 243, 242 };
	pfverts[242] = { 1348, 764, 241, 765 };
	pfverts[243] = { 915, 244, 245, 241 };
	pfverts[244] = { 932, -1, 243, 762 };
	pfverts[245] = { 68, 243, -1, 765 };
	pfverts[246] = { 797, 644, 250, 1013 };
	pfverts[247] = { 826, 833, 834, 158 };
	pfverts[248] = { 1183, 1176, 701, 825 };
	pfverts[249] = { 678, 798, 812, 670 };
	pfverts[250] = { 246, 796, 64, 62 };
	pfverts[251] = { 809, 152, 778, 802 };
	pfverts[252] = { 131, 63, 157, 798 };
	pfverts[253] = { 849, 914, 254, 332 };
	pfverts[254] = { -1, 253, 850, 256 };
	pfverts[255] = { -1, 860, -1, 850 };
	pfverts[256] = { 916, -1, -1, 254 };
	pfverts[257] = { 1193, 874, 259, 877 };
	pfverts[258] = { 903, 904, 161, 875 };
	pfverts[259] = { 1190, 162, -1, 257 };
	pfverts[260] = { 899, 839, 263, 465 };
	pfverts[261] = { 264, 466, -1, 164 };
	pfverts[262] = { 265, 471, -1, 467 };
	pfverts[263] = { 469, 260, 264, 265 };
	pfverts[264] = { 901, 263, -1, 261 };
	pfverts[265] = { 841, 262, -1, 263 };
	pfverts[266] = { 1343, 1330, -1, 719 };
	pfverts[267] = { 166, -1, 71, 42 };
	pfverts[268] = { 953, 1013, 270, 1084 };
	pfverts[269] = { 271, 1085, -1, 169 };
	pfverts[270] = { 1087, 268, 271, 65 };
	pfverts[271] = { 955, 270, -1, 269 };
	pfverts[272] = { 967, 830, 1025, -1 };
	pfverts[273] = { 979, 990, 275, 1245 };
	pfverts[274] = { -1, 275, -1, 225 };
	pfverts[275] = { -1, 273, 58, 274 };
	pfverts[276] = { 171, 987, 989, 984 };
	pfverts[277] = { 961, 538, 541, 280 };
	pfverts[278] = { 541, -1, 281, 966 };
	pfverts[279] = { 207, 543, 541, 563 };
	pfverts[280] = { 60, 277, -1, 561 };
	pfverts[281] = { 543, 539, 55, 278 };
	pfverts[282] = { 1023, 968, 285, 777 };
	pfverts[283] = { 286, 778, -1, 176 };
	pfverts[284] = { 287, 783, -1, 779 };
	pfverts[285] = { 781, 282, 286, 287 };
	pfverts[286] = { 1024, 285, -1, 283 };
	pfverts[287] = { 971, 284, -1, 285 };
	pfverts[288] = { 1051, 1189, 291, 339 };
	pfverts[289] = { 290, 1213, 179, 292 };
	pfverts[290] = { 84, 1081, 289, 293 };
	pfverts[291] = { 343, 288, 179, 293 };
	pfverts[292] = { 1075, 1072, 289, 1052 };
	pfverts[293] = { 1195, 70, 290, 291 };
	pfverts[294] = { 1098, 1036, 297, 941 };
	pfverts[295] = { 186, 956, 298, 1121 };
	pfverts[296] = { 85, 1122, 298, 54 };
	pfverts[297] = { 944, 294, 298, 943 };
	pfverts[298] = { 295, 1099, 296, 297 };
	pfverts[299] = { 985, 1232, 1242, 1233 };
	pfverts[300] = { 1154, 1084, 303, 644 };
	pfverts[301] = { 1165, 304, 1155, 1105 };
	pfverts[302] = { 86, 36, 304, 305 };
	pfverts[303] = { 647, 300, 194, 305 };
	pfverts[304] = { 1182, 302, 194, 301 };
	pfverts[305] = { 1089, 697, 302, 303 };
	pfverts[306] = { 1204, 1139, 309, 928 };
	pfverts[307] = { 199, 878, 310, 1227 };
	pfverts[308] = { 87, 1226, 310, 311 };
	pfverts[309] = { 864, 306, 310, 311 };
	pfverts[310] = { 307, 1205, 308, 309 };
	pfverts[311] = { 1145, 1330, 308, 309 };
	pfverts[312] = { -1, 393, 341, 1076 };
	pfverts[313] = { -1, -1, 686, 37 };
	pfverts[314] = { 126, 322, 317, 204 };
	pfverts[315] = { 1234, 327, 319, 606 };
	pfverts[316] = { 329, 88, 318, 319 };
	pfverts[317] = { 314, 1235, 318, 319 };
	pfverts[318] = { 324, 1128, 316, 317 };
	pfverts[319] = { 315, 610, 316, 317 };
	pfverts[320] = { 230, 968, -1, 1240 };
	pfverts[321] = { 88, 13, 325, 324 };
	pfverts[322] = { 314, 631, 323, 1240 };
	pfverts[323] = { 322, -1, 51, 324 };
	pfverts[324] = { 318, -1, 323, 321 };
	pfverts[325] = { 334, -1, 321, 332 };
	pfverts[326] = { 159, 334, 329, 88 };
	pfverts[327] = { 315, 1129, 331, 844 };
	pfverts[328] = { 1131, 89, 330, 331 };
	pfverts[329] = { 326, 316, 330, 331 };
	pfverts[330] = { 336, 1346, 328, 329 };
	pfverts[331] = { 327, 848, 328, 329 };
	pfverts[332] = { 253, -1, 325, 620 };
	pfverts[333] = { 89, 1126, 337, 859 };
	pfverts[334] = { 326, 858, 338, 325 };
	pfverts[335] = { 863, 337, 160, 90 };
	pfverts[336] = { -1, 330, 337, 338 };
	pfverts[337] = { -1, 333, 335, 336 };
	pfverts[338] = { 334, -1, 160, 336 };
	pfverts[339] = { 222, 366, 343, 288 };
	pfverts[340] = { 1052, 364, 344, 392 };
	pfverts[341] = { 219, 69, 411, 312 };
	pfverts[342] = { 1076, 395, 411, -1 };
	pfverts[343] = { 339, 291, 344, 345 };
	pfverts[344] = { 368, 340, 18, 343 };
	pfverts[345] = { 69, 902, 18, 343 };
	pfverts[346] = { 1251, 1302, 544, 1290 };
	pfverts[347] = { 546, 1248, 1257, 1249 };
	pfverts[348] = { 473, 1336, 373, 1337 };
	pfverts[349] = { 474, 1257, 1250, 1249 };
	pfverts[350] = { 178, 380, 355, 891 };
	pfverts[351] = { 478, 1289, 356, 548 };
	pfverts[352] = { 479, 549, 354, 1198 };
	pfverts[353] = { 1048, 92, 354, 355 };
	pfverts[354] = { 551, 352, 353, 356 };
	pfverts[355] = { 353, 382, 356, 350 };
	pfverts[356] = { 351, 354, 355, 480 };
	pfverts[357] = { 1296, 566, 363, 1339 };
	pfverts[358] = { 562, 372, 1297, 77 };
	pfverts[359] = { 563, 363, 567, 537 };
	pfverts[360] = { 212, 93, 1291, 361 };
	pfverts[361] = { 360, 383, 564, 537 };
	pfverts[362] = { 95, 1262, 1296, 1260 };
	pfverts[363] = { 357, 119, 359, 564 };
	pfverts[364] = { 96, 402, 368, 340 };
	pfverts[365] = { 513, 373, 1293, 485 };
	pfverts[366] = { 339, 888, 370, 380 };
	pfverts[367] = { 1279, 369, 1292, 370 };
	pfverts[368] = { 364, 344, 369, 370 };
	pfverts[369] = { 367, 377, 1293, 368 };
	pfverts[370] = { 366, 384, 367, 368 };
	pfverts[371] = { 407, 1295, 405, 375 };
	pfverts[372] = { 27, 1277, 1285, 358 };
	pfverts[373] = { 348, 1300, 78, 365 };
	pfverts[374] = { 377, 1279, 1262, 1298 };
	pfverts[375] = { 1356, 376, 371, 404 };
	pfverts[376] = { 27, 1285, 1299, 375 };
	pfverts[377] = { 1301, 369, 374, 402 };
	pfverts[378] = { 488, 1281, 1279, 1289 };
	pfverts[379] = { 92, 537, 383, 1067 };
	pfverts[380] = { 350, 1066, 384, 366 };
	pfverts[381] = { 1068, 96, 383, 384 };
	pfverts[382] = { 1278, 355, 383, 384 };
	pfverts[383] = { 361, 379, 381, 382 };
	pfverts[384] = { 380, 370, 381, 382 };
	pfverts[385] = { -1, -1, -1, 1082 };  // tetra spanning corner
	pfverts[386] = { -1, 391, 1082, 512 };
	pfverts[387] = { 18, 391, 411, 395 };
	pfverts[388] = { 390, 206, 1038, 587 };
	pfverts[389] = { 98, 122, -1, 583 };
	pfverts[390] = { 388, 401, 76, 583 };
	pfverts[391] = { 19, 394, 386, 387 };
	pfverts[392] = { 99, 1075, 395, 340 };
	pfverts[393] = { 1077, 312, 70, 1191 };
	pfverts[394] = { 412, -1, 391, 406 };
	pfverts[395] = { 392, 406, 387, 342 };
	pfverts[396] = { 76, 414, 399, -1 };
	pfverts[397] = { 1056, 1078, 400, 208 };
	pfverts[398] = { 210, 98, 400, 401 };
	pfverts[399] = { 396, 1353, 400, 401 };
	pfverts[400] = { 413, 397, 398, 399 };
	pfverts[401] = { 1038, 390, 398, 399 };
	pfverts[402] = { 364, 377, 405, 1295 };
	pfverts[403] = { -1, 1277, 1285, 1246 };
	pfverts[404] = { -1, 375, 1294, 1299 };
	pfverts[405] = { 99, 402, 1294, 371 };
	pfverts[406] = { 99, 409, 394, 395 };
	pfverts[407] = { 371, 415, 1078, 1079 };
	pfverts[408] = { -1, 1062, 1056, 1353 };
	pfverts[409] = { 415, -1, -1, 406 };
	pfverts[410] = { -1, 78, 412, 513 };
	pfverts[411] = { 341, 387, 1082, 342 };
	pfverts[412] = { 1294, 410, 394, -1 };
	pfverts[413] = { 400, 1078, 414, 1079 };
	pfverts[414] = { 396, -1, -1, 413 };
	pfverts[415] = { 99, 409, 17, 407 };
	pfverts[416] = { 224, 438, 418, 1245 };
	pfverts[417] = { 100, 440, -1, 418 };
	pfverts[418] = { 417, 416, -1, 226 };
	pfverts[419] = { 0, -1, 475, 116 };
	pfverts[420] = { 1249, 545, 544, 472 };
	pfverts[421] = { 1303, 422, 472, 1288 };
	pfverts[422] = { 101, 445, 475, 421 };
	pfverts[423] = { 191, 451, 427, 117 };
	pfverts[424] = { 548, 449, 428, 478 };
	pfverts[425] = { 549, 479, 429, 1129 };
	pfverts[426] = { 481, 102, 428, 429 };
	pfverts[427] = { 423, 550, 428, 429 };
	pfverts[428] = { 453, 424, 426, 427 };
	pfverts[429] = { 425, 1133, 426, 427 };
	pfverts[430] = { 449, 105, 434, 111 };
	pfverts[431] = { 492, 443, 435, 79 };
	pfverts[432] = { 493, 215, 436, 450 };
	pfverts[433] = { 217, 1344, 435, 436 };
	pfverts[434] = { 430, 497, 435, 436 };
	pfverts[435] = { 447, 431, 433, 434 };
	pfverts[436] = { 432, 454, 433, 434 };
	pfverts[437] = { 100, 444, 441, 555 };
	pfverts[438] = { 416, 556, 442, 451 };
	pfverts[439] = { 558, 104, 441, 442 };
	pfverts[440] = { 103, 417, 441, 442 };
	pfverts[441] = { 462, 437, 439, 440 };
	pfverts[442] = { 438, 455, 439, 440 };
	pfverts[443] = { 458, 431, 447, 101 };
	pfverts[444] = { 457, 437, 0, 459 };
	pfverts[445] = { 422, 105, 447, 448 };
	pfverts[446] = { 460, 0, 448, 459 };
	pfverts[447] = { 435, 443, 445, 107 };
	pfverts[448] = { 107, 446, 445, 462 };
	pfverts[449] = { 104, 430, 453, 424 };
	pfverts[450] = { 102, 432, 454, 1136 };
	pfverts[451] = { 423, 1135, 455, 438 };
	pfverts[452] = { 1137, 103, 454, 455 };
	pfverts[453] = { 449, 428, 454, 455 };
	pfverts[454] = { 436, 450, 452, 453 };
	pfverts[455] = { 451, 442, 452, 453 };
	pfverts[456] = { 534, 79, -1, 463 };
	pfverts[457] = { -1, -1, 464, 444 };
	pfverts[458] = { -1, 463, 461, 443 };
	pfverts[459] = { -1, 444, 446, 462 };
	pfverts[460] = { -1, -1, 461, 446 };
	pfverts[461] = { -1, 458, 107, 460 };
	pfverts[462] = { 106, 459, 448, 441 };
	pfverts[463] = { -1, 456, -1, 458 };
	pfverts[464] = { 584, -1, 100, 457 };
	pfverts[465] = { 215, 493, 469, 260 };
	pfverts[466] = { 261, 498, 470, 524 };
	pfverts[467] = { 262, 528, -1, -1 };
	pfverts[468] = { 526, 108, 470, 471 };
	pfverts[469] = { 465, 263, 470, 471 };
	pfverts[470] = { 494, 466, 468, 469 };
	pfverts[471] = { 262, 216, 468, 469 };
	pfverts[472] = { 474, 475, 420, 421 };
	pfverts[473] = { 1304, 1305, 348, 499 };
	pfverts[474] = { -1, 349, 476, 472 };
	pfverts[475] = { -1, 419, 472, 422 };
	pfverts[476] = { -1, 474, 1304, 1303 };
	pfverts[477] = { 163, 507, 481, 102 };
	pfverts[478] = { 424, 505, 482, 351 };
	pfverts[479] = { 425, 352, 483, 893 };
	pfverts[480] = { 356, 506, 482, 483 };
	pfverts[481] = { 477, 426, 482, 483 };
	pfverts[482] = { 511, 478, 480, 481 };
	pfverts[483] = { 479, 897, 480, 481 };
	pfverts[484] = { 1340, 489, 1292, 1293 };
	pfverts[485] = { 365, 499, 489, 45 };
	pfverts[486] = { 1292, 218, 510, 888 };
	pfverts[487] = { 890, 501, 498, 114 };
	pfverts[488] = { 505, 378, 94, 1292 };
	pfverts[489] = { 501, 485, 221, 484 };
	pfverts[490] = { 222, 903, 218, 888 };
	pfverts[491] = { 20, 164, 498, 508 };
	pfverts[492] = { 108, 500, 494, 431 };
	pfverts[493] = { 465, 432, 497, 507 };
	pfverts[494] = { 492, 470, 497, 495 };
	pfverts[495] = { 503, 498, 494, 496 };
	pfverts[496] = { 111, 495, 497, 508 };
	pfverts[497] = { 434, 494, 493, 496 };
	pfverts[498] = { 495, 491, 487, 466 };
	pfverts[499] = { 502, 473, 531, 485 };
	pfverts[500] = { 529, 1357, 504, 492 };
	pfverts[501] = { 487, 489, 503, 502 };
	pfverts[502] = { 530, 501, 504, 499 };
	pfverts[503] = { 495, 1340, 504, 501 };
	pfverts[504] = { 500, 503, 1305, 502 };
	pfverts[505] = { 111, 488, 510, 478 };
	pfverts[506] = { 511, 480, 891, 109 };
	pfverts[507] = { 477, 908, 509, 493 };
	pfverts[508] = { 496, 491, 509, 510 };
	pfverts[509] = { 910, 507, 508, 511 };
	pfverts[510] = { 508, 486, 511, 505 };
	pfverts[511] = { 506, 509, 510, 482 };
	pfverts[512] = { 386, 522, -1, 518 };
	pfverts[513] = { 365, 410, 45, 19 };
	pfverts[514] = { 1082, 219, -1, 522 };
	pfverts[515] = { -1, 523, -1, 516 };
	pfverts[516] = { -1, 527, 515, -1 };
	pfverts[517] = { 529, 536, 534, 108 };
	pfverts[518] = { 512, 520, 80, 91 };
	pfverts[519] = { 114, -1, 521, 46 };
	pfverts[520] = { 1335, 518, 521, 522 };
	pfverts[521] = { 110, -1, 520, 519 };
	pfverts[522] = { 512, 514, -1, 520 };
	pfverts[523] = { -1, 515, 534, 525 };
	pfverts[524] = { 466, 533, 526, -1 };
	pfverts[525] = { 523, 113, 527, 536 };
	pfverts[526] = { 468, 536, 528, 524 };
	pfverts[527] = { 516, -1, 528, 525 };
	pfverts[528] = { 467, 527, -1, 526 };
	pfverts[529] = { -1, 500, 44, 517 };
	pfverts[530] = { 114, 502, 44, 532 };
	pfverts[531] = { 532, -1, 97, 499 };
	pfverts[532] = { 46, 531, 530, -1 };
	pfverts[533] = { 114, -1, 535, 524 };
	pfverts[534] = { -1, 523, 456, 517 };
	pfverts[535] = { 44, 533, 113, 536 };
	pfverts[536] = { 517, 535, 525, 526 };
	pfverts[537] = { 207, 361, 359, 379 };
	pfverts[538] = { 561, 277, 563, 567 };
	pfverts[539] = { 582, -1, 587, 281 };
	pfverts[540] = { 588, 115, 565, 543 };
	pfverts[541] = { 941, 279, 278, 277 };
	pfverts[542] = { 561, 565, 31, 174 };
	pfverts[543] = { 281, 213, 540, 279 };
	pfverts[544] = { 346, 420, 1286, 1288 };
	pfverts[545] = { 116, 1308, 546, 420 };
	pfverts[546] = { 1307, -1, 545, 347 };
	pfverts[547] = { 946, 24, 994, 1010 };
	pfverts[548] = { 351, 575, 552, 424 };
	pfverts[549] = { 352, 425, 553, 993 };
	pfverts[550] = { 427, 117, 552, 553 };
	pfverts[551] = { 577, 354, 552, 553 };
	pfverts[552] = { 579, 548, 550, 551 };
	pfverts[553] = { 549, 997, 550, 551 };
	pfverts[554] = { 575, 120, 558, 104 };
	pfverts[555] = { 437, 569, 559, 81 };
	pfverts[556] = { 438, 224, 560, 576 };
	pfverts[557] = { 228, 121, 559, 560 };
	pfverts[558] = { 554, 439, 559, 560 };
	pfverts[559] = { 571, 555, 557, 558 };
	pfverts[560] = { 556, 581, 557, 558 };
	pfverts[561] = { 280, 1009, 538, 542 };
	pfverts[562] = { 115, 570, 566, 358 };
	pfverts[563] = { 565, 279, 359, 538 };
	pfverts[564] = { 1298, 1278, 363, 361 };
	pfverts[565] = { 540, 542, 563, 566 };
	pfverts[566] = { 574, 562, 357, 565 };
	pfverts[567] = { 580, 538, 359, 173 };
	pfverts[568] = { 1001, 118, 594, 121 };
	pfverts[569] = { 596, 555, 573, 116 };
	pfverts[570] = { 598, 1306, 572, 562 };
	pfverts[571] = { 121, 573, 559, 574 };
	pfverts[572] = { 570, 1308, 573, 574 };
	pfverts[573] = { 595, 569, 571, 572 };
	pfverts[574] = { 566, 572, 571, 120 };
	pfverts[575] = { 119, 554, 579, 548 };
	pfverts[576] = { 117, 556, 581, 1008 };
	pfverts[577] = { 551, 580, 173, 996 };
	pfverts[578] = { 1005, 1009, 1011, 168 };
	pfverts[579] = { 575, 552, 581, 580 };
	pfverts[580] = { 567, 577, 579, 1012 };
	pfverts[581] = { 560, 579, 576, 1012 };
	pfverts[582] = { 583, 588, -1, 539 };
	pfverts[583] = { 390, 389, -1, 582 };
	pfverts[584] = { 123, -1, 585, 464 };
	pfverts[585] = { 584, -1, -1, 593 };
	pfverts[586] = { 124, 588, -1, 597 };
	pfverts[587] = { 388, -1, 1080, 539 };
	pfverts[588] = { 122, 540, 586, 582 };
	pfverts[589] = { 1039, 1331, 55, 1073 };
	pfverts[590] = { 592, 594, -1, 599 };
	pfverts[591] = { 227, 123, 592, 593 };
	pfverts[592] = { 591, 229, -1, 590 };
	pfverts[593] = { 585, -1, -1, 591 };
	pfverts[594] = { 590, 124, 595, 568 };
	pfverts[595] = { 594, 573, 598, 599 };
	pfverts[596] = { 569, 123, -1, 599 };
	pfverts[597] = { 586, -1, 598, 122 };
	pfverts[598] = { 595, -1, 570, 597 };
	pfverts[599] = { 595, 596, -1, 590 };
	pfverts[600] = { 636, 51, 53, 13 };
	pfverts[601] = { 626, 1347, -1, 772 };
	pfverts[602] = { 639, 4, 22, 624 };
	pfverts[603] = { -1, 639, 1354, 1309 };
	pfverts[604] = { 88, 13, 608, 142 };
	pfverts[605] = { 728, 630, 609, 789 };
	pfverts[606] = { 729, 790, 610, 315 };
	pfverts[607] = { 792, 126, 609, 610 };
	pfverts[608] = { 604, 730, 609, 610 };
	pfverts[609] = { 634, 605, 607, 608 };
	pfverts[610] = { 606, 319, 607, 608 };
	pfverts[611] = { 630, 128, 615, 153 };
	pfverts[612] = { 803, 625, 616, 82 };
	pfverts[613] = { 804, 230, 617, 631 };
	pfverts[614] = { 232, 628, 616, 617 };
	pfverts[615] = { 611, 806, 616, 617 };
	pfverts[616] = { 4, 612, 614, 615 };
	pfverts[617] = { 613, 635, 614, 615 };
	pfverts[618] = { 129, 623, 1347, 622 };
	pfverts[619] = { 626, 32, 3, 772 };
	pfverts[620] = { 239, 632, 52, 332 };
	pfverts[621] = { 736, 127, 622, 623 };
	pfverts[622] = { 618, 621, 626, 627 };
	pfverts[623] = { 53, 636, 621, 618 };
	pfverts[624] = { 1354, 602, 627, 128 };
	pfverts[625] = { 638, 612, 4, 125 };
	pfverts[626] = { 619, 601, 622, 1332 };
	pfverts[627] = { 1332, 22, 622, 624 };
	pfverts[628] = { 614, -1, 129, 633 };
	pfverts[629] = { 640, 129, 22, 4 };
	pfverts[630] = { 127, 611, 634, 605 };
	pfverts[631] = { 126, 613, 635, 322 };
	pfverts[632] = { 735, 13, 53, 620 };
	pfverts[633] = { 51, 628, 635, 636 };
	pfverts[634] = { 630, 609, 635, 636 };
	pfverts[635] = { 617, 631, 633, 634 };
	pfverts[636] = { 600, 623, 633, 634 };
	pfverts[637] = { 761, 33, 1350, 1349 };
	pfverts[638] = { 625, -1, 640, 130 };
	pfverts[639] = { 125, 602, -1, 603 };
	pfverts[640] = { 638, -1, 641, 629 };
	pfverts[641] = { -1, 21, 640, -1 };
	pfverts[642] = { 156, -1, 130, 82 };
	pfverts[643] = { 967, 1183, 697, 1185 };
	pfverts[644] = { 246, 669, 647, 300 };
	pfverts[645] = { 135, 648, 671, 673 };
	pfverts[646] = { 699, 131, 667, 649 };
	pfverts[647] = { 644, 303, 648, 649 };
	pfverts[648] = { 194, 647, 667, 645 };
	pfverts[649] = { 1185, 62, 646, 647 };
	pfverts[650] = { 1265, 1329, 724, 660 };
	pfverts[651] = { 150, 1265, 1310, 1324 };
	pfverts[652] = { 1276, 1266, 1265, 1269 };
	pfverts[653] = { 193, 681, 658, 1015 };
	pfverts[654] = { 789, 1314, 659, 728 };
	pfverts[655] = { 790, 729, 657, 1044 };
	pfverts[656] = { 1152, 132, 657, 658 };
	pfverts[657] = { 731, 655, 656, 659 };
	pfverts[658] = { 656, 683, 659, 653 };
	pfverts[659] = { 654, 657, 658, 791 };
	pfverts[660] = { 1269, 650, 1311, 134 };
	pfverts[661] = { 740, 136, 1320, 83 };
	pfverts[662] = { 717, 233, 666, 680 };
	pfverts[663] = { 237, 133, 1316, 666 };
	pfverts[664] = { 1321, 1282, 1315, 666 };
	pfverts[665] = { 134, 1271, 1319, 30 };
	pfverts[666] = { 662, 684, 663, 664 };
	pfverts[667] = { 696, 646, 671, 648 };
	pfverts[668] = { 675, 1318, 1323, 672 };
	pfverts[669] = { 644, 797, 673, 681 };
	pfverts[670] = { 672, 1283, 795, 249 };
	pfverts[671] = { 674, 667, 672, 645 };
	pfverts[672] = { 670, 671, 668, 678 };
	pfverts[673] = { 669, 685, 795, 645 };
	pfverts[674] = { 671, 133, 1318, 138 };
	pfverts[675] = { 668, 1283, 1271, 1321 };
	pfverts[676] = { 677, 713, 1361, 710 };
	pfverts[677] = { 1320, 136, 30, 676 };
	pfverts[678] = { 249, 672, 6, 131 };
	pfverts[679] = { 1284, 1317, 5, 1283 };
	pfverts[680] = { 132, 662, 684, 1170 };
	pfverts[681] = { 653, 1169, 685, 669 };
	pfverts[682] = { 1171, 135, 684, 685 };
	pfverts[683] = { 1282, 658, 684, 685 };
	pfverts[684] = { 666, 680, 682, 683 };
	pfverts[685] = { 681, 673, 682, 683 };
	pfverts[686] = { 767, 692, 313, 40 };
	pfverts[687] = { -1, -1, -1, 825 }; // tetra spanning corner
	pfverts[688] = { -1, 698, 689, 823 };
	pfverts[689] = { 832, 701, 688, 825 };
	pfverts[690] = { 713, -1, 38, 136 };
	pfverts[691] = { 139, 234, 694, 1177 };
	pfverts[692] = { 686, 693, 1333, 1180 };
	pfverts[693] = { 38, 712, 694, 692 };
	pfverts[694] = { 235, 691, 1333, 693 };
	pfverts[695] = { 824, 1358, 699, 131 };
	pfverts[696] = { 667, 714, 700, 1182 };
	pfverts[697] = { 66, 305, 643, 35 };
	pfverts[698] = { 704, 699, 701, 688 };
	pfverts[699] = { 695, 646, 704, 698 };
	pfverts[700] = { 696, 1352, 715, 704 };
	pfverts[701] = { 1185, 698, 689, 248 };
	pfverts[702] = { 1177, 34, -1, 705 };
	pfverts[703] = { 1187, -1, -1, 711 };
	pfverts[704] = { 699, 698, 700, 1188 };
	pfverts[705] = { 714, 139, 1187, 702 };
	pfverts[706] = { 1187, -1, -1, 715 };
	pfverts[707] = { -1, 1325, 136, 11 };
	pfverts[708] = { 784, -1, 1360, 10 };
	pfverts[709] = { 1270, 1361, 1359, 710 };
	pfverts[710] = { 1318, 676, 709, 138 };
	pfverts[711] = { 712, -1, 703, 1180 };
	pfverts[712] = { 693, 139, 711, 713 };
	pfverts[713] = { 690, 676, 712, -1 };
	pfverts[714] = { 138, 705, 715, 696 };
	pfverts[715] = { 714, 700, 706, 1188 };
	pfverts[716] = { 146, 775, 720, 926 };
	pfverts[717] = { 662, 1315, 741, 745 };
	pfverts[718] = { 72, 744, 741, 146 };
	pfverts[719] = { 266, 769, 1181, 1141 };
	pfverts[720] = { 723, 73, -1, 716 };
	pfverts[721] = { 71, 43, 741, 745 };
	pfverts[722] = { 233, 746, 866, 741 };
	pfverts[723] = { 1181, 720, 1343, 746 };
	pfverts[724] = { 154, 787, 650, 1312 };
	pfverts[725] = { 787, 1275, 1313, 1309 };
	pfverts[726] = { 1269, 141, 1267, 1326 };
	pfverts[727] = { 1328, 1350, 1326, 749 };
	pfverts[728] = { 654, 752, 732, 605 };
	pfverts[729] = { 655, 606, 733, 919 };
	pfverts[730] = { 608, 142, 732, 733 };
	pfverts[731] = { 754, 657, 732, 733 };
	pfverts[732] = { 756, 728, 730, 731 };
	pfverts[733] = { 729, 923, 730, 731 };
	pfverts[734] = { 752, 145, 737, 127 };
	pfverts[735] = { 632, 239, 736, 753 };
	pfverts[736] = { 621, 735, 739, 737 };
	pfverts[737] = { 734, 758, 736, 738 };
	pfverts[738] = { 751, 737, 739, 743 };
	pfverts[739] = { 3, 736, 240, 738 };
	pfverts[740] = { 140, 748, 744, 661 };
	pfverts[741] = { 717, 718, 722, 721 };
	pfverts[742] = { 1314, 752, 1341, 1315 };
	pfverts[743] = { 738, 762, 143, 937 };
	pfverts[744] = { 750, 740, 144, 718 };
	pfverts[745] = { 165, 757, 717, 721 };
	pfverts[746] = { 238, 146, 723, 722 };
	pfverts[747] = { 771, 773, -1, 770 };
	pfverts[748] = { 776, 740, 1325, 749 };
	pfverts[749] = { 750, 727, 748, 770 };
	pfverts[750] = { 143, 749, 751, 744 };
	pfverts[751] = { 738, 1328, 145, 750 };
	pfverts[752] = { 742, 734, 756, 728 };
	pfverts[753] = { 142, 735, 758, 936 };
	pfverts[754] = { 731, 757, 165, 922 };
	pfverts[755] = { 932, 937, 940, 166 };
	pfverts[756] = { 752, 732, 758, 757 };
	pfverts[757] = { 745, 754, 756, 1334 };
	pfverts[758] = { 737, 756, 753, 937 };
	pfverts[759] = { 760, -1, 763, 67 };
	pfverts[760] = { -1, -1, -1, 759 }; // tetra spanning corner
	pfverts[761] = { 637, 67, 763, 764 };
	pfverts[762] = { 743, 148, 244, 240 };
	pfverts[763] = { 771, 759, 148, 761 };
	pfverts[764] = { 240, 761, 242, 33 };
	pfverts[765] = { -1, 14, 242, 245 };
	pfverts[766] = { -1, -1, 774, 74 };
	pfverts[767] = { 686, 768, -1, 1178 };
	pfverts[768] = { 767, -1, 73, 40 };
	pfverts[769] = { -1, 137, 1178, 719 };
	pfverts[770] = { 749, 147, 747, 776 };
	pfverts[771] = { 1350, 763, -1, 747 };
	pfverts[772] = { 601, 619, 33, 1349 };
	pfverts[773] = { 147, -1, 148, 747 };
	pfverts[774] = { 39, 766, 73, 775 };
	pfverts[775] = { 140, 74, 774, 716 };
	pfverts[776] = { 770, -1, 748, 74 };
	pfverts[777] = { 230, 804, 781, 282 };
	pfverts[778] = { 283, 251, 833, 158 };
	pfverts[779] = { 284, 822, -1, -1 };
	pfverts[780] = { 822, 149, 782, 783 };
	pfverts[781] = { 777, 285, 782, 783 };
	pfverts[782] = { 808, 833, 780, 781 };
	pfverts[783] = { 284, 231, 780, 781 };
	pfverts[784] = { 708, 26, -1, 1310 };
	pfverts[785] = { 810, 786, 48, 49 };
	pfverts[786] = { -1, 812, 785, 26 };
	pfverts[787] = { 1268, 724, 48, 725 };
	pfverts[788] = { 175, 816, 792, 126 };
	pfverts[789] = { 605, 814, 793, 654 };
	pfverts[790] = { 606, 655, 794, 1017 };
	pfverts[791] = { 659, 815, 793, 794 };
	pfverts[792] = { 788, 607, 793, 794 };
	pfverts[793] = { 818, 789, 791, 792 };
	pfverts[794] = { 790, 1021, 791, 792 };
	pfverts[795] = { 673, 799, 801, 670 };
	pfverts[796] = { 801, 250, 152, 802 };
	pfverts[797] = { 669, 246, 801, 1015 };
	pfverts[798] = { 249, 152, 252, 800 };
	pfverts[799] = { 1282, 5, 795, 1283 };
	pfverts[800] = { 835, 811, 812, 798 };
	pfverts[801] = { 797, 817, 796, 795 };
	pfverts[802] = { 820, 796, 251, 176 };
	pfverts[803] = { 149, 837, 806, 612 };
	pfverts[804] = { 777, 613, 808, 816 };
	pfverts[805] = { 813, 153, 807, 806 };
	pfverts[806] = { 803, 805, 808, 615 };
	pfverts[807] = { 809, 820, 808, 805 };
	pfverts[808] = { 782, 807, 804, 806 };
	pfverts[809] = { 807, 251, 811, 833 };
	pfverts[810] = { 812, 1317, 813, 785 };
	pfverts[811] = { 834, 809, 813, 800 };
	pfverts[812] = { 786, 810, 249, 800 };
	pfverts[813] = { 811, 805, 837, 810 };
	pfverts[814] = { 153, 5, 817, 789 };
	pfverts[815] = { 818, 791, 1015, 151 };
	pfverts[816] = { 788, 1030, 819, 804 };
	pfverts[817] = { 801, 814, 818, 820 };
	pfverts[818] = { 793, 815, 817, 819 };
	pfverts[819] = { 818, 1032, 820, 816 };
	pfverts[820] = { 802, 817, 819, 807 };
	pfverts[821] = { -1, -1, 828, -1 }; // tetra spanning corner
	pfverts[822] = { 779, 155, 827, 780 };
	pfverts[823] = { 688, 832, -1, 824 };
	pfverts[824] = { 695, 823, -1, 157 };
	pfverts[825] = { -1, 248, 689, 687 };
	pfverts[826] = { 247, 827, -1, 836 };
	pfverts[827] = { -1, 828, 826, 822 };
	pfverts[828] = { 821, 12, 827, 155 };
	pfverts[829] = { 158, -1, 1014, 838 };
	pfverts[830] = { 1185, 831, 62, 272 };
	pfverts[831] = { 832, -1, 1014, 830 };
	pfverts[832] = { 823, 689, -1, 831 };
	pfverts[833] = { 782, 809, 247, 778 };
	pfverts[834] = { 247, 811, 836, 838 };
	pfverts[835] = { 800, 157, -1, 838 };
	pfverts[836] = { 834, 826, -1, 47 };
	pfverts[837] = { 803, 47, 813, 50 };
	pfverts[838] = { 834, 835, -1, 829 };
	pfverts[839] = { 260, 854, 840, 1126 };
	pfverts[840] = { 839, -1, 853, 841 };
	pfverts[841] = { -1, 265, -1, 840 };
	pfverts[842] = { 89, 859, 846, 163 };
	pfverts[843] = { 847, 7, 892, 934 };
	pfverts[844] = { 893, 919, 848, 327 };
	pfverts[845] = { 921, 159, 847, 848 };
	pfverts[846] = { 842, 894, 847, 848 };
	pfverts[847] = { 861, 843, 845, 846 };
	pfverts[848] = { 844, 331, 845, 846 };
	pfverts[849] = { 930, 253, 852, 858 };
	pfverts[850] = { 254, 255, -1, 852 };
	pfverts[851] = { 857, 920, -1, 931 };
	pfverts[852] = { 849, 862, 850, 931 };
	pfverts[853] = { -1, 840, -1, 855 };
	pfverts[854] = { 839, 899, 856, 859 };
	pfverts[855] = { 90, 853, -1, 856 };
	pfverts[856] = { 854, 863, 898, 855 };
	pfverts[857] = { 900, 851, 861, 934 };
	pfverts[858] = { 159, 849, 862, 334 };
	pfverts[859] = { 842, 333, 863, 854 };
	pfverts[860] = { -1, 160, 90, 255 };
	pfverts[861] = { 857, 847, 862, 863 };
	pfverts[862] = { 852, 858, 160, 861 };
	pfverts[863] = { 859, 856, 335, 861 };
	pfverts[864] = { 309, 878, 928, 41 };
	pfverts[865] = { 162, -1, 16, 878 };
	pfverts[866] = { 722, 41, 928, 43 };
	pfverts[867] = { 43, 880, 882, 879 };
	pfverts[868] = { 920, 918, 7, 934 };
	pfverts[869] = { 919, 893, 873, 1198 };
	pfverts[870] = { 1196, 884, 871, 872 };
	pfverts[871] = { 892, 895, 870, 873 };
	pfverts[872] = { 870, 935, 873, 1202 };
	pfverts[873] = { 25, 871, 872, 869 };
	pfverts[874] = { 1189, 161, 257, 875 };
	pfverts[875] = { 222, 258, 874, 877 };
	pfverts[876] = { -1, 907, 927, 924 };
	pfverts[877] = { 887, 875, 257, -1 };
	pfverts[878] = { 307, -1, 865, 864 };
	pfverts[879] = { 918, 886, 867, 935 };
	pfverts[880] = { 1218, 928, 16, 867 };
	pfverts[881] = { 911, 892, 927, 907 };
	pfverts[882] = { 165, 867, 198, 935 };
	pfverts[883] = { 1220, 162, 885, 16 };
	pfverts[884] = { 927, 870, 885, 886 };
	pfverts[885] = { 883, 1219, 884, 161 };
	pfverts[886] = { 884, 879, 16, 924 };
	pfverts[887] = { 877, -1, 1083, 889 };
	pfverts[888] = { 366, 490, 486, 891 };
	pfverts[889] = { 220, -1, 887, 1342 };
	pfverts[890] = { 487, 20, 221, 110 };
	pfverts[891] = { 350, 888, 506, 909 };
	pfverts[892] = { 896, 843, 871, 881 };
	pfverts[893] = { 869, 844, 897, 479 };
	pfverts[894] = { 846, 163, 896, 897 };
	pfverts[895] = { 109, 871, 896, 897 };
	pfverts[896] = { 911, 892, 894, 895 };
	pfverts[897] = { 893, 483, 894, 895 };
	pfverts[898] = { 856, 901, -1, 900 };
	pfverts[899] = { 854, 260, 901, 908 };
	pfverts[900] = { 857, 898, -1, 907 };
	pfverts[901] = { 899, 912, 264, 898 };
	pfverts[902] = { 345, 222, 220, 223 };
	pfverts[903] = { 490, 258, 906, 909 };
	pfverts[904] = { 258, 927, -1, 906 };
	pfverts[905] = { 906, -1, 20, 164 };
	pfverts[906] = { 903, 913, 904, 905 };
	pfverts[907] = { 900, 934, 881, 876 };
	pfverts[908] = { 163, 899, 912, 507 };
	pfverts[909] = { 1, 891, 109, 903 };
	pfverts[910] = { 509, 164, 912, 913 };
	pfverts[911] = { 881, 896, 912, 913 };
	pfverts[912] = { 901, 908, 910, 911 };
	pfverts[913] = { 109, 906, 910, 911 };
	pfverts[914] = { 253, 930, 915, 239 };
	pfverts[915] = { 914, 243, 929, 916 };
	pfverts[916] = { 68, 256, -1, 915 };
	pfverts[917] = { 142, 936, 921, 159 };
	pfverts[918] = { 924, 868, 879, 925 };
	pfverts[919] = { 844, 869, 923, 729 };
	pfverts[920] = { 851, 939, 868, 925 };
	pfverts[921] = { 917, 845, 938, 923 };
	pfverts[922] = { 923, 754, 25, 938 };
	pfverts[923] = { 919, 733, 922, 921 };
	pfverts[924] = { 876, -1, 918, 886 };
	pfverts[925] = { -1, 920, 918, 42 };
	pfverts[926] = { 75, 72, -1, 716 };
	pfverts[927] = { 904, 881, 884, 876 };
	pfverts[928] = { 306, 864, 880, 866 };
	pfverts[929] = { 915, 933, -1, 932 };
	pfverts[930] = { 914, 849, 933, 936 };
	pfverts[931] = { 852, 851, -1, 933 };
	pfverts[932] = { 244, 929, -1, 755 };
	pfverts[933] = { 930, 940, 931, 929 };
	pfverts[934] = { 857, 868, 843, 907 };
	pfverts[935] = { 25, 879, 882, 872 };
	pfverts[936] = { 917, 753, 940, 930 };
	pfverts[937] = { 743, 758, 755, 1334 };
	pfverts[938] = { 921, 922, 7, 939 };
	pfverts[939] = { 940, 938, 920, 166 };
	pfverts[940] = { 936, 933, 755, 939 };
	pfverts[941] = { 541, 961, 944, 294 };
	pfverts[942] = { 956, 944, -1, 958 };
	pfverts[943] = { 297, 1331, 54, 966 };
	pfverts[944] = { 941, 297, 942, 966 };
	pfverts[945] = { 1112, 23, 961, 999 };
	pfverts[946] = { 977, 547, 59, 973 };
	pfverts[947] = { 993, 1017, 951, 1044 };
	pfverts[948] = { 167, 8, 949, 950 };
	pfverts[949] = { 59, 1019, 948, 951 };
	pfverts[950] = { 948, 957, 951, 1096 };
	pfverts[951] = { 24, 949, 950, 947 };
	pfverts[952] = { -1, 955, 998, 962 };
	pfverts[953] = { 1026, 268, 959, 960 };
	pfverts[954] = { -1, 1029, 1028, 962 };
	pfverts[955] = { 271, 964, 952, 959 };
	pfverts[956] = { 169, 295, -1, 942 };
	pfverts[957] = { 950, 24, 23, 999 };
	pfverts[958] = { 1007, 964, 998, 942 };
	pfverts[959] = { 953, 8, 955, 962 };
	pfverts[960] = { 167, 953, 8, 1113 };
	pfverts[961] = { 941, 945, 277, 1007 };
	pfverts[962] = { 952, 1016, 954, 959 };
	pfverts[963] = { 999, 1114, 8, 964 };
	pfverts[964] = { 958, 169, 963, 955 };
	pfverts[965] = { 1120, 1119, 1345, -1 };
	pfverts[966] = { 944, -1, 278, 943 };
	pfverts[967] = { 643, -1, 272, 66 };
	pfverts[968] = { 282, 982, 970, 320 };
	pfverts[969] = { -1, 970, -1, 983 };
	pfverts[970] = { 968, -1, 969, 971 };
	pfverts[971] = { -1, 287, -1, 970 };
	pfverts[972] = { 204, 986, 976, 175 };
	pfverts[973] = { 1016, 994, 61, 946 };
	pfverts[974] = { 1017, 993, 978, 1234 };
	pfverts[975] = { 995, 170, 977, 978 };
	pfverts[976] = { 972, 1018, 977, 978 };
	pfverts[977] = { 984, 946, 975, 976 };
	pfverts[978] = { 974, 1238, 975, 976 };
	pfverts[979] = { 1002, 273, 1003, 985 };
	pfverts[980] = { 1006, 1005, -1, 57 };
	pfverts[981] = { 57, 1243, 987, 984 };
	pfverts[982] = { 968, 1023, 983, 986 };
	pfverts[983] = { 982, 989, 171, 969 };
	pfverts[984] = { 61, 981, 276, 977 };
	pfverts[985] = { 170, 979, 1243, 299 };
	pfverts[986] = { 972, 1239, 989, 982 };
	pfverts[987] = { -1, 988, 276, 981 };
	pfverts[988] = { 172, -1, 56, 987 };
	pfverts[989] = { 986, 983, 172, 276 };
	pfverts[990] = { 273, 1002, 991, 224 };
	pfverts[991] = { 1004, 990, 1001, 225 };
	pfverts[992] = { 117, 1008, 995, 170 };
	pfverts[993] = { 974, 947, 997, 549 };
	pfverts[994] = { 547, 973, 168, 1000 };
	pfverts[995] = { 992, 975, 1010, 997 };
	pfverts[996] = { 997, 577, 24, 1010 };
	pfverts[997] = { 993, 553, 996, 995 };
	pfverts[998] = { 958, -1, 952, 1000 };
	pfverts[999] = { 963, 957, 945, 1007 };
	pfverts[1000] = { 1007, 994, 60, 998 };
	pfverts[1001] = { 568, 991, 228, 229 };
	pfverts[1002] = { 990, 979, 1005, 1008 };
	pfverts[1003] = { 979, 57, 1243, 58 };
	pfverts[1004] = { 118, 991, -1, 1005 };
	pfverts[1005] = { 1002, 578, 980, 1004 };
	pfverts[1006] = { -1, 980, 61, 168 };
	pfverts[1007] = { 958, 999, 961, 1000 };
	pfverts[1008] = { 992, 576, 1011, 1002 };
	pfverts[1009] = { 174, 561, 1012, 578 };
	pfverts[1010] = { 995, 996, 547, 1011 };
	pfverts[1011] = { 1008, 1012, 1010, 578 };
	pfverts[1012] = { 581, 580, 1011, 1009 };
	pfverts[1013] = { 268, 1026, 64, 246 };
	pfverts[1014] = { 157, 829, 63, 831 };
	pfverts[1015] = { 653, 797, 815, 1031 };
	pfverts[1016] = { 973, 59, 1029, 962 };
	pfverts[1017] = { 947, 974, 1021, 790 };
	pfverts[1018] = { 976, 175, 1020, 1021 };
	pfverts[1019] = { 151, 949, 1020, 1021 };
	pfverts[1020] = { 1033, 59, 1018, 1019 };
	pfverts[1021] = { 1017, 794, 1018, 1019 };
	pfverts[1022] = { 1024, 9, -1, 171 };
	pfverts[1023] = { 982, 282, 1024, 1030 };
	pfverts[1024] = { 1023, 1034, 286, 1022 };
	pfverts[1025] = { 65, 64, -1, 272 };
	pfverts[1026] = { 1013, 953, 1028, 1031 };
	pfverts[1027] = { 1028, -1, 64, 176 };
	pfverts[1028] = { 1026, 1035, 954, 1027 };
	pfverts[1029] = { 9, 1033, 954, 1016 };
	pfverts[1030] = { 175, 1023, 1034, 816 };
	pfverts[1031] = { 167, 1015, 151, 1026 };
	pfverts[1032] = { 819, 176, 1034, 1035 };
	pfverts[1033] = { 1029, 1020, 1034, 1035 };
	pfverts[1034] = { 1024, 1030, 1032, 1033 };
	pfverts[1035] = { 151, 1028, 1032, 1033 };
	pfverts[1036] = { 294, 1057, 1040, 207 };
	pfverts[1037] = { 1039, 1062, 1059, 1060 };
	pfverts[1038] = { 388, 401, 1080, 1074 };
	pfverts[1039] = { 1037, 1040, 589, 1331 };
	pfverts[1040] = { 1036, 211, 1059, 1039 };
	pfverts[1041] = { 208, 182, 1059, 1072 };
	pfverts[1042] = { 92, 1067, 1046, 185 };
	pfverts[1043] = { 1091, 1065, 1047, 1197 };
	pfverts[1044] = { 947, 655, 1148, 1198 };
	pfverts[1045] = { 1200, 178, 1047, 1048 };
	pfverts[1046] = { 1042, 1093, 1047, 1048 };
	pfverts[1047] = { 1069, 1043, 1045, 1046 };
	pfverts[1048] = { 1092, 353, 1045, 1046 };
	pfverts[1049] = { 1055, 200, 1211, 1214 };
	pfverts[1050] = { 1351, 1054, 1053, 1211 };
	pfverts[1051] = { 1210, 288, 1055, 1066 };
	pfverts[1052] = { 340, 292, 179, 182 };
	pfverts[1053] = { 181, 201, 1050, 200 };
	pfverts[1054] = { -1, 1063, 1050, 1213 };
	pfverts[1055] = { 1051, 1070, 179, 1049 };
	pfverts[1056] = { 183, 1072, 397, 408 };
	pfverts[1057] = { 1036, 1098, 1061, 1067 };
	pfverts[1058] = { 1100, 180, 1060, 1061 };
	pfverts[1059] = { 1040, 1041, 1061, 1037 };
	pfverts[1060] = { 1064, 177, 1058, 1037 };
	pfverts[1061] = { 1057, 1071, 1058, 1059 };
	pfverts[1062] = { 1064, 1037, 1072, 408 };
	pfverts[1063] = { 181, 1054, -1, 1064 };
	pfverts[1064] = { 1063, 1062, -1, 1060 };
	pfverts[1065] = { 180, 200, 1069, 1043 };
	pfverts[1066] = { 178, 1051, 1070, 380 };
	pfverts[1067] = { 1042, 379, 1071, 1057 };
	pfverts[1068] = { 381, 182, 1070, 1071 };
	pfverts[1069] = { 1065, 1047, 1070, 1071 };
	pfverts[1070] = { 1055, 1066, 1068, 1069 };
	pfverts[1071] = { 1067, 1061, 1068, 1069 };
	pfverts[1072] = { 292, 1041, 1056, 1062 };
	pfverts[1073] = { 1074, 15, 589, 1080 };
	pfverts[1074] = { 1353, 1038, -1, 1073 };
	pfverts[1075] = { 392, 183, 1081, 292 };
	pfverts[1076] = { 312, 342, -1, 1077 };
	pfverts[1077] = { 393, 1076, -1, 1081 };
	pfverts[1078] = { 397, 17, 413, 407 };
	pfverts[1079] = { 1356, 407, 413, -1 };
	pfverts[1080] = { 1038, 1073, 55, 587 };
	pfverts[1081] = { 290, 1077, -1, 1075 };
	pfverts[1082] = { 385, 411, 514, 386 };
	pfverts[1083] = { 1195, 70, 887, 1191 };
	pfverts[1084] = { 300, 1104, 1087, 268 };
	pfverts[1085] = { 269, 1103, 1088, 190 };
	pfverts[1086] = { 1120, 86, 1088, 1089 };
	pfverts[1087] = { 1084, 270, 1088, 1089 };
	pfverts[1088] = { 1106, 1085, 1086, 1087 };
	pfverts[1089] = { 1345, 305, 1086, 1087 };
	pfverts[1090] = { 167, 1113, 1094, 193 };
	pfverts[1091] = { 1147, 1111, 1095, 1043 };
	pfverts[1092] = { 1148, 1048, 1198, 1202 };
	pfverts[1093] = { 1046, 185, 1095, 1096 };
	pfverts[1094] = { 1090, 1149, 1095, 1096 };
	pfverts[1095] = { 1115, 1091, 1093, 1094 };
	pfverts[1096] = { 1148, 950, 1093, 1094 };
	pfverts[1097] = { 1111, 188, 1100, 180 };
	pfverts[1098] = { 1057, 294, 1102, 1112 };
	pfverts[1099] = { 298, 186, 1101, 1102 };
	pfverts[1100] = { 1097, 1058, 1101, 1102 };
	pfverts[1101] = { 1110, 177, 1099, 1100 };
	pfverts[1102] = { 1098, 1116, 1099, 1100 };
	pfverts[1103] = { 189, 1109, 1106, 1085 };
	pfverts[1104] = { 1084, 1154, 1108, 1113 };
	pfverts[1105] = { 301, 196, 187, 1153 };
	pfverts[1106] = { 1103, 1088, 1157, 1108 };
	pfverts[1107] = { -1, 1124, 1157, 1166 };
	pfverts[1108] = { 1104, 1117, 1156, 1106 };
	pfverts[1109] = { 1103, 1110, 1123, 186 };
	pfverts[1110] = { 1101, 1109, -1, 1124 };
	pfverts[1111] = { 1156, 1097, 1115, 1091 };
	pfverts[1112] = { 185, 1098, 1116, 945 };
	pfverts[1113] = { 1090, 960, 1117, 1104 };
	pfverts[1114] = { 963, 189, 1116, 1117 };
	pfverts[1115] = { 1111, 1095, 1116, 1117 };
	pfverts[1116] = { 1102, 1112, 1114, 1115 };
	pfverts[1117] = { 1113, 1108, 1114, 1115 };
	pfverts[1118] = { 190, 1120, -1, -1 };
	pfverts[1119] = { -1, -1, 66, 965 };
	pfverts[1120] = { -1, 1086, 1118, 965 };
	pfverts[1121] = { -1, 1125, 1122, 295 };
	pfverts[1122] = { -1, 1121, -1, 296 };
	pfverts[1123] = { -1, 190, 1109, 1125 };
	pfverts[1124] = { 188, 1110, -1, 1107 };
	pfverts[1125] = { -1, 1123, -1, 1121 };
	pfverts[1126] = { 839, 333, -1, 1134 };
	pfverts[1127] = { 102, 1136, 1131, 89 };
	pfverts[1128] = { 1346, 318, -1, 1237 };
	pfverts[1129] = { 327, 1234, 1133, 425 };
	pfverts[1130] = { 1236, 191, 1132, 1133 };
	pfverts[1131] = { 1127, 328, 1132, 1133 };
	pfverts[1132] = { 1138, 1346, 1130, 1131 };
	pfverts[1133] = { 1129, 429, 1130, 1131 };
	pfverts[1134] = { 215, 1126, -1, 1136 };
	pfverts[1135] = { 191, 1232, 1138, 451 };
	pfverts[1136] = { 1127, 450, 1137, 1134 };
	pfverts[1137] = { 452, 1136, -1, 1138 };
	pfverts[1138] = { 1135, 1137, -1, 1132 };
	pfverts[1139] = { 306, 1160, 1143, 233 };
	pfverts[1140] = { 234, 1159, 1144, 1177 };
	pfverts[1141] = { 1179, 1145, 719, 1330 };
	pfverts[1142] = { 1175, 87, 1144, 1145 };
	pfverts[1143] = { 1139, 236, 1144, 1145 };
	pfverts[1144] = { 1162, 1140, 1142, 1143 };
	pfverts[1145] = { 1141, 311, 1142, 1143 };
	pfverts[1146] = { 132, 1170, 1150, 198 };
	pfverts[1147] = { 1197, 1168, 1151, 1091 };
	pfverts[1148] = { 1152, 1096, 1044, 1092 };
	pfverts[1149] = { 1094, 193, 1151, 1152 };
	pfverts[1150] = { 1146, 1199, 1151, 1152 };
	pfverts[1151] = { 1172, 1147, 1149, 1150 };
	pfverts[1152] = { 1148, 656, 1149, 1150 };
	pfverts[1153] = { 1158, 1105, 1156, 1168 };
	pfverts[1154] = { 1104, 300, 1158, 1169 };
	pfverts[1155] = { 197, 194, 1158, 301 };
	pfverts[1156] = { 1108, 187, 1111, 1153 };
	pfverts[1157] = { 86, 1106, 187, 1107 };
	pfverts[1158] = { 1154, 1173, 1155, 1153 };
	pfverts[1159] = { 197, 1165, 1162, 1140 };
	pfverts[1160] = { 1139, 1204, 1164, 1170 };
	pfverts[1161] = { 1206, 195, 1163, 1164 };
	pfverts[1162] = { 1159, 1144, 1163, 1164 };
	pfverts[1163] = { 1167, 192, 1161, 1162 };
	pfverts[1164] = { 1160, 1174, 1161, 1162 };
	pfverts[1165] = { 1159, 301, 1167, 34 };
	pfverts[1166] = { 188, 1107, 184, 187 };
	pfverts[1167] = { 196, 1165, -1, 1163 };
	pfverts[1168] = { 195, 1153, 1172, 1147 };
	pfverts[1169] = { 193, 1154, 1173, 681 };
	pfverts[1170] = { 1146, 680, 1174, 1160 };
	pfverts[1171] = { 682, 197, 1173, 1174 };
	pfverts[1172] = { 1168, 1151, 1173, 1174 };
	pfverts[1173] = { 1158, 1169, 1171, 1172 };
	pfverts[1174] = { 1170, 1164, 1171, 1172 };
	pfverts[1175] = { 1142, 1177, -1, 1179 };
	pfverts[1176] = { 248, -1, 1184, 1352 };
	pfverts[1177] = { 702, 691, 1175, 1140 };
	pfverts[1178] = { 767, 1333, 1181, 769 };
	pfverts[1179] = { 1141, 1175, -1, 137 };
	pfverts[1180] = { -1, 37, 692, 711 };
	pfverts[1181] = { 723, 1178, 238, 719 };
	pfverts[1182] = { 304, 696, 1186, 34 };
	pfverts[1183] = { -1, 1184, 643, 248 };
	pfverts[1184] = { 1183, -1, 35, 1176 };
	pfverts[1185] = { 830, 701, 643, 649 };
	pfverts[1186] = { 36, 1352, -1, 1182 };
	pfverts[1187] = { 703, 706, 705, -1 };
	pfverts[1188] = { 1359, 704, 715, -1 };
	pfverts[1189] = { 288, 1210, 1193, 874 };
	pfverts[1190] = { 259, 1209, 1194, 203 };
	pfverts[1191] = { 1225, 393, 1083, -1 };
	pfverts[1192] = { 1224, 84, 1194, 1195 };
	pfverts[1193] = { 1189, 257, 1194, 1195 };
	pfverts[1194] = { 1212, 1190, 1192, 1193 };
	pfverts[1195] = { 1083, 293, 1192, 1193 };
	pfverts[1196] = { 1200, 870, 1219, 1 };
	pfverts[1197] = { 1043, 1217, 1201, 1147 };
	pfverts[1198] = { 352, 1044, 1092, 869 };
	pfverts[1199] = { 1150, 198, 1201, 1202 };
	pfverts[1200] = { 1196, 1045, 1201, 1202 };
	pfverts[1201] = { 1221, 1197, 1199, 1200 };
	pfverts[1202] = { 1092, 872, 1199, 1200 };
	pfverts[1203] = { 1217, 201, 1206, 195 };
	pfverts[1204] = { 1160, 306, 1208, 1218 };
	pfverts[1205] = { 310, 199, 1207, 1208 };
	pfverts[1206] = { 1203, 1161, 1207, 1208 };
	pfverts[1207] = { 1216, 192, 1205, 1206 };
	pfverts[1208] = { 1204, 1222, 1205, 1206 };
	pfverts[1209] = { 202, 1215, 1212, 1190 };
	pfverts[1210] = { 1189, 1051, 1214, 2 };
	pfverts[1211] = { 1212, 1213, 1050, 1049 };
	pfverts[1212] = { 1209, 1194, 1211, 1214 };
	pfverts[1213] = { 84, 289, 1054, 1211 };
	pfverts[1214] = { 1210, 1223, 1049, 1212 };
	pfverts[1215] = { 1209, 1216, 1229, 199 };
	pfverts[1216] = { 1351, 1207, -1, 1215 };
	pfverts[1217] = { 200, 1203, 1221, 1197 };
	pfverts[1218] = { 198, 1204, 1222, 880 };
	pfverts[1219] = { 1223, 1196, 885, 2 };
	pfverts[1220] = { 883, 202, 1222, 1223 };
	pfverts[1221] = { 1217, 1201, 1222, 1223 };
	pfverts[1222] = { 1208, 1218, 1220, 1221 };
	pfverts[1223] = { 1219, 1214, 1220, 1221 };
	pfverts[1224] = { -1, 1192, 203, 1225 };
	pfverts[1225] = { -1, 1191, 1224, -1 };
	pfverts[1226] = { 308, -1, 1231, 1228 };
	pfverts[1227] = { 307, 1229, -1, 1231 };
	pfverts[1228] = { 1330, -1, 1226, -1 };
	pfverts[1229] = { -1, 1230, 1215, 1227 };
	pfverts[1230] = { 1229, 203, -1, -1 };
	pfverts[1231] = { 1226, 1227, -1, -1 };
	pfverts[1232] = { 1135, 299, -1, 1245 };
	pfverts[1233] = { 170, 299, 1236, 191 };
	pfverts[1234] = { 1129, 315, 1238, 974 };
	pfverts[1235] = { 317, 204, 1237, 1238 };
	pfverts[1236] = { 1233, 1130, 1237, 1238 };
	pfverts[1237] = { 1128, 1235, 1241, 1236 };
	pfverts[1238] = { 1234, 978, 1235, 1236 };
	pfverts[1239] = { 204, 1240, 1244, 986 };
	pfverts[1240] = { 322, 320, -1, 1239 };
	pfverts[1241] = { 1244, -1, 1237, 1242 };
	pfverts[1242] = { -1, 56, 299, 1241 };
	pfverts[1243] = { 985, 56, 1003, 981 };
	pfverts[1244] = { -1, 1241, 1239, 172 };
	pfverts[1245] = { 416, 1232, -1, 273 };
	pfverts[1246] = { 403, -1, 29, 1299 };
	pfverts[1247] = { 1277, 1254, 1248, 1307 };
	pfverts[1248] = { 1258, 29, 1247, 347 };
	pfverts[1249] = { 420, 347, 349, 1251 };
	pfverts[1250] = { 1336, 349, 1274, 1253 };
	pfverts[1251] = { 346, 1258, 1249, 1253 };
	pfverts[1252] = { 29, 1274, 28, 1259 };
	pfverts[1253] = { 1337, 1251, 1250, 28 };
	pfverts[1254] = { 1290, 1247, 1258, 1260 };
	pfverts[1255] = { -1, 78, 1336, 1274 };
	pfverts[1256] = { 1301, 78, 1294, 1259 };
	pfverts[1257] = { -1, 214, 349, 347 };
	pfverts[1258] = { 1248, 1254, 28, 1251 };
	pfverts[1259] = { 1263, 1256, 1299, 1252 };
	pfverts[1260] = { 1254, 1263, 362, 27 };
	pfverts[1261] = { 1337, 1280, 1264, 1300 };
	pfverts[1262] = { 1281, 374, 362, 1264 };
	pfverts[1263] = { 1264, 28, 1259, 1260 };
	pfverts[1264] = { 1261, 1301, 1262, 1263 };
	pfverts[1265] = { 1268, 652, 651, 650 };
	pfverts[1266] = { -1, 1324, 652, 11 };
	pfverts[1267] = { 726, 1327, 1276, -1 };
	pfverts[1268] = { 787, 1275, 1310, 1265 };
	pfverts[1269] = { 660, 652, 11, 726 };
	pfverts[1270] = { 1360, 1322, 709, 1323 };
	pfverts[1271] = { 1284, 675, 665, 1273 };
	pfverts[1272] = { 1273, 1322, 30, 1324 };
	pfverts[1273] = { 150, 1323, 1271, 1272 };
	pfverts[1274] = { 214, 1255, 1250, 1252 };
	pfverts[1275] = { -1, 1276, 1268, 725 };
	pfverts[1276] = { -1, 1267, 652, 1275 };
	pfverts[1277] = { 403, 372, 1247, 1306 };
	pfverts[1278] = { 1279, 382, 1289, 564 };
	pfverts[1279] = { 374, 367, 378, 1278 };
	pfverts[1280] = { 1337, 1287, 1302, 1261 };
	pfverts[1281] = { 1302, 95, 1262, 378 };
	pfverts[1282] = { 799, 683, 1314, 664 };
	pfverts[1283] = { 675, 679, 799, 670 };
	pfverts[1284] = { 1329, 134, 1271, 679 };
	pfverts[1285] = { 372, 1355, 403, 376 };
	pfverts[1286] = { 544, 1308, 1290, 120 };
	pfverts[1287] = { 1305, 1340, 112, 1280 };
	pfverts[1288] = { 112, 421, 544, 105 };
	pfverts[1289] = { 351, 378, 1278, 119 };
	pfverts[1290] = { 1254, 346, 1286, 95 };
	pfverts[1291] = { 360, 1297, 1298, 1295 };
	pfverts[1292] = { 486, 488, 484, 367 };
	pfverts[1293] = { 1300, 484, 369, 365 };
	pfverts[1294] = { 412, 405, 1256, 404 };
	pfverts[1295] = { 93, 1291, 402, 371 };
	pfverts[1296] = { 357, 362, 1298, 1297 };
	pfverts[1297] = { 358, 27, 1296, 1291 };
	pfverts[1298] = { 564, 1296, 374, 1291 };
	pfverts[1299] = { 376, 404, 1246, 1259 };
	pfverts[1300] = { 1301, 373, 1261, 1293 };
	pfverts[1301] = { 1300, 377, 1264, 1256 };
	pfverts[1302] = { 112, 346, 1280, 1281 };
	pfverts[1303] = { 1357, 421, 476, 1305 };
	pfverts[1304] = { 476, 473, 1336, -1 };
	pfverts[1305] = { 1287, 473, 1303, 504 };
	pfverts[1306] = { 1277, 1307, -1, 570 };
	pfverts[1307] = { 1308, 546, 1247, 1306 };
	pfverts[1308] = { 1307, 1286, 545, 572 };
	pfverts[1309] = { 49, -1, 725, 603 };
	pfverts[1310] = { 651, 1338, 784, 1268 };
	pfverts[1311] = { 1312, 141, 660, 145 };
	pfverts[1312] = { 128, 1313, 724, 1311 };
	pfverts[1313] = { 1312, 725, 1327, 1354 };
	pfverts[1314] = { 654, 5, 1282, 742 };
	pfverts[1315] = { 742, 717, 144, 664 };
	pfverts[1316] = { 663, 1320, 1321, 1318 };
	pfverts[1317] = { 679, 153, 810, 154 };
	pfverts[1318] = { 668, 1316, 710, 674 };
	pfverts[1319] = { 144, 665, 1321, 1320 };
	pfverts[1320] = { 661, 677, 1319, 1316 };
	pfverts[1321] = { 664, 1319, 675, 1316 };
	pfverts[1322] = { 1270, 1361, 1272, 10 };
	pfverts[1323] = { 1273, 1270, 6, 668 };
	pfverts[1324] = { 651, 10, 1272, 1266 };
	pfverts[1325] = { -1, 707, 1326, 748 };
	pfverts[1326] = { 726, -1, 1325, 727 };
	pfverts[1327] = { 141, 1362, 1313, 1267 };
	pfverts[1328] = { 751, 3, 141, 727 };
	pfverts[1329] = { 154, 650, 150, 1284 };
	pfverts[1330] = { 1228, 311, 1141, 266 };
	pfverts[1331] = { 1039, 589, 209, 943 };
	pfverts[1332] = { 626, 32, 21, 627 };
	pfverts[1333] = { 1178, 694, 137, 692 };
	pfverts[1334] = { 71, 757, 937, 166 };
	pfverts[1335] = { 520, 220, 110, 223 };
	pfverts[1336] = { 1304, 1255, 348, 1250 };
	pfverts[1337] = { 1280, 1253, 348, 1261 };
	pfverts[1338] = { 150, 48, 26, 1310 };
	pfverts[1339] = { 119, 95, 357, 120 };
	pfverts[1340] = { 503, 484, 94, 1287 };
	pfverts[1341] = { 742, 134, 144, 145 };
	pfverts[1342] = { 220, 219, 889, 69 };
	pfverts[1343] = { 723, 41, -1, 266 };
	pfverts[1344] = { 107, 433, -1, 106 };
	pfverts[1345] = { 1089, 66, 65, 965 };
	pfverts[1346] = { 1132, 330, -1, 1128 };
	pfverts[1347] = { 618, 1348, -1, 601 };
	pfverts[1348] = { 52, 242, -1, 1347 };
	pfverts[1349] = { -1, 772, 637, -1 };
	pfverts[1350] = { 771, -1, 637, 727 };
	pfverts[1351] = { 1050, 201, -1, 1216 };
	pfverts[1352] = { 700, -1, 1176, 1186 };
	pfverts[1353] = { 408, 399, -1, 1074 };
	pfverts[1354] = { 1362, 624, 1313, 603 };
	pfverts[1355] = { 1356, 1285, -1, 98 };
	pfverts[1356] = { 375, 1355, -1, 1079 };
	pfverts[1357] = { 101, 1303, -1, 500 };
	pfverts[1358] = { 1360, 1359, -1, 695 };
	pfverts[1359] = { 1358, 709, -1, 1188 };
	pfverts[1360] = { 1358, 6, 708, 1270 };
	pfverts[1361] = { 1322, -1, 709, 676 };
	pfverts[1362] = { -1, 1354, 32, 1327 };
  // triangular elements
	pfverts[1363] = { 1730, 1451, 1728 };
	pfverts[1364] = { 1780, 1515, 1730 };
	pfverts[1365] = { 1400, 1399, 1403 };
	pfverts[1366] = { 1533, 1651, 1403 };
	pfverts[1367] = { 1505, 1523, 1506 };
	pfverts[1368] = { 1506, 1762, 1763 };
	pfverts[1369] = { 1757, 1754, 1529 };
	pfverts[1370] = { 1494, 1492, 1754 };
	pfverts[1371] = { 1561, 1742, 1486 };
	pfverts[1372] = { 1441, 1636, 1443 };
	pfverts[1373] = { 1445, 1781, 1443 };
	pfverts[1374] = { 1606, 1592, 1593 };
	pfverts[1375] = { 1593, 1698, 1703 };
	pfverts[1376] = { 1670, 1396, 1397 };
	pfverts[1377] = { 1628, 1626, 1670 };
	pfverts[1378] = { 1510, 1513, 1677 };
	pfverts[1379] = { 1609, 1612, 1678 };
	pfverts[1380] = { 1511, 1514, 1512 };
	pfverts[1381] = { 1620, 1627, 1512 };
	pfverts[1382] = { 1402, 1401, 1646 };
	pfverts[1383] = { 1646, 1649, 1647 };
	pfverts[1384] = { 1691, 1385, 1694 };
	pfverts[1385] = { 1384, 1697, 1700 };
	pfverts[1386] = { 1395, 1393, 1387 };
	pfverts[1387] = { 1741, 1386, 1743 };
	pfverts[1388] = { 1759, 1503, 1500 };
	pfverts[1389] = { 1501, 1759, 1763 };
	pfverts[1390] = { 1497, 1391, 1392 };
	pfverts[1391] = { 1764, 1390, 1758 };
	pfverts[1392] = { 1751, 1390, 1753 };
	pfverts[1393] = { 1749, 1744, 1386 };
	pfverts[1394] = { 1738, 1740, 1395 };
	pfverts[1395] = { 1477, 1394, 1386 };
	pfverts[1396] = { 1376, 1675, 1674 };
	pfverts[1397] = { 1669, 1398, 1376 };
	pfverts[1398] = { 1616, 1673, 1397 };
	pfverts[1399] = { 1405, 1410, 1365 };
	pfverts[1400] = { 1657, 1401, 1365 };
	pfverts[1401] = { 1400, 1650, 1382 };
	pfverts[1402] = { 1652, 1382, 1403 };
	pfverts[1403] = { 1366, 1402, 1365 };
	pfverts[1404] = { 1655, 1406, 1407 };
	pfverts[1405] = { 1654, 1406, 1399 };
	pfverts[1406] = { 1654, 1405, 1404 };
	pfverts[1407] = { 1568, 1404, 1409 };
	pfverts[1408] = { 1522, 1410, 1409 };
	pfverts[1409] = { 1524, 1407, 1408 };
	pfverts[1410] = { 1408, 1533, 1399 };
	pfverts[1411] = { 1634, 1412, 1642 };
	pfverts[1412] = { 1638, 1414, 1411 };
	pfverts[1413] = { 1650, 1656, 1642 };
	pfverts[1414] = { 1412, 1511, 1649 };
	pfverts[1415] = { 1629, 1416, 1419 };
	pfverts[1416] = { 1645, 1417, 1415 };
	pfverts[1417] = { 1647, 1418, 1416 };
	pfverts[1418] = { 1628, 1419, 1417 };
	pfverts[1419] = { 1671, 1415, 1418 };
	pfverts[1420] = { 1617, 1424, 1421 };
	pfverts[1421] = { 1674, 1420, 1422 };
	pfverts[1422] = { 1673, 1423, 1421 };
	pfverts[1423] = { 1616, 1425, 1422 };
	pfverts[1424] = { 1618, 1425, 1420 };
	pfverts[1425] = { 1543, 1424, 1423 };
	pfverts[1426] = { 1613, 1427, 1428 };
	pfverts[1427] = { 1620, 1430, 1426 };
	pfverts[1428] = { 1598, 1426, 1429 };
	pfverts[1429] = { 1689, 1428, 1431 };
	pfverts[1430] = { 1510, 1431, 1427 };
	pfverts[1431] = { 1680, 1429, 1430 };
	pfverts[1432] = { 1590, 1434, 1436 };
	pfverts[1433] = { 1699, 1437, 1435 };
	pfverts[1434] = { 1714, 1438, 1432 };
	pfverts[1435] = { 1701, 1436, 1433 };
	pfverts[1436] = { 1594, 1432, 1435 };
	pfverts[1437] = { 1575, 1438, 1433 };
	pfverts[1438] = { 1713, 1434, 1437 };
	pfverts[1439] = { 1729, 1442, 1443 };
	pfverts[1440] = { 1665, 1441, 1442 };
	pfverts[1441] = { 1632, 1372, 1440 };
	pfverts[1442] = { 1587, 1440, 1439 };
	pfverts[1443] = { 1373, 1439, 1372 };
	pfverts[1444] = { 1655, 1448, 1449 };
	pfverts[1445] = { 1373, 1449, 1446 };
	pfverts[1446] = { 1779, 1445, 1447 };
	pfverts[1447] = { 1518, 1446, 1448 };
	pfverts[1448] = { 1568, 1447, 1444 };
	pfverts[1449] = { 1643, 1444, 1445 };
	pfverts[1450] = { 1708, 1451, 1704 };
	pfverts[1451] = { 1363, 1452, 1450 };
	pfverts[1452] = { 1515, 1583, 1451 };
	pfverts[1453] = { 1617, 1614, 1455 };
	pfverts[1454] = { 1538, 1455, 1456 };
	pfverts[1455] = { 1619, 1453, 1454 };
	pfverts[1456] = { 1732, 1454, 1457 };
	pfverts[1457] = { 1598, 1456, 1614 };
	pfverts[1458] = { 1564, 1460, 1461 };
	pfverts[1459] = { 1732, 1463, 1460 };
	pfverts[1460] = { 1733, 1459, 1458 };
	pfverts[1461] = { 1686, 1458, 1462 };
	pfverts[1462] = { 1688, 1461, 1463 };
	pfverts[1463] = { 1598, 1462, 1459 };
	pfverts[1464] = { 1540, 1468, 1469 };
	pfverts[1465] = { 1733, 1466, 1468 };
	pfverts[1466] = { 1564, 1467, 1465 };
	pfverts[1467] = { 1579, 1469, 1466 };
	pfverts[1468] = { 1535, 1465, 1464 };
	pfverts[1469] = { 1574, 1464, 1467 };
	pfverts[1470] = { 1575, 1473, 1472 };
	pfverts[1471] = { 1492, 1472, 1529 };
	pfverts[1472] = { 1548, 1470, 1471 };
	pfverts[1473] = { 1584, 1529, 1470 };
	pfverts[1474] = { 1747, 1475, 1478 };
	pfverts[1475] = { 1746, 1476, 1474 };
	pfverts[1476] = { 1746, 1475, 1477 };
	pfverts[1477] = { 1395, 1476, 1479 };
	pfverts[1478] = { 1544, 1474, 1480 };
	pfverts[1479] = { 1543, 1477, 1480 };
	pfverts[1480] = { 1618, 1478, 1479 };
	pfverts[1481] = { 1718, 1482, 1484 };
	pfverts[1482] = { 1559, 1485, 1481 };
	pfverts[1483] = { 1745, 1748, 1742 };
	pfverts[1484] = { 1573, 1481, 1742 };
	pfverts[1485] = { 1482, 1548, 1745 };
	pfverts[1486] = { 1371, 1747, 1489 };
	pfverts[1487] = { 1489, 1536, 1561 };
	pfverts[1488] = { 1544, 1777, 1489 };
	pfverts[1489] = { 1487, 1486, 1488 };
	pfverts[1490] = { 1739, 1492, 1491 };
	pfverts[1491] = { 1741, 1490, 1493 };
	pfverts[1492] = { 1490, 1471, 1370 };
	pfverts[1493] = { 1555, 1491, 1494 };
	pfverts[1494] = { 1755, 1493, 1370 };
	pfverts[1495] = { 1496, 1762, 1498 };
	pfverts[1496] = { 1761, 1497, 1495 };
	pfverts[1497] = { 1390, 1499, 1496 };
	pfverts[1498] = { 1495, 1523, 1524 };
	pfverts[1499] = { 1522, 1524, 1497 };
	pfverts[1500] = { 1581, 1388, 1582 };
	pfverts[1501] = { 1389, 1760, 1503 };
	pfverts[1502] = { 1694, 1582, 1504 };
	pfverts[1503] = { 1388, 1501, 1757 };
	pfverts[1504] = { 1692, 1502, 1757 };
	pfverts[1505] = { 1517, 1367, 1507 };
	pfverts[1506] = { 1368, 1520, 1367 };
	pfverts[1507] = { 1515, 1505, 1508 };
	pfverts[1508] = { 1581, 1507, 1520 };
	pfverts[1509] = { 1640, 1513, 1514 };
	pfverts[1510] = { 1430, 1512, 1378 };
	pfverts[1511] = { 1380, 1648, 1414 };
	pfverts[1512] = { 1381, 1380, 1510 };
	pfverts[1513] = { 1509, 1612, 1378 };
	pfverts[1514] = { 1638, 1509, 1380 };
	pfverts[1515] = { 1507, 1452, 1364 };
	pfverts[1516] = { 1782, 1519, 0 };
	pfverts[1517] = { 1780, 1518, 1505 };
	pfverts[1518] = { 1447, 1521, 1517 };
	pfverts[1519] = { 1516, 1527, 0 };
	pfverts[1520] = { 1759, 1508, 1506 };
	pfverts[1521] = { 1518, 1568, 1523 };
	pfverts[1522] = { 1531, 1408, 1499 };
	pfverts[1523] = { 1521, 1498, 1367 };
	pfverts[1524] = { 1499, 1498, 1409 };
	pfverts[1525] = { 1526, 1545, 0 };
	pfverts[1526] = { 1525, 1528, 1766 };
	pfverts[1527] = { 1519, 1567, 1528 };
	pfverts[1528] = { 0, 1527, 1526 };
	pfverts[1529] = { 1369, 1471, 1473 };
	pfverts[1530] = { 1553, 1532, 1750 };
	pfverts[1531] = { 1522, 1751, 1653 };
	pfverts[1532] = { 1569, 1653, 1530 };
	pfverts[1533] = { 1410, 1653, 1366 };
	pfverts[1534] = { 1752, 1755, 1554 };
	pfverts[1535] = { 1734, 1468, 1536 };
	pfverts[1536] = { 1777, 1535, 1487 };
	pfverts[1537] = { 1778, 0, 1539 };
	pfverts[1538] = { 1776, 1541, 1454 };
	pfverts[1539] = { 1537, 0, 1550 };
	pfverts[1540] = { 1464, 1572, 1561 };
	pfverts[1541] = { 1538, 1544, 1619 };
	pfverts[1542] = { 1616, 1571, 1557 };
	pfverts[1543] = { 1557, 1479, 1425 };
	pfverts[1544] = { 1541, 1488, 1478 };
	pfverts[1545] = { 1525, 1775, 1547 };
	pfverts[1546] = { 1552, 0, 1547 };
	pfverts[1547] = { 1545, 1546, 0 };
	pfverts[1548] = { 1472, 1739, 1485 };
	pfverts[1549] = { 1778, 1550, 1563 };
	pfverts[1550] = { 1539, 1551, 1549 };
	pfverts[1551] = { 0, 1552, 1550 };
	pfverts[1552] = { 1546, 1562, 1551 };
	pfverts[1553] = { 1530, 1554, 1556 };
	pfverts[1554] = { 1534, 1555, 1553 };
	pfverts[1555] = { 1554, 1493, 1737 };
	pfverts[1556] = { 1553, 1558, 1569 };
	pfverts[1557] = { 1543, 1542, 1738 };
	pfverts[1558] = { 1556, 1736, 1570 };
	pfverts[1559] = { 1482, 1560, 1575 };
	pfverts[1560] = { 1717, 1713, 1559 };
	pfverts[1561] = { 1487, 1540, 1371 };
	pfverts[1562] = { 1552, 1771, 1563 };
	pfverts[1563] = { 1549, 1562, 1595 };
	pfverts[1564] = { 1576, 1466, 1458 };
	pfverts[1565] = { 1778, 1596, 0 };
	pfverts[1566] = { 1567, 1630, 1631 };
	pfverts[1567] = { 1527, 1566, 1767 };
	pfverts[1568] = { 1521, 1448, 1407 };
	pfverts[1569] = { 1532, 1556, 1622 };
	pfverts[1570] = { 1558, 1571, 1624 };
	pfverts[1571] = { 1542, 1621, 1570 };
	pfverts[1572] = { 1540, 1574, 1573 };
	pfverts[1573] = { 1484, 1572, 1716 };
	pfverts[1574] = { 1572, 1469, 1726 };
	pfverts[1575] = { 1559, 1437, 1470 };
	pfverts[1576] = { 1564, 1580, 1577 };
	pfverts[1577] = { 1576, 1601, 1686 };
	pfverts[1578] = { 1601, 1580, 1720 };
	pfverts[1579] = { 1467, 1722, 1580 };
	pfverts[1580] = { 1576, 1579, 1578 };
	pfverts[1581] = { 1508, 1500, 1583 };
	pfverts[1582] = { 1693, 1500, 1502 };
	pfverts[1583] = { 1581, 1710, 1452 };
	pfverts[1584] = { 1473, 1699, 1692 };
	pfverts[1585] = { 1586, 1665, 1589 };
	pfverts[1586] = { 1585, 1704, 1587 };
	pfverts[1587] = { 1442, 1728, 1586 };
	pfverts[1588] = { 1663, 1593, 1589 };
	pfverts[1589] = { 1585, 1588, 1705 };
	pfverts[1590] = { 1712, 1432, 1591 };
	pfverts[1591] = { 1590, 1594, 1599 };
	pfverts[1592] = { 1600, 1696, 1374 };
	pfverts[1593] = { 1588, 1374, 1375 };
	pfverts[1594] = { 1436, 1696, 1591 };
	pfverts[1595] = { 1563, 1597, 1596 };
	pfverts[1596] = { 1565, 1595, 0 };
	pfverts[1597] = { 1770, 0, 1595 };
	pfverts[1598] = { 1457, 1428, 1463 };
	pfverts[1599] = { 1591, 1600, 1604 };
	pfverts[1600] = { 1592, 1605, 1599 };
	pfverts[1601] = { 1577, 1578, 1602 };
	pfverts[1602] = { 1684, 1601, 1603 };
	pfverts[1603] = { 1711, 1604, 1602 };
	pfverts[1604] = { 1599, 1608, 1603 };
	pfverts[1605] = { 1600, 1607, 1609 };
	pfverts[1606] = { 1374, 1660, 1607 };
	pfverts[1607] = { 1605, 1606, 1610 };
	pfverts[1608] = { 1604, 1609, 1676 };
	pfverts[1609] = { 1605, 1379, 1608 };
	pfverts[1610] = { 1607, 1611, 1612 };
	pfverts[1611] = { 1659, 1640, 1610 };
	pfverts[1612] = { 1610, 1513, 1379 };
	pfverts[1613] = { 1426, 1614, 1615 };
	pfverts[1614] = { 1613, 1457, 1453 };
	pfverts[1615] = { 1613, 1617, 1675 };
	pfverts[1616] = { 1542, 1423, 1398 };
	pfverts[1617] = { 1615, 1453, 1420 };
	pfverts[1618] = { 1480, 1424, 1619 };
	pfverts[1619] = { 1541, 1618, 1455 };
	pfverts[1620] = { 1427, 1626, 1381 };
	pfverts[1621] = { 1571, 1669, 1623 };
	pfverts[1622] = { 1569, 1624, 1651 };
	pfverts[1623] = { 1621, 1625, 1624 };
	pfverts[1624] = { 1570, 1623, 1622 };
	pfverts[1625] = { 1629, 1644, 1623 };
	pfverts[1626] = { 1620, 1675, 1377 };
	pfverts[1627] = { 1628, 1648, 1381 };
	pfverts[1628] = { 1377, 1418, 1627 };
	pfverts[1629] = { 1625, 1672, 1415 };
	pfverts[1630] = { 1782, 0, 1566 };
	pfverts[1631] = { 1768, 1566, 0 };
	pfverts[1632] = { 1668, 1633, 1441 };
	pfverts[1633] = { 1668, 1635, 1632 };
	pfverts[1634] = { 1411, 1636, 1635 };
	pfverts[1635] = { 1637, 1634, 1633 };
	pfverts[1636] = { 1643, 1372, 1634 };
	pfverts[1637] = { 1635, 1667, 1639 };
	pfverts[1638] = { 1412, 1639, 1514 };
	pfverts[1639] = { 1637, 1641, 1638 };
	pfverts[1640] = { 1611, 1509, 1641 };
	pfverts[1641] = { 1658, 1640, 1639 };
	pfverts[1642] = { 1411, 1413, 1643 };
	pfverts[1643] = { 1636, 1642, 1449 };
	pfverts[1644] = { 1625, 1645, 1652 };
	pfverts[1645] = { 1416, 1644, 1646 };
	pfverts[1646] = { 1383, 1645, 1382 };
	pfverts[1647] = { 1417, 1383, 1648 };
	pfverts[1648] = { 1627, 1647, 1511 };
	pfverts[1649] = { 1383, 1650, 1414 };
	pfverts[1650] = { 1401, 1413, 1649 };
	pfverts[1651] = { 1622, 1652, 1366 };
	pfverts[1652] = { 1644, 1402, 1651 };
	pfverts[1653] = { 1531, 1532, 1533 };
	pfverts[1654] = { 1405, 1406, 1657 };
	pfverts[1655] = { 1444, 1656, 1404 };
	pfverts[1656] = { 1413, 1657, 1655 };
	pfverts[1657] = { 1400, 1654, 1656 };
	pfverts[1658] = { 1641, 1659, 1664 };
	pfverts[1659] = { 1611, 1660, 1658 };
	pfverts[1660] = { 1606, 1662, 1659 };
	pfverts[1661] = { 1666, 1663, 1665 };
	pfverts[1662] = { 1660, 1663, 1664 };
	pfverts[1663] = { 1588, 1661, 1662 };
	pfverts[1664] = { 1658, 1662, 1667 };
	pfverts[1665] = { 1440, 1661, 1585 };
	pfverts[1666] = { 1667, 1661, 1668 };
	pfverts[1667] = { 1637, 1664, 1666 };
	pfverts[1668] = { 1633, 1666, 1632 };
	pfverts[1669] = { 1621, 1397, 1672 };
	pfverts[1670] = { 1377, 1376, 1671 };
	pfverts[1671] = { 1419, 1670, 1672 };
	pfverts[1672] = { 1629, 1669, 1671 };
	pfverts[1673] = { 1422, 1674, 1398 };
	pfverts[1674] = { 1421, 1673, 1396 };
	pfverts[1675] = { 1626, 1615, 1396 };
	pfverts[1676] = { 1608, 1679, 1685 };
	pfverts[1677] = { 1378, 1678, 1680 };
	pfverts[1678] = { 1379, 1677, 1679 };
	pfverts[1679] = { 1676, 1678, 1681 };
	pfverts[1680] = { 1677, 1681, 1431 };
	pfverts[1681] = { 1679, 1680, 1683 };
	pfverts[1682] = { 1684, 1687, 1686 };
	pfverts[1683] = { 1681, 1690, 1685 };
	pfverts[1684] = { 1602, 1685, 1682 };
	pfverts[1685] = { 1676, 1683, 1684 };
	pfverts[1686] = { 1577, 1682, 1461 };
	pfverts[1687] = { 1682, 1690, 1688 };
	pfverts[1688] = { 1462, 1689, 1687 };
	pfverts[1689] = { 1429, 1690, 1688 };
	pfverts[1690] = { 1683, 1689, 1687 };
	pfverts[1691] = { 1702, 1384, 1692 };
	pfverts[1692] = { 1584, 1691, 1504 };
	pfverts[1693] = { 1582, 1695, 1707 };
	pfverts[1694] = { 1384, 1695, 1502 };
	pfverts[1695] = { 1706, 1693, 1694 };
	pfverts[1696] = { 1592, 1594, 1698 };
	pfverts[1697] = { 1385, 1706, 1703 };
	pfverts[1698] = { 1375, 1696, 1701 };
	pfverts[1699] = { 1584, 1433, 1702 };
	pfverts[1700] = { 1385, 1701, 1702 };
	pfverts[1701] = { 1698, 1435, 1700 };
	pfverts[1702] = { 1691, 1699, 1700 };
	pfverts[1703] = { 1375, 1697, 1705 };
	pfverts[1704] = { 1586, 1705, 1450 };
	pfverts[1705] = { 1589, 1703, 1704 };
	pfverts[1706] = { 1695, 1708, 1697 };
	pfverts[1707] = { 1693, 1709, 1710 };
	pfverts[1708] = { 1450, 1706, 1709 };
	pfverts[1709] = { 1710, 1708, 1707 };
	pfverts[1710] = { 1583, 1707, 1709 };
	pfverts[1711] = { 1603, 1719, 1712 };
	pfverts[1712] = { 1590, 1711, 1714 };
	pfverts[1713] = { 1560, 1715, 1438 };
	pfverts[1714] = { 1712, 1715, 1434 };
	pfverts[1715] = { 1721, 1713, 1714 };
	pfverts[1716] = { 1573, 1727, 1718 };
	pfverts[1717] = { 1560, 1718, 1723 };
	pfverts[1718] = { 1481, 1716, 1717 };
	pfverts[1719] = { 1711, 1720, 1721 };
	pfverts[1720] = { 1578, 1722, 1719 };
	pfverts[1721] = { 1715, 1719, 1723 };
	pfverts[1722] = { 1579, 1724, 1720 };
	pfverts[1723] = { 1717, 1721, 1725 };
	pfverts[1724] = { 1722, 1726, 1725 };
	pfverts[1725] = { 1723, 1724, 1727 };
	pfverts[1726] = { 1574, 1724, 1727 };
	pfverts[1727] = { 1716, 1726, 1725 };
	pfverts[1728] = { 1587, 1729, 1363 };
	pfverts[1729] = { 1439, 1731, 1728 };
	pfverts[1730] = { 1364, 1363, 1731 };
	pfverts[1731] = { 1781, 1730, 1729 };
	pfverts[1732] = { 1456, 1459, 1735 };
	pfverts[1733] = { 1465, 1734, 1460 };
	pfverts[1734] = { 1535, 1735, 1733 };
	pfverts[1735] = { 1776, 1732, 1734 };
	pfverts[1736] = { 1558, 1737, 1738 };
	pfverts[1737] = { 1555, 1740, 1736 };
	pfverts[1738] = { 1557, 1736, 1394 };
	pfverts[1739] = { 1548, 1490, 1743 };
	pfverts[1740] = { 1737, 1741, 1394 };
	pfverts[1741] = { 1491, 1740, 1387 };
	pfverts[1742] = { 1484, 1483, 1371 };
	pfverts[1743] = { 1387, 1744, 1739 };
	pfverts[1744] = { 1393, 1745, 1743 };
	pfverts[1745] = { 1483, 1485, 1744 };
	pfverts[1746] = { 1476, 1475, 1749 };
	pfverts[1747] = { 1486, 1748, 1474 };
	pfverts[1748] = { 1483, 1749, 1747 };
	pfverts[1749] = { 1393, 1746, 1748 };
	pfverts[1750] = { 1530, 1751, 1752 };
	pfverts[1751] = { 1531, 1392, 1750 };
	pfverts[1752] = { 1534, 1750, 1753 };
	pfverts[1753] = { 1752, 1392, 1756 };
	pfverts[1754] = { 1369, 1760, 1370 };
	pfverts[1755] = { 1534, 1494, 1756 };
	pfverts[1756] = { 1753, 1755, 1758 };
	pfverts[1757] = { 1504, 1503, 1369 };
	pfverts[1758] = { 1756, 1760, 1391 };
	pfverts[1759] = { 1520, 1389, 1388 };
	pfverts[1760] = { 1754, 1501, 1758 };
	pfverts[1761] = { 1496, 1764, 1762 };
	pfverts[1762] = { 1761, 1368, 1495 };
	pfverts[1763] = { 1764, 1389, 1368 };
	pfverts[1764] = { 1391, 1763, 1761 };
	pfverts[1765] = { 0, 1769, 1768 };
	pfverts[1766] = { 1526, 1767, 1774 };
	pfverts[1767] = { 1567, 1768, 1766 };
	pfverts[1768] = { 1631, 1765, 1767 };
	pfverts[1769] = { 0, 1773, 1765 };
	pfverts[1770] = { 1597, 1771, 1772 };
	pfverts[1771] = { 1562, 1775, 1770 };
	pfverts[1772] = { 0, 1770, 1773 };
	pfverts[1773] = { 1769, 1772, 1774 };
	pfverts[1774] = { 1766, 1773, 1775 };
	pfverts[1775] = { 1545, 1774, 1771 };
	pfverts[1776] = { 1777, 1538, 1735 };
	pfverts[1777] = { 1536, 1488, 1776 };
	pfverts[1778] = { 1549, 1565, 1537 };
	pfverts[1779] = { 1446, 1780, 1781 };
	pfverts[1780] = { 1517, 1364, 1779 };
	pfverts[1781] = { 1373, 1779, 1731 };
	pfverts[1782] = { 1630, 1516, 0 };
	pfverts[1783] = { 1809, 1823 };
	pfverts[1784] = { 1837, 1824 };
	pfverts[1785] = { 1848, 1786 };
	pfverts[1786] = { 1787, 1785 };
	pfverts[1787] = { 1786, 1849 };
	pfverts[1788] = { 1810, 1849 };
	pfverts[1789] = { 1790, 1815 };
	pfverts[1790] = { 1791, 1789 };
	pfverts[1791] = { 1790, 1829 };
	pfverts[1792] = { -30, 1823 };
	pfverts[1793] = { 1842, 1814 };
	pfverts[1794] = { 1851, 1825 };
	pfverts[1795] = { 1843, 1830 };
	pfverts[1796] = { 1809, 1820 };
	pfverts[1797] = { 1798, 1844 };
	pfverts[1798] = { 1822, 1797 };
	pfverts[1799] = { 1800, 1801 };
	pfverts[1800] = { 1799, 1819 };
	pfverts[1801] = { 1807, 1799 };
	pfverts[1802] = { 1822, 1850 };
	pfverts[1803] = { 1804, 1805 };
	pfverts[1804] = { 1803, 1806 };
	pfverts[1805] = { 1803, 1818 };
	pfverts[1806] = { 1804, 1808 };
	pfverts[1807] = { 1808, 1801 };
	pfverts[1808] = { 1807, 1806 };
	pfverts[1809] = { 1783, 1796 };
	pfverts[1810] = { 1811, 1788 };
	pfverts[1811] = { 1810, 1837 };
	pfverts[1812] = { -30, 1843 };
	pfverts[1813] = { 1816, 1838 };
	pfverts[1814] = { 1848, 1793 };
	pfverts[1815] = { 1789, 1847 };
	pfverts[1816] = { 1817, 1813 };
	pfverts[1817] = { 1816, 1824 };
	pfverts[1818] = { 1850, 1805 };
	pfverts[1819] = { 1800, 1830 };
	pfverts[1820] = { 1821, 1796 };
	pfverts[1821] = { 1820, 1836 };
	pfverts[1822] = { 1802, 1798 };
	pfverts[1823] = { 1792, 1783 };
	pfverts[1824] = { 1817, 1784 };
	pfverts[1825] = { 1841, 1794 };
	pfverts[1826] = { 1831, 1842 };
	pfverts[1827] = { 1828, 1835 };
	pfverts[1828] = { 1827, -30 };
	pfverts[1829] = { 1791, 1844 };
	pfverts[1830] = { 1819, 1795 };
	pfverts[1831] = { 1826, 1839 };
	pfverts[1832] = { 1833, 1839 };
	pfverts[1833] = { 1832, 1845 };
	pfverts[1834] = { 1835, -30 };
	pfverts[1835] = { 1834, 1827 };
	pfverts[1836] = { 1821, -30 };
	pfverts[1837] = { 1784, 1811 };
	pfverts[1838] = { 1852, 1813 };
	pfverts[1839] = { 1831, 1832 };
	pfverts[1840] = { 1841, -30 };
	pfverts[1841] = { 1840, 1825 };
	pfverts[1842] = { 1826, 1793 };
	pfverts[1843] = { 1795, 1812 };
	pfverts[1844] = { 1797, 1829 };
	pfverts[1845] = { 1846, 1833 };
	pfverts[1846] = { 1845, -30 };
	pfverts[1847] = { 1853, 1815 };
	pfverts[1848] = { 1814, 1785 };
	pfverts[1849] = { 1787, 1788 };
	pfverts[1850] = { 1802, 1818 };
	pfverts[1851] = { 1852, 1794 };
	pfverts[1852] = { 1851, 1838 };
	pfverts[1853] = { -30, 1847 };
	pfverts[1854] = { 1855, 1870 };
	pfverts[1855] = { 1854, 1861 };
	pfverts[1856] = { 1867, 1863 };
	pfverts[1857] = { 1870, 1869 };
	pfverts[1858] = { 1859, 1862 };
	pfverts[1859] = { 1858, 1871 };
	pfverts[1860] = { 1861, 1863 };
	pfverts[1861] = { 1855, 1860 };
	pfverts[1862] = { 1858, 1864 };
	pfverts[1863] = { 1860, 1856 };
	pfverts[1864] = { 1862, 1868 };
	pfverts[1865] = { 1869, 1871 };
	pfverts[1866] = { 1867, 1868 };
	pfverts[1867] = { 1866, 1856 };
	pfverts[1868] = { 1866, 1864 };
	pfverts[1869] = { 1857, 1865 };
	pfverts[1870] = { 1854, 1857 };
	pfverts[1871] = { 1859, 1865 };

  // checking whether pfverts has the expected size
  size_t n_pfvert_entries{0};
  for ( auto et : etypes )
    n_pfvert_entries += CSMP_ElementSpecifications::NodesPerElementOfType( et );

  size_t n_pfvert_entries_actual{0};
  for ( auto pt : pfverts )
    n_pfvert_entries_actual += pt.size();
  assert( n_pfvert_entries_actual == n_pfvert_entries ); // all good

  vset.AddPfverts( pfverts.begin(), pfverts.end());
  
  
  // boundary flags
  const int8_t bf(IRREGULAR);
  // side boundaries to start with
  vector<int8_t>  bflags{bf,bf,0,0,0,bf,0,bf,0,bf,bf,0,0,0,0,0,bf,bf,0,bf,0,bf,0,0,0,0,0,0,0,0,0,0,0,0,0,0,bf,0,bf,0,bf,0,bf,bf,bf,0,bf,0,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,bf,bf,bf,bf,0,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,0,bf,0,0,bf,bf,bf,bf,bf,bf,bf,0,0,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,bf,0,0,0,0,bf,0,0,bf,bf,0,0,0,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,bf,bf,bf,0,0,0,0,bf,bf,0,bf,bf,0,bf,bf,bf,bf,0,0,bf,0,bf,bf,bf,bf,0,0,0,0,bf,bf,0,0,bf,bf,0,0,0,0,0,0,0,bf,bf,bf,0,0,0,bf,bf,bf,bf,bf,0,bf,bf,bf,0,0,bf,0,0,0,0,0,0,0,0,0,bf,bf,bf,bf,0,0,0,bf,0,0,0,bf,bf,bf,0,0,0,bf,0,bf,bf,bf,bf,0,bf,bf,bf,0,bf,bf,bf,bf,0,0,0,0,bf,0,0,0,0,0,0,0,bf,bf,bf,bf,0,0,0,bf,0,0,0,bf,bf,bf,0,0,0,bf,0,bf,bf,bf,0,bf,bf,0,0,0,0,0,bf,0,0,0,bf,bf,0,bf,0,bf,bf,0,bf,0,bf,bf,0,bf,bf,0,bf,bf,0,bf,0,0,bf,0,bf,bf,bf,bf,bf,bf,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf};

  // INTERNAL boundaries inferred from lower-dimensional elements that are not part of the side boundaries
  for ( const auto& eit : plist ) {
       // identifying interior elements
       bool at_boundary{false};
       for ( const auto& n : eit )
         if ( bflags[n] != NOT ) {
              at_boundary = true;
              break;
           }
       if ( !at_boundary )
         for ( const auto& n : eit )
           bflags[n] = INTERNAL;
    }

  vset.AddBFlags( bflags.begin(), bflags.end() );

  // geometry flags
  vector<int8_t>  gflags( bflags.size(), MESH_VERTEX );
  for ( size_t i{0U}; i<bflags.size(); i++ ) {
       if (      bflags[i] == IRREGULAR ) gflags[i] = EXTERIOR_SURFACE;
       else if ( bflags[i] == INTERNAL )  gflags[i] = INTERIOR_SURFACE;
    }
  
  vset.AddBREP_Flags( gflags.begin(), gflags.end() );

  
// MATERIAL PROPERTIES

  const int32_t   material_identifier{1};
  vector<int32_t> pmtrl( vset.Elements(), material_identifier );
  vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

  PropertyData elmt_nums( ELEMENT, SCALAR, 3U );
  elmt_nums.Reserve( vset.Elements() );

  for ( size_t i{0U}; i<vset.Elements(); ++i )
    pushBack( elmt_nums, makeScalar( ANY, static_cast<double>(i) ) );
  vset.AddData( "element number", elmt_nums );

  PropertyData node_nums( NODE, SCALAR, 3U );
  node_nums.Reserve( vset.Vertices() );

  for ( size_t i{0U}; i<vset.Vertices(); ++i )
    pushBack( node_nums, makeScalar( ANY, static_cast<double>(i) ) );
  vset.AddData( "node number", node_nums );

} // end create_FracBox



} // end csmp

/**
    USAGE EXAMPLE:

//Initialize Model:

VSet<3U>    mesh_container;
create_Pyramid_VSet(mesh_container, true);
Model<3U>  model3D( mesh_container, true ); 
model3D.AddFaceTopology();
mesh_container.Erase();

//Initialize properties (depends on what you need...):
//model configuration
model3D.InputPropertyValue ( "permeability", 1.0e-12 );
model3D.InputPropertyValue ( "porosity", 0.25 );
model3D.InputPropertyValue ( "storativity", 1.0e-9 );
model3D.InputPropertyValue ( "fluid volume source", 0. );
model3D.InputPropertyValue ( "nodal fluid volume source", 0. );
model3D.InputPropertyValue ( "concentration", 0. );
model3D.InputPropertyValue ( "diffusivity", 1.0e-50 );

//remember to input boundary conditions as required

*/
