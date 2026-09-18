// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
    vecNeighbors[1]=1;
    vecNeighbors[2]=3;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[0]=vecNeighbors;
    //element 1
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=2;
    vecNeighbors[3]=0;
    deqElementNeighbors[1]=vecNeighbors;
    //element 2
    vecNeighbors[0]=1;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=3;
    deqElementNeighbors[2]=vecNeighbors;
    //element 3
    vecNeighbors[0]=0;
    vecNeighbors[1]=2;
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

//    vset.Out();
    cout <<"\n"<<"create_QuadraticQUadrilateral_VSet: model 'unnamed': done."<< endl;
    
    return vset;
    
} // end create_QuadraticQUadrilateral_VSet




/**
    Unit square. Only corners! - to test neighbor flagging.
      
    3               2
     o--------o
     | \     1   o       element 1
     |   \        o
     |     \      o
     |       \    o       element 0
     |    0   \  o
     o--------o
    0                1
 */
void create_2_Triangle_VSet(VSet<2U>& vset )
{
  	IsoparametricLinearTriangle iso_tria;
  	
    const size_t   n_nodes(4), n_elmts{2};
    const uint32_t nodes_per_elmt{3}, nbors_per_elmt{3};
  	vset.Resize( nodes_per_elmt, nbors_per_elmt, ISOPARAMETRIC_LINEAR_TRIANGLE, n_nodes, n_elmts );
 
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int8_t> vecElementTypes = { ISOPARAMETRIC_LINEAR_TRIANGLE, ISOPARAMETRIC_LINEAR_TRIANGLE };
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );
 
  	//-----------------------NODES
  	//define nodes
  	deque<double> px(n_nodes);
  	deque<double> py(n_nodes);
  	deque<double> pz(n_nodes);
  	
  	px[0]=0;   py[0]=0.;    pz[0]=0.;
  	px[1]=1.;  py[1]=0.;    pz[1]=0.;
  	px[2]=1.;  py[2]=1.;    pz[2]=0.;
  	px[3]=0.;  py[3]=1.;    pz[3]=0.;
  	
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
    //define quadrilateral elements (elements 0->26), assign nodes per element
    deque< vector<size_t> > deqElements = { {0, 1, 3}, {1, 2, 3} };
    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors = { {1,LEFT_OUTSIDE,BOTTOM_OUTSIDE}, {TOP_OUTSIDE,0,RIGHT_OUTSIDE} };
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 0, CNR1 );
  	vset.BFlag( 1, CNR2 );
  	vset.BFlag( 2, CNR3 );
  	vset.BFlag( 3, CNR4 );
  	
    //vset.Out();
    cout <<"\n"<<"create_2_Triangle_VSet: model '2Triangle_VSet': done."<< endl;
    
} // end create_2_Triangle_VSet




/**
        No neighbors, but (negative) boundary flags.
*/
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
  	 
    deqElements[0][0]= 0;
    deqElements[0][1]= 1;
    deqElements[0][2]= 2;
    deqElements[0][3]= 3;
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
  	vset.BFlag( 0, CNR1 );
  	vset.BFlag( 1, CNR2 );
  	vset.BFlag( 2, CNR3 );
  	vset.BFlag( 3, CNR4 );
  	
    //vset.Out();
    cout <<"\n"<<"create_1Square_VSet: model 'unnamed': done."<< endl;
    
} // end create_1Square_VSet





/**
    Model TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object at the box boundary.
    Model is rectangle shaped
    
    @code
              [Top]            [Top]
         0 (0,1) ------- 1 (1,1) ------- 2 (2,1)
         |  \            |               |
         |    \ [Elmt0]  |               |
         |      \        |   [Elmt2]     |
  [Left] | [Elmt1]       |               | [Right]
         |        \      |               |
         |          \    |               |
         3 (0,0) ------- 4 (1,0) ------- 5 (2,0)
              [Bottom]       [Bottom]
    @endcode

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
    vector<int8_t> bflags = { CNR4, U, CNR3, CNR1, B, CNR2 };
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
    for ( size_t i{0U}; i<vset.Elements(); ++i ) pushBack( elmt_nums, makeScalar( ANY, static_cast<double>(i) ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t i{0U}; i<vset.Vertices(); ++i ) pushBack( node_nums, makeScalar( ANY, static_cast<double>(i) ) );
    vset.AddData( "node number", node_nums );
    // permeability
    PropertyData perm( ELEMENT, SCALAR, 2U );
    perm.Reserve( vset.Elements() );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("TRIA"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-15 ) );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("MIXED"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-14 ) );
    for ( size_t i{0U}; i<mesh_topology.CellsWithinDomain("LINE"); ++i ) pushBack( perm, makeScalar( ANY, 1.0e-12 ) );
    vset.AddData( "permeability", perm );

 //   vset.Out();
    cout <<"\n"<<"create_SimplePolyElement2DModel: model 'unnamed': done."<< endl;
    
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
    
    cout <<"\n"<<"create_TrianglePatch_VSet: model 'TrianglePatch': done."<< endl;
//    vset.Out();
    
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
  	px[9]=2.8;   py[9]=3.;    bflags[9] = INTERNAL;   gflags[9] = INTERIOR_POINT;
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
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end() );

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
      
     cout <<"\n"<<"create_MeshPatchWithLineElements_VSet: model 'un-named': done."<< endl;
   // vset.Out();
    
    return mesh_topology;
    
} // end create_MeshPatchWithLineElements_VSet





    /// 2D rectangular model the two halfs of each are offset from one another
void create_Disconnected2D_VSet(VSet<2U>& vset)
{
    // ----------------------- ELEMENT TYPES
    static constexpr size_t n_cells = 5;
    static constexpr size_t n_nodes = 12;
    static const deque<int8_t> vecElementTypes = {
        ISOPARAMETRIC_LINEAR_QUADRILATERAL, // Element 0
        ISOPARAMETRIC_LINEAR_QUADRILATERAL, // Element 1
        ISOPARAMETRIC_LINEAR_TRIANGLE,      // Element 2
        ISOPARAMETRIC_LINEAR_TRIANGLE,      // Element 3
        ISOPARAMETRIC_LINEAR_TRIANGLE       // Element 4 (corrected to triangle)
    };
    static const deque<uint32_t> npes = {4, 4, 3, 3, 3}; // Nodes per element
    static const deque<uint32_t> epes = {4, 4, 3, 3, 3}; // Neighbors per element
    assert(vecElementTypes.size() == n_cells && "Error: Element types size mismatch.");
    assert(npes.size() == n_cells && epes.size() == n_cells && "Error: npes/epes size mismatch.");

    vset.Resize(vecElementTypes, npes, epes, n_nodes, 0, 0);
    vset.AddElementTypes(vecElementTypes.begin(), vecElementTypes.end());

    // ----------------------- NODES
    static const deque<double> px = {1.0, 4.0, 1.0, 4.0, 1.0, 4.0, 5.0, 7.0, 5.0, 7.0, 5.0, 7.0};
    static const deque<double> py = {5.0, 5.0, 3.0, 3.0, 1.0, 1.0, 5.0, 5.0, 3.0, 3.0, 1.0, 1.0};
    deque<double> pz(n_nodes, 0.0);
    assert(px.size() == n_nodes && py.size() == n_nodes && pz.size() == n_nodes && "Error: Node coordinates size mismatch.");

    vset.AddXYZ(px, py, pz);

    // ----------------------- NODE BOUNDARY FLAGS
    static const vector<int8_t> bflags = {
        CNR4, TOP, LEFT, INTERNAL, CNR1, BOTTOM, // Nodes 0-5
        TOP, CNR3, INTERNAL, RIGHT, BOTTOM, CNR2  // Nodes 6-11
    };
    assert(bflags.size() == n_nodes && "Error: Boundary flags size mismatch.");

    vset.AddBFlags(bflags.begin(), bflags.end());

    // ----------------------- TOPOTYPE NODE FLAGS
    static const vector<int8_t> gflags = {
        MESH_VERTEX, MESH_VERTEX, MESH_VERTEX, EXTERIOR_POINT, // Nodes 0-3
        MESH_VERTEX, MESH_VERTEX, MESH_VERTEX, MESH_VERTEX,    // Nodes 4-7
        EXTERIOR_POINT, MESH_VERTEX, MESH_VERTEX, MESH_VERTEX  // Nodes 8-11
    };
    assert(gflags.size() == n_nodes && "Error: Topological flags size mismatch.");

    vset.AddBREP_Flags(gflags.begin(), gflags.end());

    // ----------------------- ELEMENTS
    static const deque<vector<size_t>> plist = {
        {2, 3, 1, 0}, // Quad 0: (1,3), (4,3), (4,5), (1,5)
        {4, 5, 3, 2}, // Quad 1: (1,1), (4,1), (4,3), (1,3)
        {6, 9, 7},    // Tri 2: (5,5), (7,3), (7,5)
        {6, 8, 9},    // Tri 3: (5,5), (5,3), (7,3)
        {10, 11, 9}   // Tri 4: (5,1), (7,1), (7,3)
    };
    assert(plist.size() == n_cells && "Error: Elements size mismatch.");

    vset.AddPlist(plist.begin(), plist.end());

    // ----------------------- NEIGHBORS
    static const deque<vector<int64_t>> pfverts = {
        {1, -INTERNAL, -TOP, -LEFT},        // Quad 0: to quad 1, internal, top, left
        {-BOTTOM, -INTERNAL, 0, -LEFT},     // Quad 1: bottom, internal, quad 0, left
        {4, -RIGHT, -TOP},                  // Tri 2: to tri 4, right, top
        {-RIGHT, -TOP, 4},                  // Tri 3: right, top, tri 4
        {-BOTTOM, -RIGHT, 3}                // Tri 4: bottom, right, tri 3
    };
    assert(pfverts.size() == n_cells && "Error: Neighbors size mismatch.");

    vset.AddPfverts(pfverts.begin(), pfverts.end());

    // ----------------------- MATERIALS
    static const vector<int32_t> pmtrl = {1, 1, 2, 2, 2};
    assert(pmtrl.size() == n_cells && "Error: Materials size mismatch.");

    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());

    // ----------------------- ELEMENT & NODE NUMBERS
    PropertyData elmt_nums(ELEMENT, SCALAR, 2U);
    elmt_nums.Reserve(vset.Elements());
    for (size_t n = 0; n < vset.Elements(); ++n) {
        pushBack( elmt_nums, makeScalar(ANY, static_cast<double>(n)) );
    }
    vset.AddData("element number", elmt_nums);

    PropertyData node_nums(NODE, SCALAR, 2U);
    node_nums.Reserve(vset.Vertices());
    for (size_t n = 0; n < vset.Vertices(); ++n) {
        pushBack( node_nums, makeScalar(ANY, static_cast<double>(n)) );
    }
    vset.AddData("node number", node_nums);

    // ----------------------- CONSISTENCY CHECKS
    for (size_t i = 0; i < n_cells; ++i) {
        assert(vset.ElementType(i) == vecElementTypes[i] && "Error: Element type mismatch.");
        assert(plist[i].size() == npes[i] && "Error: Incorrect number of nodes for element.");
        assert(pfverts[i].size() == epes[i] && "Error: Incorrect number of neighbors for element.");
    }
    for (size_t i = 0; i < n_nodes; ++i) {
        assert(px[i] >= 1.0 && px[i] <= 7.0 && py[i] >= 1.0 && py[i] <= 5.0 && pz[i] == 0.0 &&
               "Error: Node coordinates out of range.");
    }
    // Verify boundary flags
    assert(bflags[0] == CNR4 && bflags[4] == CNR1 && bflags[7] == CNR3 && bflags[11] == CNR2 &&
           "Error: Corner boundary flags incorrect.");
    assert(bflags[1] == TOP && bflags[6] == TOP && "Error: Top boundary flags incorrect.");
    assert(bflags[5] == BOTTOM && bflags[10] == BOTTOM && "Error: Bottom boundary flags incorrect.");
    assert(bflags[2] == LEFT && bflags[9] == RIGHT && "Error: Left/right boundary flags incorrect.");
    assert(bflags[3] == INTERNAL && bflags[8] == INTERNAL && "Error: Internal boundary flags incorrect.");

    cout << "create_Disconnected2D_VSet: model 'Disconnected2D': done." << endl;
}

/* FORMER VERSION
void create_Disconnected2D_VSet( VSet<2U>& vset )
 {
    //--------------------------ELEMENT TYPES
  	//add element types
    const CSMP_FEM_TYPE T(ISOPARAMETRIC_LINEAR_TRIANGLE), Q(ISOPARAMETRIC_LINEAR_QUADRILATERAL);
    deque<int8_t> vecElementTypes = { Q,Q,T,T,Q };
    const int n_cells{5};
  	assert( vecElementTypes.size() == n_cells );
    const size_t     n_nodes(12); // number of nodes
  	deque<uint32_t>  npes{ 4, 4, 3, 3, 4 };  // number of nodes per element
    deque<uint32_t>  epes{ 4, 4, 3, 3, 4 };  // nbors per element
    // 0 = no faces nor interfaces
  	vset.Resize( vecElementTypes, npes, epes, n_nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES (12)
  	//define node coordinates
  	deque<double> px = { 1., 4., 1., 4., 1., 4., 5., 7., 5., 7., 5., 7. };
    assert( px.size() == n_nodes );
  	deque<double> py = { 5., 5., 3., 3., 1., 1., 5., 5., 3., 3., 1., 1. };
    assert( py.size() == n_nodes );
  	deque<double> pz(n_nodes,0.);
  	  	  	
  	vset.AddXYZ( px, py, pz );

  	//--------------------------NODE BOUNDARY FLAGS (12)
    const BOX_BOUNDARY B{BOTTOM}, R{RIGHT}, U{TOP}, L{LEFT}, I{INTERNAL};
    vector<int8_t> bflags = { CNR4, U, L, I, CNR1, B, U, CNR3, I, R, B, CNR2 };
    assert( bflags.size() == n_nodes );
    
    vset.AddBFlags( bflags.begin(), bflags.end() );
    
    //--------------------------TOPOTYPE NODE FLAGS (12)
    const TOPOTYPE v{MESH_VERTEX}, e{EXTERIOR_POINT};
    //                        0  1  2  3  4  5  6  7  8  9 10 11
    vector<int8_t> gflags = { e, e, e, v, e, e, e, e, v, e, e, e };
    assert( gflags.size() == n_nodes );
    vset.AddBREP_Flags( gflags.begin(), gflags.end() );

  	//--------------------------ELEMENTS (5)
  	// define nodes per element
    deque< vector<size_t> >  plist = { {2,3,1,0}, {4,5,3,2}, {6,9,7}, {6,8,9}, {10,11,9,8} };
    assert( plist.size() == n_cells );
    vset.AddPlist( plist.begin(), plist.end() );

     //---------------------------------NEIGHBORS
    //define neighbors per element
    deque<vector<int64_t> > pfverts = { {1,I,U,L}, {B,I,0,L}, {4,2,I}, {R,U,3}, {B,R,3,I} };
    assert( pfverts.size() == n_cells );
    vset.AddPfverts( pfverts.begin(), pfverts.end() );

     //---------------------------------MATERIALS
    const size_t n_elements{5};
    vector<int32_t> pmtrl = { 1,1,2,2,2 };
    assert( pmtrl.size() == n_elements );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
     //---------------------------------ELEMENT & NODE NUMBERS
    PropertyData elmt_nums( ELEMENT, SCALAR, 2U );
    elmt_nums.Reserve( vset.Elements() );
    for ( size_t n{0U}; n<vset.Elements(); ++n ) pushBack( elmt_nums, makeScalar( ANY, n ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t n{0U}; n<vset.Vertices(); ++n ) pushBack( node_nums, makeScalar( ANY, n ) );
    vset.AddData( "node number", node_nums );
    
     //---------------------------------PERMEABILITY
    PropertyData perm( ELEMENT, SCALAR, 2U );
    perm.Reserve( vset.Elements() );
    for ( size_t n{0U}; n<5; ++n ) pushBack( perm, makeScalar( ANY, 1.0e-13 ) );
    vset.AddData( "permeability", perm );

     cout <<"\n"<<"create_Disconnected2D_VSet: model 'un-named': done."<< endl;
//    vset.Out();

 } // end create_Disconnected2D_VSet
*/





/// model  SPLIT22_BASIC  with box boundaries (Faces) and one through-going and one internal crossing split boundary
ModelTopology  create_BoundarySplitBoundaryPatch( VSet<2U>& vset )
 {
    //--------------------------ELEMENT TYPES (3)
  	//add element types
    const CSMP_FEM_TYPE T(ISOPARAMETRIC_LINEAR_TRIANGLE), Q(ISOPARAMETRIC_LINEAR_QUADRILATERAL), P(ISOPARAMETRIC_LINEAR_BAR);
    deque<int8_t> vecElementTypes = { Q,Q,Q,Q,T,Q,Q,Q,Q,T,Q,Q,Q,Q,T,Q,Q,Q,Q, // elements (16 quads + 3 triangles)
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
    const int8_t B{BOTTOM_OUTSIDE}, R{RIGHT_OUTSIDE}, U{TOP_OUTSIDE}, L{LEFT_OUTSIDE}, I{REGION_BOUNDARY}, N{0},
                 C1{CNR_MIN}, C2{CNR_X}, C3{CNR_XY}, C4{CNR_Y};
    //                          0 1 2  3  4 5 6 7  8 9 1011121314 151617181920 212223242526 27282930  31 32 33  34
    vector<int8_t> bflags = {  C4,U,U,C3, L,N,N,R, L,I, L,I,I,N,R, L,I,I,I,I,R, L,I,I,I,I,R, L,N,I,R, C1, B, B, C2 };
    assert( bflags.size() == n_nodes );
    vset.AddBFlags( bflags.begin(), bflags.end() );
    
    
    //--------------------------TOPOTYPE NODE FLAGS (35)
    const TOPOTYPE v{MESH_VERTEX}, i{INTERIOR_POINT}, e{EXTERIOR_POINT}, p{PERIMETER_POINT}, m{PERIMETER_LINE}, l{INTERIOR_LINE}, x{EXTERIOR_LINE};
    //                        0 1 2 3  4 5 6 7  8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
    vector<int8_t> gflags = { e,x,x,e, x,v,v,x, x,p, x, m, m, v, x, e, m, i, i, m, e, e, m, i, i, m, e, x, v, p, x, e, x, x, e };
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

     //---------------------------------NEIGHBORS                                                   was {5,10,7,4}
    //define neighbors per element          0          1          2          3         4         5           6           7          8        9
    deque<vector<int64_t> > pfverts = { {L,3,1,U}, {0,4,2,U}, {1,7,R,U}, {L,5,4,0}, {6,1,3}, {L,8,I,3}, {I,10,7,4}, {6,11,R,2}, {L,I,9,5}, {I,I,8}, // element neighbors
                                        {I,I,11,6}, {10,I,R,7}, {L,16,13,I}, {17,I,I,12}, {15,I,I}, {14,18,R,I}, {L,B,17,12}, {16,B,18,13}, {17,B,R,15},
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
    for ( size_t n{0U}; n<vset.Elements(); ++n ) pushBack( elmt_nums, makeScalar( ANY, static_cast<double>(n) ) );
    vset.AddData( "element number", elmt_nums );
    // node numbers
    PropertyData node_nums( NODE, SCALAR, 2U );
    node_nums.Reserve( vset.Vertices() );
    for ( size_t n{0U}; n<vset.Vertices(); ++n ) pushBack( node_nums, makeScalar( ANY, static_cast<double>(n) ) );
    vset.AddData( "node number", node_nums );
    // permeability
    PropertyData perm( ELEMENT, SCALAR, 2U );
    perm.Reserve( vset.Elements() );
    for ( size_t n{0U}; n<mesh_topology.CellsWithinDomain("upper"); ++n ) pushBack( perm, makeScalar( ANY, 1.0e-13 ) );
    for ( size_t n{0U}; n<mesh_topology.CellsWithinDomain("lower"); ++n ) pushBack( perm, makeScalar( ANY, 1.0e-12 ) );
    vset.AddData( "permeability", perm );

     cout <<"\n"<<"create_BoundarySplitBoundaryPatch: model 'SPLIT22_BASIC': done."<< endl;
//    vset.Out();
    
    return mesh_topology;
    
 } // end create_BoundarySplitBoundaryPatch






/**
        No neighbors, but (negative) boundary flags.
*/
void create_1Hexahedron_VSet(VSet<3U>& vset, bool bSkewed)
{
    IsoparametricLinearHexahedron iso_hexahedron;
    constexpr size_t nodes = 8;

    // Initialize VSet
    vset.SingleElementType(iso_hexahedron.ElementType());
    vset.Resize(iso_hexahedron.Nodes(), iso_hexahedron.Neighbors(),
                iso_hexahedron.ElementType(), nodes, 1);

    // Define node coordinates using initializer lists
    deque<double> px = {0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    deque<double> py = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    deque<double> pz = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};

    // Apply random perturbations if skewed
    if (bSkewed)
    {
        for (size_t i = 0; i < nodes; ++i)
        {
            px[i] += (rand() % 2000) * PERTURBATION;
            py[i] += (rand() % 2000) * PERTURBATION;
            pz[i] += (rand() % 2000) * PERTURBATION;
        }
    }

    // Load nodes
    vset.AddXYZ(px, py, pz);

    // Define element type
    vector<int8_t> vecElementTypes = {ISOPARAMETRIC_LINEAR_HEXAHEDRON};
    vset.AddElementTypes(vecElementTypes.begin(), vecElementTypes.end());

    // Define single hexahedron element
//    deque<vector<size_t>> deqElements = {{0, 1, 2, 3, 4, 5, 6, 7}};
    deque<vector<size_t>> deqElements = {{ 4, 5, 1, 0, 7, 6, 2, 3 }};
    vset.AddPlist(deqElements.begin(), deqElements.end());

    // Define neighbors (BOTTOM=0, FRONT=1, RIGHT=2, BACK=3, LEFT=4, TOP=5)
    deque<vector<int64_t>> deqElementNeighbors = {
        {BOTTOM_OUTSIDE, FRONT_OUTSIDE, RIGHT_OUTSIDE, BACK_OUTSIDE, LEFT_OUTSIDE, TOP_OUTSIDE}
    };
    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());

    // Assign node boundary flags
    vset.ResizeBFlags();
    vset.BFlag(0, CNR1);
    vset.BFlag(1, CNR2);
    vset.BFlag(2, CNR3);
    vset.BFlag(3, CNR4);
    vset.BFlag(4, CNR5);
    vset.BFlag(5, CNR6);
    vset.BFlag(6, CNR7);
    vset.BFlag(7, CNR8);

    // Assign material ID
    vector<int32_t> pmtrl(vset.Elements(), 1);
    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());

    cout << "\ncreate_1Hexahedron_VSet: model 'un-named': done.\n";
}




/**
    "Rubik-Cube"
    @test SKM 10/7/2024 - nbor connectivity and node flags are consistent with CSMP conventions
*/
void create_Hexahedra_VSet(VSet<3U>& vset, bool bSkewed)
{
    // 27 hexahedral elements on a 4x4x4 node grid (3x3x3 elements)
    const size_t iNrOfElements(27);
    
    IsoparametricLinearHexahedron iso_hexahedron;
    
    // Node dimensions of box (x right, y up, z forward)
    const size_t iDim_k(4); // k - height (z)
    const size_t iDim_j(4); // j - width (y)
    const size_t iDim_i(4); // i - length (x)
    const size_t iDim_k2(iDim_k*iDim_k); // k^2
    const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1)); // (k-1)^2
    
    // Total nodes: 64 (4x4x4)
    size_t nodes(iDim_i*iDim_j*iDim_k);

    //------------------------CREATE VSET
    vset.SingleElementType(iso_hexahedron.ElementType());
    vset.Resize(iso_hexahedron.Nodes(), // 8 nodes per hexahedron
                iso_hexahedron.Neighbors(), // 6 neighbors per hexahedron
                iso_hexahedron.ElementType(),
                nodes, iNrOfElements);
    
    //-----------------------NODES
    // Define node coordinates (0-based indexing)
    deque<double> px(nodes);
    deque<double> py(nodes);
    deque<double> pz(nodes);
  
    for(size_t k{0}; k < iDim_k; k++) // z (back to front)
        for(size_t j{0}; j < iDim_j; j++) // y (bottom to top)
            for(size_t i{0}; i < iDim_i; i++) // x (left to right)
            {
                size_t iNode = k*iDim_k2 + j*iDim_j + i;
                if(bSkewed)
                {
                    px[iNode] = i + (rand()%2000)*PERTURBATION;
                    py[iNode] = j + (rand()%2000)*PERTURBATION;
                    pz[iNode] = k + (rand()%2000)*PERTURBATION;
                }
                else
                {
                    px[iNode] = i;
                    py[iNode] = j;
                    pz[iNode] = k;
                }
            }
    
    // Load nodes
    vset.AddXYZ(px, py, pz);
    vset.ResizeBFlags();
    
    //--------------------------ELEMENTS
    // Define element types
    vector<int8_t> vecElementTypes(1, ISOPARAMETRIC_LINEAR_HEXAHEDRON);
    vset.AddElementTypes(vecElementTypes.begin(), vecElementTypes.end());

    // Define hexahedron elements (0 to 26), assign 8 nodes per element (0-based)
    deque<vector<size_t>> deqElements(iNrOfElements);
    for(size_t k{0}; k < iDim_k-1; k++) // z
        for(size_t j{0}; j < iDim_j-1; j++) // y
            for(size_t i{0}; i < iDim_i-1; i++) // x
            {
                const size_t iElement = iDim_km1_2*k + (iDim_j-1)*j + i;
                deqElements[iElement].resize(8);
                // Corrected node order to match CSMP face conventions
                deqElements[iElement][0] = iDim_k2*(k+1) + iDim_j*j + i;          // (i,j,k+1)
                deqElements[iElement][1] = iDim_k2*(k+1) + iDim_j*j + i + 1;      // (i+1,j,k+1)
                deqElements[iElement][2] = iDim_k2*k + iDim_j*j + i + 1;          // (i+1,j,k)
                deqElements[iElement][3] = iDim_k2*k + iDim_j*j + i;              // (i,j,k)
                deqElements[iElement][4] = iDim_k2*(k+1) + iDim_j*(j+1) + i;      // (i,j+1,k+1)
                deqElements[iElement][5] = iDim_k2*(k+1) + iDim_j*(j+1) + i + 1;  // (i+1,j+1,k+1)
                deqElements[iElement][6] = iDim_k2*k + iDim_j*(j+1) + i + 1;      // (i+1,j+1,k)
                deqElements[iElement][7] = iDim_k2*k + iDim_j*(j+1) + i;          // (i,j+1,k)
            }
    
    vset.AddPlist(deqElements.begin(), deqElements.end());

    //---------------------------------NEIGHBORS
    // Define neighbors (BOTTOM=0 (y=0), FRONT=1 (z=max), RIGHT=2 (x=max),
    //                  BACK=3 (z=0), LEFT=4 (x=0), TOP=5 (y=max))
    deque<vector<int64_t>> deqElementNeighbors(iNrOfElements);
    for(size_t k{0}; k < iDim_k-1; k++) // z
        for(size_t j{0}; j < iDim_j-1; j++) // y
            for(size_t i{0}; i < iDim_i-1; i++) // x
            {
                const size_t iElement = iDim_km1_2*k + (iDim_j-1)*j + i;
                deqElementNeighbors[iElement].resize(6);
                // Face 0: BOTTOM (y=0, xz-plane)
                deqElementNeighbors[iElement][0] = (j==0 ? BOTTOM_OUTSIDE : static_cast<int64_t>(iDim_km1_2*k + (iDim_j-1)*(j-1) + i));
                // Face 1: FRONT (z=max, xy-plane)
                deqElementNeighbors[iElement][1] = (k==iDim_k-2 ? FRONT_OUTSIDE : static_cast<int64_t>(iDim_km1_2*(k+1) + (iDim_j-1)*j + i));
                // Face 2: RIGHT (x=max, yz-plane)
                deqElementNeighbors[iElement][2] = (i==iDim_i-2 ? RIGHT_OUTSIDE : static_cast<int64_t>(iDim_km1_2*k + (iDim_j-1)*j + i + 1));
                // Face 3: BACK (z=0, xy-plane)
                deqElementNeighbors[iElement][3] = (k==0 ? BACK_OUTSIDE : static_cast<int64_t>(iDim_km1_2*(k-1) + (iDim_j-1)*j + i));
                // Face 4: LEFT (x=0, yz-plane)
                deqElementNeighbors[iElement][4] = (i==0 ? LEFT_OUTSIDE : static_cast<int64_t>(iDim_km1_2*k + (iDim_j-1)*j + i - 1));
                // Face 5: TOP (y=max, xz-plane)
                deqElementNeighbors[iElement][5] = (j==iDim_j-2 ? TOP_OUTSIDE : static_cast<int64_t>(iDim_km1_2*k + (iDim_j-1)*(j+1) + i));
            }
    
    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());
    
    //-----------------------------------------------------NODE BOUNDARIES
    // Define boundary flags for nodes
    for(size_t k{0}; k < iDim_k; k++) // z
        for(size_t j{0}; j < iDim_j; j++) // y
            for(size_t i{0}; i < iDim_i; i++) // x
            {
                int8_t bBoundary = NOT;
                if(k==0)
                {
                    if(j==0)
                    {
                        if(i==0) bBoundary = CNR1; // (0,0,0)
                        else if(i==iDim_i-1) bBoundary = CNR2; // (3,0,0)
                        else bBoundary = EDGE1; // x-axis, y=0, z=0
                    }
                    else if(j==iDim_j-1)
                    {
                        if(i==0) bBoundary = CNR4; // (0,3,0)
                        else if(i==iDim_i-1) bBoundary = CNR3; // (3,3,0)
                        else bBoundary = EDGE3; // x-axis, y=3, z=0
                    }
                    else
                    {
                        if(i==0) bBoundary = EDGE4; // y-axis, x=0, z=0
                        else if(i==iDim_i-1) bBoundary = EDGE2; // y-axis, x=3, z=0
                        else bBoundary = BACK_OUTSIDE; // z=0 plane
                    }
                }
                else if(k==iDim_k-1)
                {
                    if(j==0)
                    {
                        if(i==0) bBoundary = CNR5; // (0,0,3)
                        else if(i==iDim_i-1) bBoundary = CNR6; // (3,0,3)
                        else bBoundary = EDGE9; // x-axis, y=0, z=3
                    }
                    else if(j==iDim_j-1)
                    {
                        if(i==0) bBoundary = CNR8; // (0,3,3)
                        else if(i==iDim_i-1) bBoundary = CNR7; // (3,3,3)
                        else bBoundary = EDGE11; // x-axis, y=3, z=3
                    }
                    else
                    {
                        if(i==0) bBoundary = EDGE12; // y-axis, x=0, z=3
                        else if(i==iDim_i-1) bBoundary = EDGE10; // y-axis, x=3, z=3
                        else bBoundary = FRONT_OUTSIDE; // z=3 plane
                    }
                }
                else
                {
                    if(j==0)
                    {
                        if(i==0) bBoundary = EDGE5; // z-axis, x=0, y=0
                        else if(i==iDim_i-1) bBoundary = EDGE6; // z-axis, x=3, y=0
                        else bBoundary = BOTTOM_OUTSIDE; // y=0 plane
                    }
                    else if(j==iDim_j-1)
                    {
                        if(i==0) bBoundary = EDGE8; // z-axis, x=0, y=3
                        else if(i==iDim_i-1) bBoundary = EDGE7; // z-axis, x=3, y=3
                        else bBoundary = TOP_OUTSIDE; // y=3 plane
                    }
                    else
                    {
                        if(i==0) bBoundary = LEFT_OUTSIDE; // x=0 plane
                        else if(i==iDim_i-1) bBoundary = RIGHT_OUTSIDE; // x=3 plane
                        else ; // internal node, no boundary
                    }
                }
    
                if(bBoundary != NOT)
                {
                    const size_t iNode = iDim_k2*k + (iDim_j)*j + i;
                    vset.BFlag(iNode, bBoundary);
                }
            }
    
    //-------------------------MATERIALS
    // Assign material ID 1 to all elements
    vector<int32_t> pmtrl(vset.Elements(), 1);
    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());

    cout << "\n" << "create_Hexahedra_VSet: model 'un-named': done." << endl;
}







void create_Square_VSet( VSet<2U>& vset, size_t size_sides, double dimension, bool bSkewed )
{
  if(size_sides==1)
    create_1Square_VSet( vset, dimension, bSkewed );
  else
    create_SlitRectangle_VSet( vset, size_sides, size_sides, dimension, dimension, 0, bSkewed );
}




void create_SlitRectangle_VSet( VSet<2U>& vset, size_t x_dimension, size_t y_dimension,
                                double x_length, double y_length, size_t depth_of_slit, bool bSkewed )
{
    assert(x_dimension > 0);
    assert(y_dimension > 0);
    assert(depth_of_slit < x_dimension);
    assert(x_length > 1e-7);
    assert(y_length > 1e-7);

    IsoparametricLinearQuadrilateral iso_quadrilateral;

    // Node dimensions of box
    size_t iDim_i = x_dimension + 1; // x-direction nodes
    size_t iDim_j = y_dimension + 1; // y-direction nodes
    size_t iNrOfElements = (iDim_i - 1) * (iDim_j - 1);
    size_t nodes = iDim_i * iDim_j; // Total number of nodes initially

    cout << "create_SlitRectangle_VSet:\n";
    cout << "\tDim i: " << iDim_i << " Dim j: " << iDim_j
              << " nr of elements: " << iNrOfElements << " nodes: " << nodes << endl;

    vset.Resize(iso_quadrilateral.Nodes(),
                iso_quadrilateral.Neighbors(),
                iso_quadrilateral.ElementType(),
                nodes, iNrOfElements);

    // ----------------------- NODES
    deque<double> px(nodes);
    deque<double> py(nodes);
    deque<double> pz(nodes);

    double delta_x = x_length / static_cast<double>(x_dimension);
    double delta_y = y_length / static_cast<double>(y_dimension);

    for (size_t j = 0; j < iDim_j; ++j) { // y
        for (size_t i = 0; i < iDim_i; ++i) { // x
            size_t idx = iDim_i * j + i;
            if (bSkewed) {
                px[idx] = i * delta_x + (rand() % 2000) * PERTURBATION;
                py[idx] = j * delta_y + (rand() % 2000) * PERTURBATION;
                pz[idx] = 0.0;
            } else {
                px[idx] = i * delta_x;
                py[idx] = j * delta_y;
                pz[idx] = 0.0;
            }
        }
    }

    // ----------------------- ELEMENTS
    deque<vector<size_t>> deqElements(iNrOfElements);
    for (size_t j = 0; j < iDim_j - 1; ++j) { // y
        for (size_t i = 0; i < iDim_i - 1; ++i) { // x
            size_t iElement = (iDim_i - 1) * j + i;
            deqElements[iElement].resize(4);
            deqElements[iElement][0] = iDim_i * j + i;         // Bottom-left
            deqElements[iElement][1] = iDim_i * j + i + 1;     // Bottom-right
            deqElements[iElement][2] = iDim_i * (j + 1) + i + 1; // Top-right
            deqElements[iElement][3] = iDim_i * (j + 1) + i;     // Top-left
        }
    }

    // ----------------------- NEIGHBORS
    deque<vector<int64_t>> deqElementNeighbors(iNrOfElements);
    for (size_t j = 0; j < iDim_j - 1; ++j) { // y
        for (size_t i = 0; i < iDim_i - 1; ++i) { // x
            size_t iElement = (iDim_i - 1) * j + i;
            deqElementNeighbors[iElement].resize(4);

            // Face 0: BOTTOM (connects to element below, j-1)
            deqElementNeighbors[iElement][0] = (j == 0) ? BOTTOM_OUTSIDE : static_cast<int64_t>((iDim_i - 1) * (j - 1) + i);
            // Face 2: RIGHT (connects to element right, i+1)
            deqElementNeighbors[iElement][1] = (i == iDim_i - 2) ? RIGHT_OUTSIDE : static_cast<int64_t>((iDim_i - 1) * j + i + 1);
            // Face 4: TOP (connects to element above, j+1)
            deqElementNeighbors[iElement][2] = (j == iDim_j - 2) ? TOP_OUTSIDE : static_cast<int64_t>((iDim_i - 1) * (j + 1) + i);
            // Face 6: LEFT (connects to element left, i-1)
            deqElementNeighbors[iElement][3] = (i == 0) ? LEFT_OUTSIDE : static_cast<int64_t>((iDim_i - 1) * j + i - 1);

            bool over_slit = (j == y_dimension / 2 && i >= x_dimension - depth_of_slit);
            bool under_slit = (j == y_dimension / 2 - 1 && i >= x_dimension - depth_of_slit);

            if (over_slit) {
                deqElementNeighbors[iElement][1] = static_cast<int64_t>(IRREGULAR_OUTSIDE); // Right face on slit
            } else if (under_slit) {
                deqElementNeighbors[iElement][2] = static_cast<int64_t>(IRREGULAR_OUTSIDE); // Top face on slit
            }
        }
    }

    // ----------------------- NODE BOUNDARIES
    vector<int8_t> bflags(nodes, NOT);
    for (size_t j = 0; j < iDim_j; ++j) { // y
        for (size_t i = 0; i < iDim_i; ++i) { // x
            size_t iNode = iDim_i * j + i;
            int8_t& bBoundary = bflags[iNode];

            // Corners
            if (i == 0 && j == 0) {
                bBoundary = CNR1; // (0,0,0)
            } else if (i == iDim_i - 1 && j == 0) {
                bBoundary = CNR2; // (max,0,0)
            } else if (i == iDim_i - 1 && j == iDim_j - 1) {
                bBoundary = CNR3; // (max,max,0)
            } else if (i == 0 && j == iDim_j - 1) {
                bBoundary = CNR4; // (0,max,0)
            }
            // Edges (excluding corners)
            else if (j == 0) {
                bBoundary = BOTTOM_OUTSIDE; // y=0
            } else if (j == iDim_j - 1) {
                bBoundary = TOP_OUTSIDE; // y=max
            } else if (i == 0) {
                bBoundary = LEFT_OUTSIDE; // x=0
            } else if (i == iDim_i - 1) {
                bBoundary = RIGHT_OUTSIDE; // x=max
            }
        }
    }

    // ----------------------- SLIT NODES
    cout << "\tIntroduce Slit Nodes..." << endl;
    for (size_t j = 0; j < iDim_j - 1; ++j) { // y
        for (size_t i = 0; i < iDim_i - 1; ++i) { // x
            bool over_slit = (j == y_dimension / 2 && i >= x_dimension - depth_of_slit);
            bool under_slit = (j == y_dimension / 2 - 1 && i >= x_dimension - depth_of_slit);
            size_t iElement = (iDim_i - 1) * j + i;

            if (over_slit) {
                if (i > x_dimension - depth_of_slit) { // Node 0 (bottom-left)
                    size_t node_number = px.size();
                    px.push_back(px[iDim_i * j + i]);
                    py.push_back(py[iDim_i * j + i]);
                    pz.push_back(pz[iDim_i * j + i]);
                    deqElements[iElement][0] = node_number;
                    bflags.push_back(IRREGULAR_OUTSIDE);
                }
                // Node 1 (bottom-right)
                size_t node_number = px.size();
                px.push_back(px[iDim_i * j + i + 1]);
                py.push_back(py[iDim_i * j + i + 1]);
                pz.push_back(pz[iDim_i * j + i + 1]);
                deqElements[iElement][1] = node_number;
                bflags.push_back(IRREGULAR_OUTSIDE);
            }

            if (under_slit) {
                // Node 2 (top-right)
                size_t node_number = px.size();
                px.push_back(px[iDim_i * (j + 1) + i + 1]);
                py.push_back(py[iDim_i * (j + 1) + i + 1]);
                pz.push_back(pz[iDim_i * (j + 1) + i + 1]);
                deqElements[iElement][2] = node_number;
                bflags.push_back(IRREGULAR_OUTSIDE);

                if (i > x_dimension - depth_of_slit) { // Node 3 (top-left)
                    node_number = px.size();
                    px.push_back(px[iDim_i * (j + 1) + i]);
                    py.push_back(py[iDim_i * (j + 1) + i]);
                    pz.push_back(pz[iDim_i * (j + 1) + i]);
                    deqElements[iElement][3] = node_number;
                    bflags.push_back(IRREGULAR_OUTSIDE);
                }
            }
        }
    }

    // Load nodes and boundary flags
    vset.AddXYZ(px, py, pz);
    for (size_t i = 0; i < bflags.size(); ++i) {
        if (bflags[i] != NOT) {
            vset.BFlag(i, bflags[i]);
        }
    }
    vset.AddPlist(deqElements.begin(), deqElements.end());
    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());

    // ----------------------- MATERIALS
    vector<int32_t> pmtrl(vset.Elements(), 1);
    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());

    cout << "splitRectangle_VSet: model 'un-named': done." << endl;
}







/**
    3D model, which is a cube of hexahedra with six pyramid elements in the middle.
    
    Creates:
    - 32 elements (26 hex + 6 pyramids)
    - 64 nodes
    
    @test boundary flags not correct yet
*/
/* FIXED create_Pyramid_Hexa_VSet() as follows: */
void create_Pyramid_Hexa_VSet(VSet<3U>& vset, bool bSkewed)
{
    static constexpr int32_t iNrOfElements = 32; // 26 hexahedrons + 6 pyramids
    constexpr double EPSILON = 1e-6; // Tolerance for boundary checks

    IsoparametricLinearHexahedron iso_hexahedron;
    IsoparametricLinearPyramid iso_pyramid;

    // Node dimensions of box
    static constexpr size_t iDim_i = 4; // x - length
    static constexpr size_t iDim_j = 4; // y - width
    static constexpr size_t iDim_k = 4; // z - height
    static constexpr size_t iDim_k2 = iDim_k * iDim_j; // k * j
    static constexpr size_t iDim_km1_2 = (iDim_k - 1) * (iDim_j - 1); // (k-1) * (j-1)
    static constexpr size_t iPyramidsPlacement = 13; // Center element index

    // Total nodes: 4x4x4 grid (64) + 1 barycenter
    static constexpr size_t nodes = iDim_i * iDim_j * iDim_k + 1;

    // Initialize element properties
    std::deque<uint32_t> npes(iNrOfElements);
    std::deque<uint32_t> epes(iNrOfElements);
    std::deque<int8_t> etypes(iNrOfElements);

    for (size_t iElement = 0; iElement < 26; ++iElement) {
        npes[iElement] = iso_hexahedron.Nodes(); // 8 nodes
        epes[iElement] = iso_hexahedron.Neighbors(); // 6 neighbors
        etypes[iElement] = iso_hexahedron.ElementType();
    }
    for (size_t iElement = 26; iElement < iNrOfElements; ++iElement) {
        npes[iElement] = iso_pyramid.Nodes(); // 5 nodes
        epes[iElement] = iso_pyramid.Neighbors(); // 5 neighbors
        etypes[iElement] = iso_pyramid.ElementType();
    }

    vset.Resize(etypes, npes, epes, nodes, 0, 0);
    vset.AddElementTypes(etypes.begin(), etypes.end());

    // ----------------------- MATERIALS
    static constexpr int32_t material_id = 5;
    std::vector<int32_t> materials(iNrOfElements, material_id);
    vset.AddPmtrl(materials.begin(), materials.end());

    // ----------------------- NODES
    std::deque<double> px(nodes);
    std::deque<double> py(nodes);
    std::deque<double> pz(nodes);

    // Initialize random number generator for perturbations
    std::mt19937 rng(12345); // Fixed seed for reproducibility
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t k = 0; k < iDim_k; ++k) { // z
        for (size_t j = 0; j < iDim_j; ++j) { // y
            for (size_t i = 0; i < iDim_i; ++i) { // x
                size_t idx = k * iDim_k2 + j * iDim_j + i;
                if (bSkewed) {
                    px[idx] = static_cast<double>(i) + dist(rng) * PERTURBATION;
                    py[idx] = static_cast<double>(j) + dist(rng) * PERTURBATION;
                    pz[idx] = static_cast<double>(k) + dist(rng) * PERTURBATION;
                } else {
                    px[idx] = static_cast<double>(i);
                    py[idx] = static_cast<double>(j);
                    pz[idx] = static_cast<double>(k);
                }
            }
        }
    }

    // Barycenter node (index 64)
    px[64] = 1.5 + (bSkewed ? dist(rng) * PERTURBATION : 0.0);
    py[64] = 1.5 + (bSkewed ? dist(rng) * PERTURBATION : 0.0);
    pz[64] = 1.5 + (bSkewed ? dist(rng) * PERTURBATION : 0.0);

    vset.AddXYZ(px, py, pz);
    vset.ResizeBFlags();

    // ----------------------- ELEMENTS
    std::deque<std::vector<size_t>> deqElements(iNrOfElements);
    for (size_t k = 0; k < iDim_k - 1; ++k) { // z
        for (size_t j = 0; j < iDim_j - 1; ++j) { // y
            for (size_t i = 0; i < iDim_i - 1; ++i) { // x
                size_t iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i;
                if (i == 1 && j == 1 && k == 1) continue; // Skip center for pyramids
                if (iElement > iPyramidsPlacement) --iElement;

                deqElements[iElement].resize(8);
                deqElements[iElement][0] = k * iDim_k2 + j * iDim_j + i;         // Bottom-left-back (i,j,k)
                deqElements[iElement][1] = k * iDim_k2 + j * iDim_j + i + 1;     // Bottom-right-back (i+1,j,k)
                deqElements[iElement][2] = k * iDim_k2 + (j + 1) * iDim_j + i + 1; // Bottom-right-front (i+1,j+1,k)
                deqElements[iElement][3] = k * iDim_k2 + (j + 1) * iDim_j + i;     // Bottom-left-front (i,j+1,k)
                deqElements[iElement][4] = (k + 1) * iDim_k2 + j * iDim_j + i;     // Top-left-back (i,j,k+1)
                deqElements[iElement][5] = (k + 1) * iDim_k2 + j * iDim_j + i + 1; // Top-right-back (i+1,j,k+1)
                deqElements[iElement][6] = (k + 1) * iDim_k2 + (j + 1) * iDim_j + i + 1; // Top-right-front (i+1,j+1,k+1)
                deqElements[iElement][7] = (k + 1) * iDim_k2 + (j + 1) * iDim_j + i;     // Top-left-front (i,j+1,k+1)
            }
        }
    }

    // Pyramid elements (26–31) with adjusted node orders for consistent outward normals
    // Pyramid 26
    deqElements[26] = {25, 21, 22, 26, 64};
    // Pyramid 27
    deqElements[27] = {38, 22, 21, 37, 64};
    // Pyramid 28
    deqElements[28] = {26, 22, 38, 42, 64};
    // Pyramid 29
    deqElements[29] = {41, 25, 26, 42, 64};
    // Pyramid 30
    deqElements[30] = {37, 21, 25, 41, 64};
    // Pyramid 31
    deqElements[31] = {38, 37, 41, 42, 64};

    vset.AddPlist(deqElements.begin(), deqElements.end());

    // ---------------------------------NEIGHBORS (Algorithmic Detection)
    std::deque<std::vector<int64_t>> deqElementNeighbors(iNrOfElements);

    // Define face nodes for each element type
    const std::vector<std::vector<size_t>> hexaFaceNodes = {
        {0, 3, 2, 1}, // Face 0
        {0, 1, 5, 4}, // Face 1
        {1, 2, 6, 5}, // Face 2
        {2, 3, 7, 6}, // Face 3
        {0, 4, 7, 3}, // Face 4
        {4, 5, 6, 7}  // Face 5
    };
    const std::vector<std::vector<size_t>> pyramidFaceNodes = {
        {0, 1, 4},    // Face 0
        {1, 2, 4},    // Face 1
        {2, 3, 4},    // Face 2
        {0, 4, 3},    // Face 3
        {0, 3, 2, 1}  // Face 4
    };

    // Helper function to get sorted face nodes (for matching)
    auto getSortedFaceNodes = [](const std::vector<size_t>& elementNodes, const std::vector<size_t>& faceIndices) {
        std::vector<size_t> faceNodes;
        for (size_t idx : faceIndices) {
            faceNodes.push_back(elementNodes[idx]);
        }
        std::sort(faceNodes.begin(), faceNodes.end());
        return faceNodes;
    };

    // Helper function to get face nodes in original order
    auto getFaceNodes = [](const std::vector<size_t>& elementNodes, const std::vector<size_t>& faceIndices) {
        std::vector<size_t> faceNodes;
        for (size_t idx : faceIndices) {
            faceNodes.push_back(elementNodes[idx]);
        }
        return faceNodes;
    };

    // Initialize neighbor vectors
    for (size_t e = 0; e < iNrOfElements; ++e) {
        deqElementNeighbors[e].resize(epes[e], NOT);
    }

    // Compute neighbors
    for (size_t e = 0; e < iNrOfElements; ++e) {
        bool isPyramid = (etypes[e] == iso_pyramid.ElementType());
        const auto& faceNodesRef = isPyramid ? pyramidFaceNodes : hexaFaceNodes;
        size_t numFaces = isPyramid ? 5 : 6;

        for (size_t f = 0; f < numFaces; ++f) {
            // Get sorted face nodes for current element's face
            std::vector<size_t> faceNodes = getSortedFaceNodes(deqElements[e], faceNodesRef[f]);
            std::vector<size_t> faceNodesOrdered = getFaceNodes(deqElements[e], faceNodesRef[f]);

            // Search for matching face in other elements
            bool found = false;
            for (size_t e2 = 0; e2 < iNrOfElements; ++e2) {
                if (e2 == e) continue;
                bool isPyramid2 = (etypes[e2] == iso_pyramid.ElementType());
                const auto& faceNodesRef2 = isPyramid2 ? pyramidFaceNodes : hexaFaceNodes;
                size_t numFaces2 = isPyramid2 ? 5 : 6;

                for (size_t f2 = 0; f2 < numFaces2; ++f2) {
                    std::vector<size_t> faceNodes2 = getSortedFaceNodes(deqElements[e2], faceNodesRef2[f2]);
                    if (faceNodes == faceNodes2) {
                        // Check if nodes match in reverse order (outward normals), accounting for cyclic permutations
                        std::vector<size_t> faceNodes2Ordered = getFaceNodes(deqElements[e2], faceNodesRef2[f2]);
                        std::vector<size_t> faceNodesRev = faceNodesOrdered;
                        std::reverse(faceNodesRev.begin(), faceNodesRev.end());

                        // Double the second ordered list for cyclic check
                        std::vector<size_t> doubled = faceNodes2Ordered;
                        doubled.insert(doubled.end(), faceNodes2Ordered.begin(), faceNodes2Ordered.end());

                        bool opposite = false;
                        size_t len = faceNodesRev.size();
                        for (size_t s = 0; s < faceNodes2Ordered.size(); ++s) {
                            if (std::equal(faceNodesRev.begin(), faceNodesRev.end(), doubled.begin() + s)) {
                                opposite = true;
                                break;
                            }
                        }

                        if (opposite) {
                            deqElementNeighbors[e][f] = static_cast<int64_t>(e2);
                            found = true;
                            break;
                        }
                    }
                }
                if (found) break;
            }

            // If no neighbor found, assign boundary flag
            if (!found) {
                // Compute face centroid to determine boundary
                double cx = 0, cy = 0, cz = 0;
                for (size_t node : faceNodesOrdered) {
                    cx += px[node];
                    cy += py[node];
                    cz += pz[node];
                }
                cx /= faceNodesOrdered.size();
                cy /= faceNodesOrdered.size();
                cz /= faceNodesOrdered.size();

                if (std::abs(cy) < EPSILON) deqElementNeighbors[e][f] = BOTTOM_OUTSIDE;
                else if (std::abs(cy - (iDim_j - 1)) < EPSILON) deqElementNeighbors[e][f] = TOP_OUTSIDE;
                else if (std::abs(cx) < EPSILON) deqElementNeighbors[e][f] = LEFT_OUTSIDE;
                else if (std::abs(cx - (iDim_i - 1)) < EPSILON) deqElementNeighbors[e][f] = RIGHT_OUTSIDE;
                else if (std::abs(cz) < EPSILON) deqElementNeighbors[e][f] = BACK_OUTSIDE;
                else if (std::abs(cz - (iDim_k - 1)) < EPSILON) deqElementNeighbors[e][f] = FRONT_OUTSIDE;
                else deqElementNeighbors[e][f] = INTERNAL; // Internal boundary
            }
        }
    }

    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());

    // ----------------------- NODE BOUNDARY FLAGS
    std::vector<int8_t> bflags(nodes, NOT);
    for (size_t k = 0; k < iDim_k; ++k) { // z
        for (size_t j = 0; j < iDim_j; ++j) { // y
            for (size_t i = 0; i < iDim_i; ++i) { // x
                size_t iNode = k * iDim_k2 + j * iDim_j + i;
                int8_t bBoundary = NOT;

                if (bSkewed) {
                    if (std::abs(px[iNode]) < EPSILON) {
                        if (std::abs(py[iNode]) < EPSILON && std::abs(pz[iNode]) < EPSILON) bBoundary = CNR1;
                        else if (std::abs(py[iNode]) < EPSILON && std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR5;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON && std::abs(pz[iNode]) < EPSILON) bBoundary = CNR4;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON && std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR8;
                        else if (std::abs(py[iNode]) < EPSILON) bBoundary = EDGE5;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON) bBoundary = EDGE8;
                        else if (std::abs(pz[iNode]) < EPSILON) bBoundary = EDGE4;
                        else if (std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = EDGE12;
                        else bBoundary = LEFT_OUTSIDE;
                    } else if (std::abs(px[iNode] - (iDim_i - 1)) < EPSILON) {
                        if (std::abs(py[iNode]) < EPSILON && std::abs(pz[iNode]) < EPSILON) bBoundary = CNR2;
                        else if (std::abs(py[iNode]) < EPSILON && std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR6;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON && std::abs(pz[iNode]) < EPSILON) bBoundary = CNR3;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON && std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR7;
                        else if (std::abs(py[iNode]) < EPSILON) bBoundary = EDGE6;
                        else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON) bBoundary = EDGE7;
                        else if (std::abs(pz[iNode]) < EPSILON) bBoundary = EDGE2;
                        else if (std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = EDGE10;
                        else bBoundary = RIGHT_OUTSIDE;
                    } else if (std::abs(py[iNode]) < EPSILON) {
                        if (std::abs(pz[iNode]) < EPSILON) bBoundary = EDGE1;
                        else if (std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = EDGE9;
                        else bBoundary = BOTTOM_OUTSIDE;
                    } else if (std::abs(py[iNode] - (iDim_j - 1)) < EPSILON) {
                        if (std::abs(pz[iNode]) < EPSILON) bBoundary = EDGE3;
                        else if (std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = EDGE11;
                        else bBoundary = TOP_OUTSIDE;
                    } else if (std::abs(pz[iNode]) < EPSILON) {
                        bBoundary = BACK_OUTSIDE;
                    } else if (std::abs(pz[iNode] - (iDim_k - 1)) < EPSILON) {
                        bBoundary = FRONT_OUTSIDE;
                    }
                } else {
                    // Back (z=0)
                    if (k == 0) {
                        if (j == 0) { // Edge 1
                            if (i == 0) bBoundary = CNR1;
                            else if (i == iDim_i - 1) bBoundary = CNR2;
                            else bBoundary = EDGE1;
                        } else if (j == iDim_j - 1) { // Edge 3
                            if (i == 0) bBoundary = CNR4;
                            else if (i == iDim_i - 1) bBoundary = CNR3;
                            else bBoundary = EDGE3;
                        } else { // j in middle
                            if (i == 0) bBoundary = EDGE4;
                            else if (i == iDim_i - 1) bBoundary = EDGE2;
                            else bBoundary = BACK_OUTSIDE;
                        }
                    }
                    // Front (z=max)
                    else if (k == iDim_k - 1) {
                        if (j == 0) { // Edge 9
                            if (i == 0) bBoundary = CNR5;
                            else if (i == iDim_i - 1) bBoundary = CNR6;
                            else bBoundary = EDGE9;
                        } else if (j == iDim_j - 1) { // Edge 11
                            if (i == 0) bBoundary = CNR8;
                            else if (i == iDim_i - 1) bBoundary = CNR7;
                            else bBoundary = EDGE11;
                        } else { // j in middle
                            if (i == 0) bBoundary = EDGE12;
                            else if (i == iDim_i - 1) bBoundary = EDGE10;
                            else bBoundary = FRONT_OUTSIDE;
                        }
                    }
                    // Middle k
                    else {
                        if (j == 0) { // Bottom (y=0)
                            if (i == 0) bBoundary = EDGE5;
                            else if (i == iDim_i - 1) bBoundary = EDGE6;
                            else bBoundary = BOTTOM_OUTSIDE;
                        } else if (j == iDim_j - 1) { // Top (y=max)
                            if (i == 0) bBoundary = EDGE8;
                            else if (i == iDim_i - 1) bBoundary = EDGE7;
                            else bBoundary = TOP_OUTSIDE;
                        } else { // j in middle
                            if (i == 0) bBoundary = LEFT_OUTSIDE;
                            else if (i == iDim_i - 1) bBoundary = RIGHT_OUTSIDE;
                        }
                    }
                }
                bflags[iNode] = bBoundary;
            }
        }
    }
    // Barycenter node is internal
    bflags[64] = NOT;

    // Apply boundary flags
    for (size_t i = 0; i < bflags.size(); ++i) {
        if (bflags[i] != NOT) {
            vset.BFlag(i, bflags[i]);
        }
    }

    // ----------------------- CONSISTENCY CHECKS
    for (size_t i = 0; i < 26; ++i) {
        assert(etypes[i] == iso_hexahedron.ElementType() && "Error: Element type mismatch for hexahedrons.");
        assert(npes[i] == 8 && "Error: Hexahedron does not have 8 nodes.");
    }
    for (size_t i = 26; i < iNrOfElements; ++i) {
        assert(etypes[i] == iso_pyramid.ElementType() && "Error: Element type mismatch for pyramids.");
        assert(npes[i] == 5 && "Error: Pyramid does not have 5 nodes.");
    }
    for (size_t k = 0; k < iDim_k; ++k) {
        for (size_t j = 0; j < iDim_j; ++j) {
            for (size_t i = 0; i < iDim_i; ++i) {
                size_t idx = k * iDim_k2 + j * iDim_j + i;
                if (bSkewed) {
                    assert(px[idx] >= static_cast<double>(i) && px[idx] <= static_cast<double>(i) + PERTURBATION &&
                           "Error: Skewed px coordinate out of range.");
                    assert(py[idx] >= static_cast<double>(j) && py[idx] <= static_cast<double>(j) + PERTURBATION &&
                           "Error: Skewed py coordinate out of range.");
                    assert(pz[idx] >= static_cast<double>(k) && pz[idx] <= static_cast<double>(k) + PERTURBATION &&
                           "Error: Skewed pz coordinate out of range.");
                } else {
                    assert(px[idx] == static_cast<double>(i) && py[idx] == static_cast<double>(j) &&
                           pz[idx] == static_cast<double>(k) && "Error: Non-skewed coordinate mismatch.");
                }
            }
        }
    }
    // Barycenter check
    assert(px[64] >= 1.5 && px[64] <= 1.5 + PERTURBATION && "Error: Barycenter px out of expected range.");
    assert(py[64] >= 1.5 && py[64] <= 1.5 + PERTURBATION && "Error: Barycenter py out of expected range.");
    assert(pz[64] >= 1.5 && pz[64] <= 1.5 + PERTURBATION && "Error: Barycenter pz out of expected range.");

    std::cout << "create_Pyramid_Hexa_VSet: model 'un-named': done." << std::endl;
}









/**
    Generates 24 hexahedra + 6 prism elements.
    The model can be distorted on demand.
    
    @note model comes with the correct box boundary flags.
    
    TODO: boundary labelling for BACK and FRONT are not correct
*/
void create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed) {
    constexpr size_t iNrOfElements = 30;
    constexpr double EPSILON = PERTURBATION * 0.5 + 1e-6;

    IsoparametricLinearHexahedron iso_hexahedron;
    IsoparametricLinearPrism iso_prism;

    constexpr size_t iDim_i = 4;
    constexpr size_t iDim_j = 4;
    constexpr size_t iDim_k = 4;
    constexpr size_t iDim_k2 = iDim_i * iDim_j;
    const size_t nodes = iDim_i * iDim_j * iDim_k;

    deque<uint32_t> npes(iNrOfElements);
    deque<uint32_t> epes(iNrOfElements);
    deque<int8_t> etypes(iNrOfElements);

    size_t iElement = 0;
    for (size_t k = 0; k < iDim_k - 1; ++k) {
        for (size_t j = 0; j < iDim_j - 1; ++j) {
            for (size_t i = 0; i < iDim_i - 1; ++i) {
                if (i == 1 && j == 1) {
                    npes[iElement] = iso_prism.Nodes();
                    epes[iElement] = iso_prism.Neighbors();
                    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
                    ++iElement;

                    npes[iElement] = iso_prism.Nodes();
                    epes[iElement] = iso_prism.Neighbors();
                    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
                    ++iElement;
                    continue;
                }
                npes[iElement] = iso_hexahedron.Nodes();
                epes[iElement] = iso_hexahedron.Neighbors();
                etypes[iElement] = ISOPARAMETRIC_LINEAR_HEXAHEDRON;
                ++iElement;
            }
        }
    }
    vset.Resize(etypes, npes, epes, nodes, 0, 0);
    vset.AddElementTypes(etypes.begin(), etypes.end());

    // Nodes
    deque<double> px(nodes), py(nodes), pz(nodes);
    mt19937 rng(12345);
    uniform_real_distribution<double> dist(0.0, 1.0);
    for (size_t k = 0; k < iDim_k; ++k) {
        for (size_t j = 0; j < iDim_j; ++j) {
            for (size_t i = 0; i < iDim_i; ++i) {
                const size_t iNode = k * iDim_k2 + iDim_j * j + i;
                px[iNode] = i + (bSkewed ? dist(rng) * PERTURBATION : 0.0);
                py[iNode] = j + (bSkewed ? dist(rng) * PERTURBATION : 0.0);
                pz[iNode] = k + (bSkewed ? dist(rng) * PERTURBATION : 0.0);
            }
        }
    }
    vset.AddXYZ(px, py, pz);
    vset.ResizeBFlags();

    // Elements
    deque<vector<size_t>> deqElements(iNrOfElements);
    iElement = 0;
    for (size_t k = 0; k < iDim_k - 1; ++k) {         // z - axis
        for (size_t j = 0; j < iDim_j - 1; ++j) {     // y
            for (size_t i = 0; i < iDim_i - 1; ++i) { // x
                if (i == 1 && j == 1) {
                    // nodes of a cubic cell
                    const size_t node1 = k * iDim_k2 + iDim_j * j + i;
                    const size_t node2 = k * iDim_k2 + iDim_j * j + i + 1;
                    const size_t node3 = k * iDim_k2 + iDim_j * (j + 1) + i + 1;
                    const size_t node4 = k * iDim_k2 + iDim_j * (j + 1) + i;
                    const size_t node5 = (k + 1) * iDim_k2 + iDim_j * j + i;
                    const size_t node6 = (k + 1) * iDim_k2 + iDim_j * j + i + 1;
                    const size_t node7 = (k + 1) * iDim_k2 + iDim_j * (j + 1) + i + 1;
                    const size_t node8 = (k + 1) * iDim_k2 + iDim_j * (j + 1) + i;

                    // prism elements (wrong face numbering)
                    // - SKM corrected and tested 5/12/25 (clockwise because bottom and top of prism are flipped)
                    // bottom-face z is < top-face z
                    deqElements[iElement] = { node1, node2, node4, node5, node6, node8 };
                    ++iElement;

//                    deqElements[iElement] = {node6, node7, node8, node2, node3, node4};
                    // bottom-face z is < top-face z
                    deqElements[iElement] = { node2, node3, node4, node6, node7, node8 };
                    ++iElement;
                    continue;
                }

                // hexahedra
                deqElements[iElement] = {
                    k * iDim_k2 + iDim_j * j + i,
                    k * iDim_k2 + iDim_j * j + i + 1,
                    k * iDim_k2 + iDim_j * (j + 1) + i + 1,
                    k * iDim_k2 + iDim_j * (j + 1) + i,
                    (k + 1) * iDim_k2 + iDim_j * j + i,
                    (k + 1) * iDim_k2 + iDim_j * j + i + 1,
                    (k + 1) * iDim_k2 + iDim_j * (j + 1) + i + 1,
                    (k + 1) * iDim_k2 + iDim_j * (j + 1) + i
                };
                ++iElement;
            }
        }
    }
    vset.AddPlist(deqElements.begin(), deqElements.end());


    // Neighbor connectivity
    // =====================
    deque<vector<int64_t>> deqElementNeighbors(iNrOfElements);

    const vector<vector<size_t>> hexaFaceNodes = {
        {0, 3, 2, 1}, // bottom - counter-clockwise numbering
        {0, 1, 5, 4}, // front
        {1, 2, 6, 5}, // right
        {2, 3, 7, 6}, // back
        {0, 4, 7, 3}, // left
        {4, 5, 6, 7}  // top
    };
    const vector<vector<size_t>> prismFaceNodes = {
        {0, 2, 1},    // Face 0: bottom (triangle, adjusted for -z normal)
        {0, 1, 4, 3}, // Face 1: quad
        {1, 2, 5, 4}, // Face 2: quad
        {0, 3, 5, 2}, // Face 3: quad
        {3, 4, 5}     // Face 4: top (triangle, adjusted for +z normal)
    };

    auto getSortedFaceNodes = [](const vector<size_t>& elementNodes, const vector<size_t>& faceIndices) {
        vector<size_t> faceNodes;
        for (size_t idx : faceIndices) {
            faceNodes.push_back(elementNodes[idx]);
        }
        sort(faceNodes.begin(), faceNodes.end());
        return faceNodes;
    };

    auto getFaceNodes = [](const vector<size_t>& elementNodes, const vector<size_t>& faceIndices) {
        vector<size_t> faceNodes;
        for (size_t idx : faceIndices) {
            faceNodes.push_back(elementNodes[idx]);
        }
        return faceNodes;
    };

    auto computeFaceNormal = [&](const vector<size_t>& faceNodesOrdered) -> vector<double> {
        vector<double> normal(3, 0.0);
        size_t n = faceNodesOrdered.size();
        if (n < 3) return normal;
        for (size_t i = 0; i < n; ++i) {
            size_t j = (i + 1) % n;
            normal[0] += (py[faceNodesOrdered[i]] - py[faceNodesOrdered[j]]) * (pz[faceNodesOrdered[i]] + pz[faceNodesOrdered[j]]);
            normal[1] += (pz[faceNodesOrdered[i]] - pz[faceNodesOrdered[j]]) * (px[faceNodesOrdered[i]] + px[faceNodesOrdered[j]]);
            normal[2] += (px[faceNodesOrdered[i]] - px[faceNodesOrdered[j]]) * (py[faceNodesOrdered[i]] + py[faceNodesOrdered[j]]);
        }
        double mag = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        if (mag > EPSILON) {
            normal[0] /= mag;
            normal[1] /= mag;
            normal[2] /= mag;
        }
        return normal;
    };

    for (size_t e = 0; e < iNrOfElements; ++e) {
        deqElementNeighbors[e].resize(epes[e], NOT);
    }

    for (size_t e = 0; e < iNrOfElements; ++e) {
        bool isPrism = (etypes[e] == ISOPARAMETRIC_LINEAR_PRISM);
        const auto& faceNodesRef = isPrism ? prismFaceNodes : hexaFaceNodes;
        size_t numFaces = isPrism ? 5 : 6;

        for (size_t f = 0; f < numFaces; ++f) {
            vector<size_t> faceNodesOrdered = getFaceNodes(deqElements[e], faceNodesRef[f]);
            vector<size_t> faceNodes = getSortedFaceNodes(deqElements[e], faceNodesRef[f]);

            bool found = false;
            for (size_t e2 = 0; e2 < iNrOfElements; ++e2) {
                if (e2 == e) continue;
                bool isPrism2 = (etypes[e2] == ISOPARAMETRIC_LINEAR_PRISM);
                const auto& faceNodesRef2 = isPrism2 ? prismFaceNodes : hexaFaceNodes;
                size_t numFaces2 = isPrism2 ? 5 : 6;

                for (size_t f2 = 0; f2 < numFaces2; ++f2) {
                    vector<size_t> faceNodes2 = getSortedFaceNodes(deqElements[e2], faceNodesRef2[f2]);
                    if (faceNodes == faceNodes2) {
                        vector<size_t> faceNodes2Ordered = getFaceNodes(deqElements[e2], faceNodesRef2[f2]);
                        vector<size_t> faceNodesRev = faceNodesOrdered;
                        reverse(faceNodesRev.begin(), faceNodesRev.end());
                        vector<size_t> doubled = faceNodes2Ordered;
                        doubled.insert(doubled.end(), faceNodes2Ordered.begin(), faceNodes2Ordered.end());

                        bool opposite = false;
                        for (size_t s = 0; s < faceNodes2Ordered.size(); ++s) {
                            if (equal(faceNodesRev.begin(), faceNodesRev.end(), doubled.begin() + static_cast<long>(s))) {
                                opposite = true;
                                break;
                            }
                        }
                        if (opposite) {
                            deqElementNeighbors[e][f] = static_cast<int64_t>(e2);
                            found = true;
                            break;
                        }
                    }
                }
                if (found) break;
            }

            if (!found) {
                vector<double> normal = computeFaceNormal(faceNodesOrdered);
                if (abs(normal[0] + 1.0) < EPSILON && abs(normal[1]) < EPSILON && abs(normal[2]) < EPSILON) {
                    deqElementNeighbors[e][f] = LEFT_OUTSIDE;
                } else if (abs(normal[0] - 1.0) < EPSILON && abs(normal[1]) < EPSILON && abs(normal[2]) < EPSILON) {
                    deqElementNeighbors[e][f] = RIGHT_OUTSIDE;
                } else if (abs(normal[1] + 1.0) < EPSILON && abs(normal[0]) < EPSILON && abs(normal[2]) < EPSILON) {
                    deqElementNeighbors[e][f] = BOTTOM_OUTSIDE;
                } else if (abs(normal[1] - 1.0) < EPSILON && abs(normal[0]) < EPSILON && abs(normal[2]) < EPSILON) {
                    deqElementNeighbors[e][f] = TOP_OUTSIDE;
                } else if (abs(normal[2] + 1.0) < EPSILON && abs(normal[0]) < EPSILON && abs(normal[1]) < EPSILON) {
                    deqElementNeighbors[e][f] = BACK_OUTSIDE;
                } else if (abs(normal[2] - 1.0) < EPSILON && abs(normal[0]) < EPSILON && abs(normal[1]) < EPSILON) {
                    deqElementNeighbors[e][f] = FRONT_OUTSIDE;
                } else {
                    double cx = 0, cy = 0, cz = 0;
                    for (size_t node : faceNodesOrdered) {
                        cx += px[node];
                        cy += py[node];
                        cz += pz[node];
                    }
                    cx /= faceNodesOrdered.size();
                    cy /= faceNodesOrdered.size();
                    cz /= faceNodesOrdered.size();
                    if (abs(cy) < EPSILON) deqElementNeighbors[e][f] = BOTTOM_OUTSIDE;
                    else if (abs(cy - (iDim_j - 1)) < EPSILON) deqElementNeighbors[e][f] = TOP_OUTSIDE;
                    else if (abs(cx) < EPSILON) deqElementNeighbors[e][f] = LEFT_OUTSIDE;
                    else if (abs(cx - (iDim_i - 1)) < EPSILON) deqElementNeighbors[e][f] = RIGHT_OUTSIDE;
                    else if (abs(cz) < EPSILON) deqElementNeighbors[e][f] = BACK_OUTSIDE;
                    else if (abs(cz - (iDim_k - 1)) < EPSILON) deqElementNeighbors[e][f] = FRONT_OUTSIDE;
                    else deqElementNeighbors[e][f] = IRREGULAR_OUTSIDE;
                }
            }
        }
    }
    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());


    // Node boundary flags
    // ===================
    for (size_t k = 0; k < iDim_k; ++k) {
        for (size_t j = 0; j < iDim_j; ++j) {
            for (size_t i = 0; i < iDim_i; ++i) {
                const size_t iNode = k * iDim_k2 + iDim_j * j + i;
                int8_t bBoundary = NOT;

                if (bSkewed) {
                    if (abs(px[iNode]) < EPSILON) {
                        if (abs(py[iNode]) < EPSILON && abs(pz[iNode]) < EPSILON) bBoundary = CNR1;
                        else if (abs(py[iNode]) < EPSILON && abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR5;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && abs(pz[iNode]) < EPSILON) bBoundary = CNR4;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR8;
                        else if (abs(py[iNode]) < EPSILON && pz[iNode] > EPSILON && pz[iNode] < (iDim_k - 1) - EPSILON) bBoundary = EDGE5;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && pz[iNode] > EPSILON && pz[iNode] < (iDim_k - 1) - EPSILON) bBoundary = EDGE8;
                        else if (abs(pz[iNode]) < EPSILON && py[iNode] > EPSILON && py[iNode] < (iDim_j - 1) - EPSILON) bBoundary = EDGE4;
                        else if (abs(pz[iNode] - (iDim_k - 1)) < EPSILON && py[iNode] > EPSILON && py[iNode] < (iDim_j - 1) - EPSILON) bBoundary = EDGE12;
                        else bBoundary = LEFT_OUTSIDE;
                    } else if (abs(px[iNode] - (iDim_i - 1)) < EPSILON) {
                        if (abs(py[iNode]) < EPSILON && abs(pz[iNode]) < EPSILON) bBoundary = CNR2;
                        else if (abs(py[iNode]) < EPSILON && abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR6;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && abs(pz[iNode]) < EPSILON) bBoundary = CNR3;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && abs(pz[iNode] - (iDim_k - 1)) < EPSILON) bBoundary = CNR7;
                        else if (abs(py[iNode]) < EPSILON && pz[iNode] > EPSILON && pz[iNode] < (iDim_k - 1) - EPSILON) bBoundary = EDGE6;
                        else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON && pz[iNode] > EPSILON && pz[iNode] < (iDim_k - 1) - EPSILON) bBoundary = EDGE7;
                        else if (abs(pz[iNode]) < EPSILON && py[iNode] > EPSILON && py[iNode] < (iDim_j - 1) - EPSILON) bBoundary = EDGE2;
                        else if (abs(pz[iNode] - (iDim_k - 1)) < EPSILON && py[iNode] > EPSILON && py[iNode] < (iDim_j - 1) - EPSILON) bBoundary = EDGE10;
                        else bBoundary = RIGHT_OUTSIDE;
                    } else if (abs(py[iNode]) < EPSILON) {
                        if (abs(pz[iNode]) < EPSILON && px[iNode] > EPSILON && px[iNode] < (iDim_i - 1) - EPSILON) bBoundary = EDGE1;
                        else if (abs(pz[iNode] - (iDim_k - 1)) < EPSILON && px[iNode] > EPSILON && px[iNode] < (iDim_i - 1) - EPSILON) bBoundary = EDGE9;
                        else bBoundary = BOTTOM_OUTSIDE;
                    } else if (abs(py[iNode] - (iDim_j - 1)) < EPSILON) {
                        if (abs(pz[iNode]) < EPSILON && px[iNode] > EPSILON && px[iNode] < (iDim_i - 1) - EPSILON) bBoundary = EDGE3;
                        else if (abs(pz[iNode] - (iDim_k - 1)) < EPSILON && px[iNode] > EPSILON && px[iNode] < (iDim_i - 1) - EPSILON) bBoundary = EDGE11;
                        else bBoundary = TOP_OUTSIDE;
                    } else if (abs(pz[iNode]) < EPSILON) {
                        bBoundary = BACK_OUTSIDE;
                    } else if (abs(pz[iNode] - (iDim_k - 1)) < EPSILON) {
                        bBoundary = FRONT_OUTSIDE;
                    }
                } else {
                    if (i == 0 && j == 0 && k == 0) bBoundary = CNR1;
                    else if (i == iDim_i - 1 && j == 0 && k == 0) bBoundary = CNR2;
                    else if (i == iDim_i - 1 && j == iDim_j - 1 && k == 0) bBoundary = CNR3;
                    else if (i == 0 && j == iDim_j - 1 && k == 0) bBoundary = CNR4;
                    else if (i == 0 && j == 0 && k == iDim_k - 1) bBoundary = CNR5;
                    else if (i == iDim_i - 1 && j == 0 && k == iDim_k - 1) bBoundary = CNR6;
                    else if (i == iDim_i - 1 && j == iDim_j - 1 && k == iDim_k - 1) bBoundary = CNR7;
                    else if (i == 0 && j == iDim_j - 1 && k == iDim_k - 1) bBoundary = CNR8;
                    else if (j == 0 && k == 0 && i > 0 && i < iDim_i - 1) bBoundary = EDGE1;
                    else if (i == iDim_i - 1 && k == 0 && j > 0 && j < iDim_j - 1) bBoundary = EDGE2;
                    else if (j == iDim_j - 1 && k == 0 && i > 0 && i < iDim_i - 1) bBoundary = EDGE3;
                    else if (i == 0 && k == 0 && j > 0 && j < iDim_j - 1) bBoundary = EDGE4;
                    else if (i == 0 && j == 0 && k > 0 && k < iDim_k - 1) bBoundary = EDGE5;
                    else if (i == iDim_i - 1 && j == 0 && k > 0 && k < iDim_k - 1) bBoundary = EDGE6;
                    else if (i == iDim_i - 1 && j == iDim_j - 1 && k > 0 && k < iDim_k - 1) bBoundary = EDGE7;
                    else if (i == 0 && j == iDim_j - 1 && k > 0 && k < iDim_k - 1) bBoundary = EDGE8;
                    else if (j == 0 && k == iDim_k - 1 && i > 0 && i < iDim_i - 1) bBoundary = EDGE9;
                    else if (i == iDim_i - 1 && k == iDim_k - 1 && j > 0 && j < iDim_j - 1) bBoundary = EDGE10;
                    else if (j == iDim_j - 1 && k == iDim_k - 1 && i > 0 && i < iDim_i - 1) bBoundary = EDGE11;
                    else if (i == 0 && k == iDim_k - 1 && j > 0 && j < iDim_j - 1) bBoundary = EDGE12;
                    else if (k == 0 && i > 0 && i < iDim_i - 1 && j > 0 && j < iDim_j - 1) bBoundary = BACK_OUTSIDE;
                    else if (k == iDim_k - 1 && i > 0 && i < iDim_i - 1 && j > 0 && j < iDim_j - 1) bBoundary = FRONT_OUTSIDE;
                    else if (i == 0 && j > 0 && j < iDim_j - 1 && k > 0 && k < iDim_k - 1) bBoundary = LEFT_OUTSIDE;
                    else if (i == iDim_i - 1 && j > 0 && j < iDim_j - 1 && k > 0 && k < iDim_k - 1) bBoundary = RIGHT_OUTSIDE;
                    else if (j == 0 && i > 0 && i < iDim_i - 1 && k > 0 && k < iDim_k - 1) bBoundary = BOTTOM_OUTSIDE;
                    else if (j == iDim_j - 1 && i > 0 && i < iDim_i - 1 && k > 0 && k < iDim_k - 1) bBoundary = TOP_OUTSIDE;
                }
                vset.BFlag(iNode, bBoundary);
            }
        }
    }
    vset.ResizeNodes(nodes);

    // Validity check for skewed
    if (bSkewed) {
        for (size_t e = 0; e < iNrOfElements; ++e) {
            set<size_t> unique_nodes(deqElements[e].begin(), deqElements[e].end());
            if (unique_nodes.size() != deqElements[e].size()) {
                cout << "Degenerate element detected in element " << e << endl;
            }
        }
    }

    // Materials
    vector<int32_t> pmtrl(vset.Elements(), 1);
    iElement = 0;
    for (size_t k = 0; k < iDim_k - 1; ++k) {
        for (size_t j = 0; j < iDim_j - 1; ++j) {
            for (size_t i = 0; i < iDim_i - 1; ++i) {
                if (i == 1 && j == 1) {
                    pmtrl[iElement++] = 7;
                    pmtrl[iElement++] = 7;
                    continue;
                }
                pmtrl[iElement++] = 1;
            }
        }
    }
    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());

    // Element number
    PropertyData elmt_nums(ELEMENT, SCALAR, 3U);
    elmt_nums.Reserve(vset.Elements());
    for (size_t i = 0; i < vset.Elements(); ++i) {
        pushBack(elmt_nums, makeScalar(ANY, i));
    }
    vset.AddData("element number", elmt_nums);

    cout << "\ncreate_Prism_Hexa_VSet: model 'Prism_Hexa': done." << endl;
}









/**
    Prism that lies on its front side.
    No neighbors, but ony (negative) boundary flags.
*/
void create_1Prism_VSet(VSet<3U>& vset, bool bSkewed )
{
  	IsoparametricLinearPrism iso_prism;
  	
  	//elements
    uint32_t nodes(iso_prism.Nodes()); // number of nodes in single Prism mesh
  	
   // resize necessitates element type to be defined
    vset.SingleElementType(iso_prism.ElementType());
  	vset.Resize( iso_prism.Nodes(), iso_prism.Neighbors(), static_cast<int8_t>(iso_prism.ElementType()), nodes, 1 );

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
  	 for ( size_t i = 0; i < 6; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.ResizeBFlags();
  	
  	//--------------------------ELEMENTS
    //define prism element, assign nodes per element (0..n-1)
    deque<vector<size_t> > deqElements = { {0, 1, 2, 3, 4, 5} };
  	vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    deque<vector<int64_t> > deqElementNeighbors = { {BACK_OUTSIDE, BOTTOM_OUTSIDE, RIGHT_OUTSIDE, IRREGULAR, FRONT_OUTSIDE } };
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.BFlag( 0, CNR1 );
  	vset.BFlag( 1, CNR2 );
  	vset.BFlag( 2, CNR3 );
  	vset.BFlag( 3, CNR5 );
  	vset.BFlag( 4, CNR6 );
  	vset.BFlag( 5, CNR7 );
    
    //-------------------------MATERIALS
    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

    cout <<"\n"<<"create_1Prism_VSet: model 'un-named': done."<< endl;
//    vset.Out();
    
} // end create_1Prism_VSet






void create_Prism_VSet(VSet<3U>& vset, bool bSkewed)
{
    // --- Dimensions of the prism stack ---
    const size_t nx = 3; // nodes in x
    const size_t ny = 3; // nodes in y
    const size_t nz = 3; // nodes in z

    size_t nNodes = nx * ny * nz;
    size_t nCells = (nx-1)*(ny-1)*(nz-1) * 2; // two prisms per brick

    // --- Coordinates ---
    std::deque<double> px(nNodes), py(nNodes), pz(nNodes);
    for (size_t k = 0; k < nz; ++k)
    for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
        size_t idx = k*(nx*ny) + j*nx + i;
        px[idx] = i + (bSkewed ? ((rand()%1000)/1000.0)*PERTURBATION : 0.0);
        py[idx] = j + (bSkewed ? ((rand()%1000)/1000.0)*PERTURBATION : 0.0);
        pz[idx] = k + (bSkewed ? ((rand()%1000)/1000.0)*PERTURBATION : 0.0);
    }
    vset.AddXYZ(px, py, pz);

    // --- Node flags (CSMP BOX_BOUNDARY) ---
    int8_t nodeFlag = NOT;
    for (size_t k = 0; k < nz; ++k)
    for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
        size_t idx = k*(nx*ny) + j*nx + i;
        bool onX0=i==0, onX1=i==nx-1, onY0=j==0, onY1=j==ny-1, onZ0=k==0, onZ1=k==nz-1;

        if (onX0 && onY0 && onZ0) nodeFlag = CNR1;
        else if (onX1 && onY0 && onZ0) nodeFlag = CNR2;
        else if (onX1 && onY1 && onZ0) nodeFlag = CNR3;
        else if (onX0 && onY1 && onZ0) nodeFlag = CNR4;
        else if (onX0 && onY0 && onZ1) nodeFlag = CNR5;
        else if (onX1 && onY0 && onZ1) nodeFlag = CNR6;
        else if (onX1 && onY1 && onZ1) nodeFlag = CNR7;
        else if (onX0 && onY1 && onZ1) nodeFlag = CNR8;
        // Edges
        else if (onZ0 && onY0) nodeFlag = EDGE1;
        else if (onZ0 && onX1) nodeFlag = EDGE2;
        else if (onZ0 && onY1) nodeFlag = EDGE3;
        else if (onZ0 && onX0) nodeFlag = EDGE4;
        else if (onY0 && onX0) nodeFlag = EDGE5;
        else if (onY0 && onX1) nodeFlag = EDGE6;
        else if (onY1 && onX1) nodeFlag = EDGE7;
        else if (onY1 && onX0) nodeFlag = EDGE8;
        else if (onZ1 && onY0) nodeFlag = EDGE9;
        else if (onZ1 && onX1) nodeFlag = EDGE10;
        else if (onZ1 && onY1) nodeFlag = EDGE11;
        else if (onZ1 && onX0) nodeFlag = EDGE12;
        // Faces
        else if (onZ0) nodeFlag = BACK;
        else if (onZ1) nodeFlag = FRONT;
        else if (onX0) nodeFlag = LEFT;
        else if (onX1) nodeFlag = RIGHT;
        else if (onY0) nodeFlag = BOTTOM;
        else if (onY1) nodeFlag = TOP;
        else nodeFlag = NOT;
        vset.BFlag( idx, nodeFlag );
    }

    // --- Elements (two prisms per brick) ---
    std::deque<std::vector<size_t>> deqElements;
    std::vector<int8_t> vecElementTypes;
    for (size_t k=0;k<nz-1;++k)
    for (size_t j=0;j<ny-1;++j)
    for (size_t i=0;i<nx-1;++i)
    {
        size_t n000 = k*(nx*ny) + j*nx + i;
        size_t n100 = n000 + 1;
        size_t n010 = n000 + nx;
        size_t n110 = n010 + 1;
        size_t n001 = n000 + nx*ny;
        size_t n101 = n001 + 1;
        size_t n011 = n001 + nx;
        size_t n111 = n011 + 1;

        deqElements.push_back({n000,n100,n010,n001,n101,n011}); // lower prism
        deqElements.push_back({n100,n110,n010,n101,n111,n011}); // upper prism
        vecElementTypes.push_back(ISOPARAMETRIC_LINEAR_PRISM);
        vecElementTypes.push_back(ISOPARAMETRIC_LINEAR_PRISM);
    }
    vset.AddPlist(deqElements.begin(), deqElements.end());
    vset.AddElementTypes(vecElementTypes.begin(), vecElementTypes.end());

    // --- Neighbor connectivity --- TODO: nbors of elmt 3, 4, and 5 are incorrect
    std::deque<std::vector<int64_t>> deqNeighbors(deqElements.size(), std::vector<int64_t>(5,-1));
    size_t nElX = nx-1, nElY = ny-1, nElZ = nz-1;

    for (size_t k=0;k<nElZ;++k)
    for (size_t j=0;j<nElY;++j)
    for (size_t i=0;i<nElX;++i)
    {
        size_t base = (k*nElY*nElX + j*nElX + i) * 2;
        size_t lower = base;
        size_t upper = base + 1;

        // --- Lower prism neighbors ---
        // face 0: back
        deqNeighbors[lower][0] = (k==0) ? BACK_OUTSIDE : static_cast<int64_t>(lower - 2*nElX*nElY);
        // face 1: bottom
        deqNeighbors[lower][1] = (j==0) ? BOTTOM_OUTSIDE : static_cast<int64_t>(lower - 2*nElX);
        // face 2: right
        deqNeighbors[lower][2] = (i==nElX-1) ? RIGHT_OUTSIDE : static_cast<int64_t>(lower + 2);
        // face 3: left
        deqNeighbors[lower][3] = (i==0) ? LEFT_OUTSIDE : static_cast<int64_t>(lower - 2);
        // face 4: front (top neighbor)
        deqNeighbors[lower][4] = static_cast<int64_t>(upper);

        // --- Upper prism neighbors ---
        // face 0: back
        deqNeighbors[upper][0] = (k==0) ? BACK_OUTSIDE : static_cast<int64_t>(upper - 2*nElX*nElY);
        // face 1: top
        deqNeighbors[upper][1] = (j==nElY-1) ? TOP_OUTSIDE : static_cast<int64_t>(upper + 2*nElX);
        // face 2: left
        deqNeighbors[upper][2] = (i==0) ? LEFT_OUTSIDE : static_cast<int64_t>(upper - 2);
        // face 3: right
        deqNeighbors[upper][3] = (i==nElX-1) ? RIGHT_OUTSIDE : static_cast<int64_t>(upper + 2);
        // face 4: bottom neighbor
        deqNeighbors[upper][4] = static_cast<int64_t>(lower);
    }
    vset.AddPfverts(deqNeighbors.begin(), deqNeighbors.end());
}



/**
       Decomposition of a hexahedron into 6 tetrahedra.
       Only 6 elements!
       Coordinate range from -1 - 1
       
       Model is suitable for testing neighbor connectivity and boundary node flagging.
       
       @test verified SKM 5/9/25
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
    // CSMP_FEM_conventions.pdf
    //               CNR  1     2     3     4     5     6     7     8
   	deque<double> px = {-1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0};
  	deque<double> py = {-1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0};
  	deque<double> pz = {-1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0};

   	//load nodes
  	vset.AddXYZ( px, py, pz );
 
    vset.ResizeBFlags();
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


    //--------------------------ELEMENT NEIGHBORS
    deque<vector<int64_t>> deqElementNeighbors{ {1, LEFT_OUTSIDE, BOTTOM_OUTSIDE, BACK_OUTSIDE}, // Element 0
                                                {4, 0, 2, BOTTOM_OUTSIDE},                       // Element 1
                                                {5, LEFT_OUTSIDE, FRONT_OUTSIDE, 1},             // Element 2
                                                {TOP_OUTSIDE, 4, RIGHT_OUTSIDE, BACK_OUTSIDE},   // Element 3
                                                {RIGHT_OUTSIDE, 5, 1, 3},                        // Element 4
                                                {TOP_OUTSIDE, 2, FRONT_OUTSIDE, 4}               // Element 5
                                              };
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());

    vector<int32_t> pmtrl( vset.Elements(), 1 );
    vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
    
    cout <<"\n"<<"testCreateTetra_VSet: model 'Tetra': done."<< endl;
//    vset.Out();

 } // end testCreateTetra_VSet





void create_Pyramid_VSet(VSet<3U>& vset, bool bSkewed)
{
    IsoparametricLinearPyramid iso_pyramid;
    
    // Node dimensions of box
    const int iDim_k(4); // k - height
    const int iDim_j(4); // j - width
    const int iDim_i(4); // i - length
    const int iDim_k2(iDim_k * iDim_k); // k^2
    const int iDim_km1_2((iDim_k - 1) * (iDim_k - 1)); // (k-1)^2
    
    const int iNrOfCells((iDim_i - 1) * (iDim_j - 1) * (iDim_k - 1));
    const int iNrOfElements(iNrOfCells * 6U); // 6 pyramids per cell
    
    // Total nodes: 64 (4x4x4 grid) + 27 barycenters (one per cell)
    size_t nodes(iDim_i * iDim_j * iDim_k + iNrOfCells);
    
    // Initialize VSet
    vset.SingleElementType(iso_pyramid.ElementType());
    vset.Resize(iso_pyramid.Nodes(), // 5 nodes per pyramid
                iso_pyramid.Neighbors(), // 5 neighbors per pyramid
                iso_pyramid.ElementType(),
                nodes, iNrOfElements);
    
    //-----------------------NODES
    // Define node coordinates
    deque<double> px(nodes);
    deque<double> py(nodes);
    deque<double> pz(nodes);
    
    // Initialize random number generator for perturbations
    mt19937 rng(12345); // Fixed seed for reproducibility
    uniform_real_distribution<double> dist(0.0, 1.0);
    
    // Grid nodes
    for (size_t k{0}; k < iDim_k; k++) // z
        for (size_t j{0}; j < iDim_j; j++) // y
            for (size_t i{0}; i < iDim_i; i++) // x
            {
                size_t iNode = k * iDim_k2 + j * iDim_j + i;
                if (bSkewed)
                {
                    px[iNode] = i + dist(rng) * PERTURBATION;
                    py[iNode] = j + dist(rng) * PERTURBATION;
                    pz[iNode] = k + dist(rng) * PERTURBATION;
                }
                else
                {
                    px[iNode] = i;
                    py[iNode] = j;
                    pz[iNode] = k;
                }
            }
    
    // Barycenter nodes
    for (size_t k{0ul}; k < iDim_k - 1; k++) // z
        for (size_t j{0ul}; j < iDim_j - 1; j++) // y
            for (size_t i{0ul}; i < iDim_i - 1; i++) // x
            {
                size_t iNode = k * iDim_km1_2 + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k;
                if (bSkewed)
                {
                    px[iNode] = i + 0.5 + dist(rng) * PERTURBATION;
                    py[iNode] = j + 0.5 + dist(rng) * PERTURBATION;
                    pz[iNode] = k + 0.5 + dist(rng) * PERTURBATION;
                }
                else
                {
                    px[iNode] = i + 0.5;
                    py[iNode] = j + 0.5;
                    pz[iNode] = k + 0.5;
                }
            }
    
    // Load nodes
    vset.AddXYZ(px, py, pz);
    vset.ResizeBFlags();
    
    //--------------------------ELEMENTS
    const bool verbose{false};
    if (verbose)
        cout << "\n" << "create_Pyramid_VSet: printing coordinates of selected pyramid elements.";
    
    // Define element types
    vector<int8_t> vecElementTypes(1, ISOPARAMETRIC_LINEAR_PYRAMID);
    vset.AddElementTypes(vecElementTypes.begin(), vecElementTypes.end());
    
    // Define pyramid elements (0 to 161), assign nodes per element
    deque<vector<size_t>> deqElements(iNrOfElements);
    for (size_t k{0}; k < iDim_k - 1; k++) // z
        for (size_t j{0}; j < iDim_j - 1; j++) // y
            for (size_t i{0}; i < iDim_i - 1; i++) // x
            {
                // Pyramid 0
                size_t iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * k + iDim_j * j + i; // (i,j,k)
                deqElements[iElement][1] = iDim_k2 * k + iDim_j * j + i + 1; // (i+1,j,k)
                deqElements[iElement][2] = iDim_k2 * k + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k)
                deqElements[iElement][3] = iDim_k2 * k + iDim_j * (j + 1) + i; // (i,j+1,k)
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
                
                // Pyramid 1
                iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * (k + 1) + iDim_j * j + i; // (i,j,k+1)
                deqElements[iElement][1] = iDim_k2 * (k + 1) + iDim_j * j + i + 1; // (i+1,j,k+1)
                deqElements[iElement][2] = iDim_k2 * k + iDim_j * j + i + 1; // (i+1,j,k)
                deqElements[iElement][3] = iDim_k2 * k + iDim_j * j + i; // (i,j,k)
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
                
                // Pyramid 2
                if (verbose) cout << "Element:(i=" << i << ",j=" << j << ",k=" << k << ")" << iDim_km1_2 * k + (iDim_j - 1) * j + i + 27 * 2 << endl;
                iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * k + iDim_j * j + i + 1; // (i+1,j,k)
                if (verbose) cout << "\tCoord [" << iElement << "][0]: " << deqElements[iElement][0] << endl;
                deqElements[iElement][1] = iDim_k2 * (k + 1) + iDim_j * j + i + 1; // (i+1,j,k+1)
                if (verbose) cout << "\tCoord [" << iElement << "][1]: " << deqElements[iElement][1] << endl;
                deqElements[iElement][2] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k+1)
                if (verbose) cout << "\tCoord [" << iElement << "][2]: " << deqElements[iElement][2] << endl;
                deqElements[iElement][3] = iDim_k2 * k + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k)
                if (verbose) cout << "\tCoord [" << iElement << "][3]: " << deqElements[iElement][3] << endl;
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
                if (verbose) cout << "\tCoord [" << iElement << "][4]: " << deqElements[iElement][4] << endl;
                
                // Pyramid 3
                iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * k + iDim_j * (j + 1) + i; // (i,j+1,k)
                deqElements[iElement][1] = iDim_k2 * k + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k)
                deqElements[iElement][2] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k+1)
                deqElements[iElement][3] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i; // (i,j+1,k+1)
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
                
                // Pyramid 4
                iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * (k + 1) + iDim_j * j + i; // (i,j,k+1)
                deqElements[iElement][1] = iDim_k2 * k + iDim_j * j + i; // (i,j,k)
                deqElements[iElement][2] = iDim_k2 * k + iDim_j * (j + 1) + i; // (i,j+1,k)
                deqElements[iElement][3] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i; // (i,j+1,k+1)
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
                
                // Pyramid 5 (reoriented for correct face 1 and face 3)
                iElement = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5;
                deqElements[iElement].resize(5);
                deqElements[iElement][0] = iDim_k2 * (k + 1) + iDim_j * j + i; // (i,j,k+1)
                deqElements[iElement][1] = iDim_k2 * (k + 1) + iDim_j * j + i + 1; // (i+1,j,k+1)
                deqElements[iElement][2] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i + 1; // (i+1,j+1,k+1)
                deqElements[iElement][3] = iDim_k2 * (k + 1) + iDim_j * (j + 1) + i; // (i,j+1,k+1)
                deqElements[iElement][4] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iDim_i * iDim_j * iDim_k; // barycenter
            }
    
    vset.AddPlist(deqElements.begin(), deqElements.end());
    
    //---------------------------------NEIGHBORS
    // Define neighbors
    deque<vector<int64_t>> deqElementNeighbors(iNrOfElements);
    for (int64_t k{0}; k < iDim_k - 1; k++) // z
        for (int64_t j{0}; j < iDim_j - 1; j++) // y
            for (int64_t i{0}; i < iDim_i - 1; i++) // x
            {
                // Pyramid 0
                size_t iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4; // face 3
                deqElementNeighbors[iElement][4] = (k == 0 ? BACK_OUTSIDE : iDim_km1_2 * (k - 1) + (iDim_j - 1) * j + i + iNrOfCells * 5); // face 4 (base)
                
                // Pyramid 1
                iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 0; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4; // face 3
                deqElementNeighbors[iElement][4] = (j == 0 ? BOTTOM_OUTSIDE : iDim_km1_2 * k + (iDim_j - 1) * (j - 1) + i + iNrOfCells * 3); // face 4 (base)
                
                // Pyramid 2
                iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 0; // face 3
                deqElementNeighbors[iElement][4] = (i == iDim_i - 2 ? RIGHT_OUTSIDE : iDim_km1_2 * k + (iDim_j - 1) * j + i + 1 + iNrOfCells * 4); // face 4 (base)
                
                // Pyramid 3
                iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 0; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4; // face 3
                deqElementNeighbors[iElement][4] = (j == iDim_j - 2 ? TOP_OUTSIDE : iDim_km1_2 * k + (iDim_j - 1) * (j + 1) + i + iNrOfCells * 1); // face 4 (base)
                
                // Pyramid 4
                iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 0; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5; // face 3
                deqElementNeighbors[iElement][4] = (i == 0 ? LEFT_OUTSIDE : iDim_km1_2 * k + (iDim_j - 1) * j + i - 1 + iNrOfCells * 2); // face 4 (base)
                
                // Pyramid 5
                iElement = static_cast<size_t>(iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 5);
                deqElementNeighbors[iElement].resize(5);
                deqElementNeighbors[iElement][0] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 1; // face 0
                deqElementNeighbors[iElement][1] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 2; // face 1
                deqElementNeighbors[iElement][2] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 3; // face 2
                deqElementNeighbors[iElement][3] = iDim_km1_2 * k + (iDim_j - 1) * j + i + iNrOfCells * 4; // face 3
                deqElementNeighbors[iElement][4] = (k == iDim_k - 2 ? FRONT_OUTSIDE : iDim_km1_2 * (k + 1) + (iDim_j - 1) * j + i + iNrOfCells * 0); // face 4 (base)
            }
    
    vset.AddPfverts(deqElementNeighbors.begin(), deqElementNeighbors.end());
    
    //-----------------------------------------------------NODE BOUNDARIES
    // Define boundary flags for nodes
    for (int k{0}; k < iDim_k; k++) // z
        for (int j{0}; j < iDim_j; j++) // y
            for (int i{0}; i < iDim_i; i++) // x
            {
                int8_t bBoundary = NOT;
                
                if (k == 0)
                {
                    if (j == 0)
                    {
                        if (i == 0) bBoundary = CNR1;
                        else if (i == (iDim_i - 1)) bBoundary = CNR2;
                        else bBoundary = EDGE1;
                    }
                    else if (j == (iDim_j - 1))
                    {
                        if (i == 0) bBoundary = CNR4;
                        else if (i == (iDim_i - 1)) bBoundary = CNR3;
                        else bBoundary = EDGE3;
                    }
                    else
                    {
                        if (i == 0) bBoundary = EDGE4;
                        else if (i == (iDim_i - 1)) bBoundary = EDGE2;
                        else bBoundary = BACK_OUTSIDE;
                    }
                }
                else if (k == (iDim_k - 1))
                {
                    if (j == 0)
                    {
                        if (i == 0) bBoundary = CNR5;
                        else if (i == (iDim_i - 1)) bBoundary = CNR6;
                        else bBoundary = EDGE9;
                    }
                    else if (j == (iDim_j - 1))
                    {
                        if (i == 0) bBoundary = CNR8;
                        else if (i == (iDim_i - 1)) bBoundary = CNR7;
                        else bBoundary = EDGE11;
                    }
                    else
                    {
                        if (i == 0) bBoundary = EDGE12;
                        else if (i == (iDim_i - 1)) bBoundary = EDGE10;
                        else bBoundary = FRONT_OUTSIDE;
                    }
                }
                else
                {
                    if (j == 0)
                    {
                        if (i == 0) bBoundary = EDGE5;
                        else if (i == (iDim_i - 1)) bBoundary = EDGE6;
                        else bBoundary = BOTTOM_OUTSIDE;
                    }
                    else if (j == (iDim_j - 1))
                    {
                        if (i == 0) bBoundary = EDGE8;
                        else if (i == (iDim_i - 1)) bBoundary = EDGE7;
                        else bBoundary = TOP_OUTSIDE;
                    }
                    else
                    {
                        if (i == 0) bBoundary = LEFT_OUTSIDE;
                        else if (i == (iDim_i - 1)) bBoundary = RIGHT_OUTSIDE;
                        else ; // do nothing: no boundary
                    }
                }
                
                if (bBoundary != NOT)
                {
                    auto iNode = iDim_k2 * k + iDim_j * j + i;
                    vset.BFlag( static_cast<size_t>(iNode), bBoundary);
                }
            }
    
    // Explicitly set boundary flags for barycenter nodes
    for (size_t iNode = iDim_i * iDim_j * iDim_k; iNode < nodes; ++iNode)
    {
        vset.BFlag(iNode, NOT);
    }
    
    //-------------------------MATERIALS
    vector<int32_t> pmtrl(vset.Elements(), 1);
    vset.AddPmtrl(pmtrl.begin(), pmtrl.end());
    
    vset.EstablishZeroBasedNumbering();
    
    cout << "\n" << "create_Pyramid_VSet: model 'Pyra': done." << endl;
    
} // end create_Pyramid_VSet





// helper Function to generate the corner points for a 3D corner-point grid
static  vector<array<Point<3U>,8>> generateCornerPointGrid( double dx, double dy, double dz )
  {
      // Define the dimensions of the grid
      constexpr int I = 6; // Number of cells in the I-direction
      constexpr int J = 8; // Number of cells in the J-direction
      constexpr int K = 4; // Number of cells in the K-direction
      
      vector<array<Point<3U>, 8>> grid;

      // Loop over each cell in the grid
      for (int k = 0; k < K; ++k) {
          for (int j = 0; j < J; ++j) {
              for (int i = 0; i < I; ++i) {
                  // Calculate the coordinates of the 8 corner points of the current cell
                  array<Point<3U>,8> corners = {
                      Point<3U>{ i * dx,     j * dy,     k * dz     }, // (0, 0, 0)
                      Point<3U>{ (i+1) * dx, j * dy,     k * dz     }, // (1, 0, 0)
                      Point<3U>{ (i+1) * dx, (j+1) * dy, k * dz     }, // (1, 1, 0)
                      Point<3U>{ i * dx,     (j+1) * dy, k * dz     }, // (0, 1, 0)
                      Point<3U>{ i * dx,     j * dy,     (k+1) * dz }, // (0, 0, 1)
                      Point<3U>{ (i+1) * dx, j * dy,     (k+1) * dz }, // (1, 0, 1)
                      Point<3U>{ (i+1) * dx, (j+1) * dy, (k+1) * dz }, // (1, 1, 1)
                      Point<3U>{ i * dx,     (j+1) * dy, (k+1) * dz }  // (0, 1, 1)
                  };

                  // Add the corners to the grid
                  grid.push_back(corners);
              }
          }
      }

    return grid;
      
 } // end generateCornerPointGrid



/**
    Basic hexahedral grid for mesh-modifcation tests
    
    TODO: the VSet does not get initialised yet

*/
static void create_CornerPointGrid_6i_8j_4k( VSet<3U>&  )
 {
     // Define the cell dimensions
     double dx = 1.0;
     double dy = 1.0;
     double dz = 1.0;

     // Generate the corner-point grid
     vector<array<Point<3U>,8>> grid = generateCornerPointGrid(dx, dy, dz);

      // Print the corner points for each cell
      for (size_t cellIndex = 0; cellIndex < grid.size(); ++cellIndex) {
          cout << "Cell " << cellIndex << " corner points:\n";
          for (const auto& corner : grid[cellIndex]) {
              cout << "(" << corner[0] << ", " << corner[1] << ", " << corner[2] << ")\n";
          }
          cout << "\n";
      }
    cout <<"\n"<<"create_CornerPointGrid_6i_8j_4k: model 'un-named': done."<< endl;

 } // end create_CornerPointGrid




/**
      From ANSYS meshed FracBox with boundary surfaces and several tetrahedra
      that span model corners so that their rectification can be explored.
      
      @attention model does not contain Neighbor information so it is created via VData functionality.
      
      @author SKM
      @date 2/2/22
*/
ModelTopology create_FracBox( VSet<3U>& vset )
 {
    set<string> femTypes_matrix{"ISOPARAMETRIC_LINEAR_TETRAHEDRON"},
                femTypes_surfaces{"ISOPARAMETRIC_LINEAR_TRIANGLE"},
                femTypes_lines{"ISOPARAMETRIC_LINEAR_BAR"};

    // MODEL TOPOLOGY
    const bool isoparametric = true;
    ModelTopology topology( isoparametric );
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

  // boundary flags
  const int8_t bf(IRREGULAR);
  // side boundaries to start with
  vector<int8_t>  bflags{bf,bf,0,0,0,bf,0,bf,0,bf,bf,0,0,0,0,0,bf,bf,0,bf,0,bf,0,0,0,0,0,0,0,0,0,0,0,0,0,0,bf,0,bf,0,bf,0,bf,bf,bf,0,bf,0,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,bf,bf,bf,bf,0,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,0,bf,0,0,bf,bf,bf,bf,bf,bf,bf,0,0,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,0,bf,0,0,0,0,bf,0,0,bf,bf,0,0,0,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,0,bf,bf,bf,0,0,0,0,bf,bf,0,bf,bf,0,bf,bf,bf,bf,0,0,bf,0,bf,bf,bf,bf,0,0,0,0,bf,bf,0,0,bf,bf,0,0,0,0,0,0,0,bf,bf,bf,0,0,0,bf,bf,bf,bf,bf,0,bf,bf,bf,0,0,bf,0,0,0,0,0,0,0,0,0,bf,bf,bf,bf,0,0,0,bf,0,0,0,bf,bf,bf,0,0,0,bf,0,bf,bf,bf,bf,0,bf,bf,bf,0,bf,bf,bf,bf,0,0,0,0,bf,0,0,0,0,0,0,0,bf,bf,bf,bf,0,0,0,bf,0,0,0,bf,bf,bf,0,0,0,bf,0,bf,bf,bf,0,bf,bf,0,0,0,0,0,bf,0,0,0,bf,bf,0,bf,0,bf,bf,0,bf,0,bf,bf,0,bf,bf,0,bf,bf,0,bf,0,0,bf,0,bf,bf,bf,bf,bf,bf,bf,0,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf,bf};


  // BOX_BOUNDARY flags
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

  // PFVERTS - neighbor connectivity is created automatically
  vset.EstablishElementConnectivity3D();

  // TOPOTYPE geometry flags
  vset.InitialiseNodeTopologyIdentifiers();
  

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
  
  return topology;

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
