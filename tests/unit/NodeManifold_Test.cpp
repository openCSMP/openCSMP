#include "NodeManifold_Test.h"
#include "Node.h"
#include "NodeManifold.h"
#include "plf_colony.h"

using namespace std;

namespace csmp{

void NodeManifold_Test::Test_Classify2D()
 {
    LocalVariables nvars;

    // CONSTRUCTOR arguments: size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT, TOPOTYPE = MESH_VERTEX );
    Node<2> n1( 0, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_LINE );
    Node<2> n2( 1, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_LINE );
    Node<2> n3( 2, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_LINE );
    Node<2> n4( 3, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_LINE );

    NodeManifold<2> nmf( n1, n2 );
    n1.Assign( nmf );
    n2.Assign( nmf );

    // 0. Standard Split Boundary (Clean)
    n1.Attribute(PERIMETER_LINE);
    n2.Attribute(PERIMETER_LINE);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY );
    
    // 1. Stand Alone ManifoldNode
    n1.Attribute(PERIMETER_POINT);
    n2.Attribute(PERIMETER_POINT);
    _test( nmf.Classify() == ManifoldType::STAND_ALONE );
    // same
    n1.Attribute(INTERIOR_POINT);
    n2.Attribute(INTERIOR_POINT);
    _test( nmf.Classify() == ManifoldType::STAND_ALONE );

    // 2. Split Boundary with Internal Line Mesh (Equidimensional Inclusion)
    // split boundary tip inside of the model domain (point)
    n1.Attribute(INTERIOR_LINE);
    n2.Attribute(PERIMETER_POINT);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );

    // 3. Split Boundary with Internal Mesh (Lower-Dim Inclusion)
    nmf.Add(&n3);
    n1.Attribute(PERIMETER_LINE);
    n2.Attribute(PERIMETER_LINE);
    n3.Attribute(INTERIOR_LINE);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );

    // 4. Split Boundary Termination (Model Boundary)
    // against tip of internal boundary
    n1.Attribute(PERIMETER_LINE);
    n2.Attribute(PERIMETER_LINE);
    n3.Attribute(INTERIOR_POINT);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );
    nmf.Remove(&n3);
    n1.Attribute(EXTERIOR_POINT);
    n2.Attribute(EXTERIOR_POINT);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );
    // with internal line
    nmf.Add(&n3);
    n1.Attribute(EXTERIOR_POINT);
    n2.Attribute(EXTERIOR_POINT);
    n3.Attribute(EXTERIOR_POINT); 
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );

    // 5. Split Boundary Crossing (X-Junction); logic: any crossing must be a point
    nmf.Add(&n4);
    n1.Attribute(INTERIOR_POINT); n2.Attribute(INTERIOR_POINT);
    n3.Attribute(INTERIOR_POINT); n4.Attribute(INTERIOR_POINT);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_CROSSING );

    // 6. Split Boundary Termination against internal boundary (points against line)
    nmf.Remove(&n4);
    n1.Attribute(INTERIOR_LINE);
    n2.Attribute(INTERIOR_POINT);
    n3.Attribute(INTERIOR_POINT);
    _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_CROSSING ); // FAIL

    // 7. Multi-Split Boundary Crossing (Star Junction)
    // 3 Split boundaries crossing at a single point (6 mediated point branches)
    Node<2> n5( 4, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_POINT );
    Node<2> n6( 5, Point<2>(3.,3.), nvars, INTERNAL, INTERIOR_POINT );
    nmf.Add(&n4); nmf.Add(&n5); nmf.Add(&n6);
    
    n1.Attribute(INTERIOR_POINT); n2.Attribute(INTERIOR_POINT);
    n3.Attribute(INTERIOR_POINT); n4.Attribute(INTERIOR_POINT);
    // n5 and n6 already initialized as INTERIOR_POINT
    
    _test( nmf.Classify() == ManifoldType::MULTI_SB_CROSSING );

 } // end TestClassify2D





void NodeManifold_Test::Test_Classify3D()
 {
    LocalVariables nvars;

    // Standard nodes for 3D tests
    Node<3> n1(0, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);
    Node<3> n2(1, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);
    Node<3> n3(2, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);
    Node<3> n4(3, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);
    Node<3> n5(4, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);
    Node<3> n6(5, Point<3>(3,3,3), nvars, INTERNAL, INTERIOR_SURFACE);

    NodeManifold<3> nmf(n1, n2);
    n1.Assign(nmf); n2.Assign(nmf);

    // 0. Standard Clean SplitBoundary
    // Exactly two branches of perimeter surfaces.
    n1.Attribute(PERIMETER_SURFACE);
    n2.Attribute(PERIMETER_SURFACE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY);

    // 1. True STAND_ALONE Case
    // Topologically collocated but no split logic (e.g., standard internal boundary)
    n1.Attribute(PERIMETER_POINT);
    n2.Attribute(PERIMETER_POINT);
    _test(nmf.Classify() == ManifoldType::STAND_ALONE);

    // 2. T-Intersection: SplitBoundary hitting a Normal Boundary
    // 2 Perimeter Lines (Split) + 1 Interior Lines (Boundary)
    nmf.Add(&n3);
    n1.Attribute(PERIMETER_LINE);
    n2.Attribute(PERIMETER_LINE);
    n3.Attribute(INTERIOR_LINE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END);

    // 3. T-Intersection: SplitBoundary hitting a SplitBoundary (End)
    // Here, the split boundary terminates exactly at the line/point of another SB.
    // 2 Perimeter Surfaces (The Split) + 1 Perimeter Line (The edge/tip of another)
    n1.Attribute(PERIMETER_SURFACE);
    n2.Attribute(PERIMETER_SURFACE);
    n3.Attribute(PERIMETER_LINE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END);
    n3.Attribute(INTERIOR_LINE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END);

    // 4. SplitBoundary containing INTERIOR_SURFACE mesh (Equidim Inclusion)
    // 2 Perimeter Surfaces + 1 Interior Surface (NOT a T-junction, but an internal region)
    // Logic: In CSMP, if n3 was part of an internal mesh within the split,
    // we distinguish it from Case 2 by mesh connectivity or specific flagging.
    // Assuming here that 2 P_SURF + 1 I_SURF is an inclusion if specified:
    n1.Attribute(PERIMETER_SURFACE);
    n2.Attribute(PERIMETER_SURFACE);
    n3.Attribute(INTERIOR_SURFACE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH);

    // 5. Crossing SplitBoundaries (X-Junction in 3D)
    // 4 Perimeter Surface branches crossing at a line.
    nmf.Add(&n4);
    n1.Attribute(PERIMETER_SURFACE); n2.Attribute(PERIMETER_SURFACE);
    n3.Attribute(PERIMETER_SURFACE); n4.Attribute(PERIMETER_SURFACE);
    _test(nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_CROSSING);

    // 6. Multi-SB Crossing (Star Junction)
    // 6 branches of surface-matched nodes meeting.
    nmf.Add(&n5); nmf.Add(&n6);
    n1.Attribute(PERIMETER_SURFACE); n2.Attribute(PERIMETER_SURFACE);
    n3.Attribute(PERIMETER_SURFACE); n4.Attribute(PERIMETER_SURFACE);
    n5.Attribute(PERIMETER_SURFACE); n6.Attribute(PERIMETER_SURFACE);
    _test(nmf.Classify() == ManifoldType::MULTI_SB_CROSSING);

    // 7. Split Termination against Model Exterior
    NodeManifold<3> nmf2(n1, n2);
    n1.Assign(nmf2); n2.Assign(nmf2);
    n1.Attribute(EXTERIOR_SURFACE);
    n2.Attribute(EXTERIOR_SURFACE);
    _test(nmf2.Classify() == ManifoldType::SPLIT_BOUNDARY_END);

} // end TestClassify3D






void NodeManifold_Test::run()
{
  Test_Classify2D();
  Test_Classify3D();
  
  LocalVariables nvars;

  // CONSTRUCTOR arguments: size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT, TOPOTYPE = MESH_VERTEX );
  Node<3> n1( 0, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n2( 1, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n3( 2, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n4( 3, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );

  // NODE NEIGHBORS (mutual 1-2-3 connection ignored)
  set<Node<3>*> nbor_nodes_n1{&n2,&n3}, nbor_nodes_n2{&n1,&n3}, nbor_nodes_n3{&n1,&n2};
  n1.Assign( nbor_nodes_n1 );
  n2.Assign( nbor_nodes_n2 );
  n3.Assign( nbor_nodes_n3 );
  
  
  // CONSTRUCTING FIRST MANIFOLD
  NodeManifold<3> nmf( n1, n2 );
  //              ^^^^^^^^^^^^^
  n1.Assign( nmf );
  // ^^^^^^
  n2.Assign( nmf );
  _test( nmf.Branches() == 2 );
  //         ^^^^^^^^
  
  // adding the other nodes
  nmf.Add( &n3 );
  //  ^^^^^^^^^^
  n3.Assign( nmf );
  
  nmf.Add( &n4 );
  _test( nmf.Branches() == 4 );
  n4.Assign( nmf );
  
  if ( verbose_ ) nmf.Out();
  for ( uint32_t i{0}; i<nmf.Branches(); ++i )
    _test( nmf.N(i) != nullptr );
    
  pair<vector<vector<Node<3>*>>,bool> incorrect_node_connections = nmf.InterConnectedMemberNodes();
  //                                                                   ^^^^^^^^^^^^^^^^^^^^^^^^^
  _test( incorrect_node_connections.second == true );
  // there must be two clusters with 2 interconnected nodes each
  _test( incorrect_node_connections.first[0].size() == 2 );
  _test( incorrect_node_connections.first[1].size() == 2 );
 
   nmf.Remove( &n3 );
   //  ^^^^^^
  _test( nmf.Branches() == 3 );
  _test( n3.IsManifold() == false );

  // reassign test flags to 3-noded manifold
  n1.Attribute(PERIMETER_SURFACE);
  n2.Attribute(PERIMETER_SURFACE);
  n4.Attribute(INTERIOR_SURFACE);
  if ( verbose_ ) nmf.Out();
  //                  ^^^
  // testing 'ManifoldType' classification as a splitboundary with an internal mesh
  // DIAGNOSTICS                         three-surface nodes that are not interconnected
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );
  //          ^^^^^^^^^^
  _test( consistencyCheck<3>( nmf, verbose_ ) == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );

  // remaining nodes 0,1,3
  _test( nmf.N(0)->Idx() == 0 );
  _test( nmf.N(1)->Idx() == 1 );
  _test( nmf.N(2)->Idx() == 3 );
  
  // other configurations for classification
  n1.Attribute(PERIMETER_LINE); n2.Attribute(PERIMETER_LINE); n4.Attribute(PERIMETER_LINE);
  _test( nmf.Classify() == ManifoldType::STAND_ALONE );

  // only 2 nodes in manifold
  nmf.Remove( &n4 );
  assert( nmf.Branches() == 2 );

  n1.Attribute(INTERIOR_LINE); n2.Attribute(INTERIOR_LINE);
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY );

  n1.Attribute(INTERIOR_POINT); n2.Attribute(INTERIOR_POINT);
  _test( nmf.Classify() == ManifoldType::STAND_ALONE );

  // SB termination against model exterior
  n1.Attribute(EXTERIOR_SURFACE); n2.Attribute(EXTERIOR_SURFACE);
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );

  n1.Attribute(EXTERIOR_LINE); n2.Attribute(EXTERIOR_LINE);
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );
  
  n1.Attribute(EXTERIOR_POINT); n2.Attribute(EXTERIOR_POINT);
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_END );
  
  _test( nmf.AreNodesCollocated() == true );
  //         ^^^^^^^^^^^^^^^^^^
   
  n1.Coordinate( Point<3>(1.,1.,1.) );
  _test( nmf.AreNodesCollocated() == false );


// CONSTRUCTING 4-NODE MANIFOLD
 
  // constructing from plf_colony of nodes
  // -------------------------------------
  vector<size_t> manifold_nodes{ 0, 1, 2 };
  plf::colony<Node<3>> nodes;
  nodes.insert( n1 );
  nodes.insert( n2 );
  nodes.insert( n3 );
  
  // colony-based constructor
  NodeManifold<3> nmf2( nodes, manifold_nodes );
  //              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  nmf2.AssignManifoldToMemberNodes();   // must be called after the manifold has been constructed
  //   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  _test( nmf2.Branches() == 3 );
  _test( nmf2.N(2)->Idx() == 2 );
  //          ^^^^
  // assigning node attributes
  next(nodes.begin(),0)->Attribute(PERIMETER_SURFACE);
  next(nodes.begin(),1)->Attribute(PERIMETER_SURFACE);
  next(nodes.begin(),2)->Attribute(INTERIOR_SURFACE);
  if ( verbose_ ) nmf2.Out();
  //                   ^^^
  // DIAGNOSTICS                         three-surface nodes that are not interconnected
  _test( nmf2.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );
  //          ^^^^^^^^^^

  // DATA OUTPUT
  /// outputs manifold state to data structure used to initialise VData
  pair<vector<size_t>,ManifoldType> node_manifold_data = nmf2.Data();
  //                                                          ^^^^
  _test( node_manifold_data.second == nmf2.Classify() );

  // TODO: test SortByVariableValue()



} // end run

} // csmp



