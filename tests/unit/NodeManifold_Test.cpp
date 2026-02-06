#include "NodeManifold_Test.h"
#include "Node.h"
#include "NodeManifold.h"
#include "plf_colony.h"

using namespace std;

namespace csmp{

void NodeManifold_Test::run()
{
  LocalVariables nvars;

  // CONSTRUCTOR arguments: size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY = NOT, TOPOTYPE = MESH_VERTEX );
  Node<3> n1( 0, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n2( 1, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n3( 2, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );
  Node<3> n4( 3, Point<3>(3.,3.,3.), nvars, INTERNAL, INTERIOR_SURFACE );


  // NODE NEIGHBORS
  set<Node<3>*> nbor_nodes_n1{&n2,&n3}, nbor_nodes_n2{&n1,&n3}, nbor_nodes_n3{&n1,&n2};
  n1.Assign( nbor_nodes_n1 );
  n2.Assign( nbor_nodes_n2 );
  n3.Assign( nbor_nodes_n3 );
  
  
  // CONSTRUCTING MANIFOLDS
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
  for ( size_t i{0}; i<nmf.Branches(); ++i )
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
  
 
  // constructing from plf_colony of nodes
  // -------------------------------------
  vector<size_t> manifold_nodes{ 0, 1, 2, 3 };
  plf::colony<Node<3>> nodes;
  nodes.insert( n1 );
  nodes.insert( n2 );
  nodes.insert( n3 );
  nodes.insert( n4 );
  
  NodeManifold<3> nmf2( nodes, manifold_nodes );
  //              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  nmf2.AssignManifoldToMemberNodes();
  //   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  _test( nmf2.Branches() == 4 );
  _test( nmf2.N(2)->Idx() == 2 );
  //          ^^^^
  if ( verbose_ ) nmf.Out();
  //                  ^^^


  // DIAGNOSTICS                         three-surface nodes that are not interconnected
  _test( nmf.Classify() == ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH );
  //         ^^^^^^^^^^
  
  // remaining nodes 0,1,3
  _test( nmf.N(0)->Idx() == 0 );
  _test( nmf.N(1)->Idx() == 1 );
  _test( nmf.N(2)->Idx() == 3 );
  
  // other configurations for classification
  n1.Attribute(PERIMETER_LINE); n2.Attribute(PERIMETER_LINE); n4.Attribute(PERIMETER_LINE);
  _test( nmf.Classify() == ManifoldType::STAND_ALONE );

  // only 2 nodes in manifold
  nmf.Remove( &n4 );

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
  _test( nmf.Classify() == ManifoldType::STAND_ALONE );
  
  // TODO: SortByVariableValue()

  _test( nmf.AreNodesCollocated() == true );
  //         ^^^^^^^^^^^^^^^^^^
   
  n1.Coordinate( Point<3>(1.,1.,1.) );
  _test( nmf.AreNodesCollocated() == false );
   
  // DATA OUTPUT
  /// outputs manifold state to data structure used to initialise VData
  pair<vector<size_t>,ManifoldType> node_manifold_data = nmf2.Data();
  //                                                          ^^^^
  _test( node_manifold_data.second == nmf2.Classify() );

}

} // csmp
