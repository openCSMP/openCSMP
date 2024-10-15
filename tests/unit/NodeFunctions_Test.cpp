//
//  NodeFunctions_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "NodeFunctions_Test.h"
#include "CSMP_definitions.h"
#include "meshManagementUtilities.h"
#include "vsetMakers.h"
#include "VTU_Interface.h"
#include "Region.h"
#include "Element.h"
#include "compareFloats.h"

#include "ModelTopology.h"
#include "NodeManifold.h"


using namespace std;

namespace csmp {

/**
custom models
uses "CSMP-1phase-variables.txt"
*/
NodeFunctions_Test::NodeFunctions_Test()
{
}



template<uint32_t dim>
bool NodeFunctions_Test::Test_nodeNeighbors( const Model<dim>& model )
 {
    const Region<dim>&  model_domain(model.Region("Model"));
    
    vector<set<size_t>> node_neighbors;
    nodeNeighbors( model_domain, node_neighbors, verbose_ );
    
    // determining typical number of node neighbors in mesh
    size_t n_neighbors{0};
    for ( auto nit : node_neighbors ) n_neighbors += nit.size();
    cout <<"\n\taverage number of neighbors per node: "<< n_neighbors / node_neighbors.size();
    
    return ( n_neighbors > 0ul );
}




/**
THE MASTER TEST FUNCTION
*/
void NodeFunctions_Test::run()
{
  _test(Test_parentElementsSharedByFace()); // OK

  // node neighbors
  {
    VSet<3>    vset;
    const bool bSkewed{false};
    create_Pyramid_Hexa_VSet( vset, bSkewed );
    Model<3> model( vset );
    model.Name("PyramidHexaPatch");

    _test(Test_nodeNeighbors(model) );
  }
  
  cout << endl;

} // end run







// tests method with the same name
bool NodeFunctions_Test::Test_parentElementsSharedByFace()
 {
   // 1. poly-element mesh test case 'Prism_Hexa'
   {
      VSet<3U> vset;
      // testing with element 13 with face 4 on the LEFT outside
      create_Prism_Hexa_VSet( vset, false );
      Model<3U>    model( vset, "CSMP-variables.txt" );
      const size_t ELMT{13}; // 13 in VSet
      Element<3>*  eptr = &(*next(model.Mesh().ElementsBegin(),ELMT));
       
      // getting an inner face in the 3D model that is not on the boundary
      vector<Node<3>*> face_nodes;
      Element<3>*      inner_eptr(nullptr);
      for ( uint32_t i{0u}; i<eptr->Neighbors(); ++i )
       if ( eptr->Neighbor(i) != nullptr ) {
            face_nodes.reserve( eptr->FE()->NodesPerFace(i) );
            for ( const auto& j : eptr->FE()->NodesOfFace(i) )
              face_nodes.push_back( eptr->N(j) );
            inner_eptr = eptr;
            break;
         }
      
      // calling the function that is being tested
      auto parents = parentElements<3>( face_nodes.begin(), face_nodes.end() );
      
      // test that the correct neighbor elements were found (inner one should be first)
      _test( parents.first.first  != nullptr );
      _test( parents.second.first != nullptr );
      _test( parents.first.first  == inner_eptr );
      _test( parents.second.first->Idx() == 3u ); // element 3 shares the face
    }
    
    // 2. test case with a single hexahedron
    {
      VSet<3U> vset;
      // testing with element 13 with face 4 on the LEFT outside
      create_1Hexahedron_VSet( vset, false );
      Model<3U>    model( vset, "CSMP-variables.txt" );
      Element<3>*  eptr = &(*model.Mesh().ElementsBegin());
      
      // testing for a face of the element that lies on outside of the model
      Element<3>*      inner_eptr(nullptr);
      vector<Node<3>*> face_nodes;
      for ( uint32_t i{0u}; i<eptr->Neighbors(); ++i )
       if ( eptr->Neighbor(i) == nullptr ) {
            face_nodes.reserve( eptr->FE()->NodesPerFace(i) );
            for ( const auto& j : eptr->FE()->NodesOfFace(i) )
              face_nodes.push_back( eptr->N(j) );
            inner_eptr = eptr;
            break;
         }

      // calling the function that is being tested
      auto parents = parentElements<3>( face_nodes.begin(), face_nodes.end() );
      
      // test that the correct neighbor elements were found (inner one should be first
      _test( parents.first.first  == inner_eptr );
      _test( parents.second.first == nullptr );
      _test( parents.first.first  == eptr );
    }
    
    return true;
    
 } // end Test_parentElementsSharedByFace


} // end csmp
