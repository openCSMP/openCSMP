//
//  vsetMakers_Test.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 29/5/2024.
//  Copyright © 2024 The University of Melbourne. All rights reserved.
//

#include "vsetMakers_Test.h"
#include "vsetMakers.h"
#include "CSMP_definitions.h"
#include "CSMP_ElementSpecifications.h"
#include "VSet.h"
#include "ModelTopology.h"

using namespace std;

namespace csmp {

// auxiliary functions
template<uint32_t dim>
vector<vector<uint32_t>> vsetMakers_Test::NodesOfElementFaces( const VSet<dim>& vset, size_t elmt )
 {
    assert( elmt < vset.Elements() );
    assert( vset.HybridElementTypeMesh() );

    const auto n_nbors = distance( vset.PfvertsBegin(elmt), vset.PfvertsEnd(elmt) );
    const auto e_type  = vset.ElementType( elmt );

    vector<vector<uint32_t>> elmt_face_nodes( n_nbors );
    cout <<"\nprintNeighboursOfElement: element "<< elmt << endl;
    cout <<"\t"<<"faces and their nodes:"<< endl;
    for ( uint32_t face{0u}; face<n_nbors; ++face ) {
          const uint32_t nodes_per_face = CSMP_ElementSpecifications::NodesPerFaceForElementOfType( e_type, face );
          vector<uint32_t> fnids( nodes_per_face );
          for ( uint32_t n{0U}; n<nodes_per_face; ++n )
            fnids[n] = CSMP_ElementSpecifications::FaceNodeForElementOfType( e_type, face, n );
          // printing the relevant information
          cout <<"\t\t"<< face <<": ";
          for ( const auto& i : fnids ) cout <<" "<< vset.Plist( elmt, i );
          cout << endl;
          elmt_face_nodes[face] = fnids;
      }
    
   return elmt_face_nodes;
    
 } // end NodesOfElementFaces
 
template vector<vector<uint32_t>> vsetMakers_Test::NodesOfElementFaces( const VSet<3>&, size_t );
template vector<vector<uint32_t>> vsetMakers_Test::NodesOfElementFaces( const VSet<2>&, size_t );





void vsetMakers_Test::run()
 {
    // element neighbor connectivity
    {
       const bool bSkewed{false};
       VSet<3U> vset;
       
       create_Tetra_VSet( vset );
       TestConnectedElementNeighborNumbering( vset );

       create_Prism_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       create_Pyramid_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       create_RubikCube_VSet( vset );
       TestConnectedElementNeighborNumbering( vset );

       create_Pyramid_Hexa_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       create_Prism_Hexa_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       ModelTopology topo;
       // removing properties that were added to the VSet before
       auto pit = vset.PropertyValuesBegin();
       while ( pit!=vset.PropertyValuesEnd() ) pit = vset.RemoveData( (*pit).first.c_str() );
       create_FracBox( topo, vset );
       TestConnectedElementNeighborNumbering( vset );
    }

    // 2D tests
    {
      VSet<2U> vset;
      // 2D test cases for all element types except simplicies
      ///  Rectangle-shaped MODEL_TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object marking the box boundary.
      const double length_of_sides{5.};
      const bool bSkewed{false};
      create_1Square_VSet( vset, length_of_sides, bSkewed );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );

      /// Rectangle shaped mixed model with 2 intersecting line element regions (USED IN UNIT TESTS)
      ModelTopology topo = create_MeshPatchWithLineElements_VSet( vset );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );
       TestConnectedElementNeighborNumbering( vset );

      /// model SPLIT22_BASIC with box boundaries (Faces) and one through-going and one internal crossing split boundary (USED IN UNIT TESTS)
      topo = create_BoundarySplitBoundaryPatch( vset );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );
       TestConnectedElementNeighborNumbering( vset );
    }
    
    // 3D test cases
    {
      VSet<3U> vset;
      const bool bSkewed{false};
      create_1Hexahedron_VSet( vset );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );
      create_1Prism_VSet( vset, bSkewed );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );
    }

 } // end run


/**
     For single element meshes only.
*/
bool vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<2U>& vset )
 {
    assert( vset.Elements() == 1 );
    size_t inconsistencies_found{0ul};
    
    vector<size_t> test_elements{ 0, 5, 12, 17 };
    for ( auto eid : test_elements ) {
         auto face_nodes = NodesOfElementFaces( vset, eid );
      }
      
    return (inconsistencies_found == 0);

} // end TestConsistencyWithCSMP_Conventions


/**
      Testing numbering of neighbor faces 
*/
bool vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<3U>& vset )
 {
    assert( vset.Elements() == 1 );
    size_t inconsistencies_found{0ul};
    
    vector<size_t> test_elements{ 0, 5, 12, 17 };
    for ( auto eid : test_elements ) {
         auto face_nodes = NodesOfElementFaces( vset, eid );
      }
      
    return (inconsistencies_found == 0);

} // end TestConsistencyWithCSMP_Conventions



    /// for the supplied vset, method checks whether the assigned neighbor Element objects are contained in the overall element range
template<uint32_t dim>
void vsetMakers_Test::TestConnectedElementNeighborNumbering(  const VSet<dim>& vset )
 {
    const size_t n_elements = vset.Elements();
    for ( size_t eidx{0ul}; eidx<n_elements; ++eidx )
      for ( uint32_t i{0u}; i<vset.PfvertsSize(eidx); ++i ) {
           auto neighbor = vset.Pfvert( eidx, i );
           // testing consistency of element type
           // testing consistency of neighbors
           if ( neighbor >= 0 ) {
                _test( neighbor < n_elements );
                if ( vset.HybridElementTypeMesh() ) {
                     _test( vset.PfvertsSize(neighbor) == CSMP_ElementSpecifications::NeighborsPerElementOfType( vset.ElementType(neighbor) ) );
                     _test( vset.PlistSize(neighbor) == CSMP_ElementSpecifications::NodesPerElementOfType( vset.ElementType(neighbor) ) );
                  }
             }
        }
 
 } // end TestConnectedElementNeighbors

template void vsetMakers_Test::TestConnectedElementNeighborNumbering(  const VSet<2>& );
template void vsetMakers_Test::TestConnectedElementNeighborNumbering(  const VSet<3>& );


} // end csmp
