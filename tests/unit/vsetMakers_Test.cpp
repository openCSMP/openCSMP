//
//  vsetMakers_Test.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 29/5/2024.
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

    const auto n_nbors = distance( vset.PfvertsBegin(elmt), vset.PfvertsEnd(elmt) );
    const auto e_type  = vset.ElementType( elmt );

    vector<vector<uint32_t>> elmt_face_nodes( static_cast<size_t>(n_nbors) );
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

       create_Hexahedra_VSet( vset );
       TestConnectedElementNeighborNumbering( vset );

       create_Pyramid_Hexa_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       create_Prism_Hexa_VSet( vset, bSkewed );
       TestConnectedElementNeighborNumbering( vset );

       
       // removing properties that were added to the VSet before
       auto pit = vset.PropertyValuesBegin();
       while ( pit!=vset.PropertyValuesEnd() ) pit = vset.RemoveData( (*pit).first.c_str() );
       ModelTopology topo = create_FracBox( vset );
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
       auto pit = vset.PropertyValuesBegin();
       while ( pit!=vset.PropertyValuesEnd() ) pit = vset.RemoveData( (*pit).first.c_str() );

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
       auto pit = vset.PropertyValuesBegin();
       while ( pit!=vset.PropertyValuesEnd() ) pit = vset.RemoveData( (*pit).first.c_str() );

      create_1Prism_VSet( vset, bSkewed );
      _test( TestConsistencyWithCSMP_Conventions( vset ) );
    }

 } // end run


/**
     For single element meshes only.
*/
template<uint32_t dim>
bool vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<dim>& vset )
 {
    size_t inconsistencies_found{0ul};
    
    // do the neighbor types of the element match
    for (size_t eidx = 0; eidx < vset.Elements(); ++eidx) {
        for (uint32_t j = 0; j < vset.PfvertsSize(eidx); ++j) {
            auto neighbor_idx = vset.Pfvert(eidx, j);
            // If the face has a neighbor
            if (neighbor_idx >= 0 ) {  // numbers below zero indicate no neighbor
                // Check the neighbors of this neighbor
              for ( uint32_t k{0u}; k < vset.PfvertsSize(static_cast<size_t>(neighbor_idx)); ++k) {
                    auto neighbor_of_neighbor_idx = vset.Pfvert(static_cast<size_t>(neighbor_idx), k);

                    // If neighbor_of_neighbor is valid and connected to the current element
                    if ( neighbor_of_neighbor_idx >= 0 && neighbor_of_neighbor_idx == static_cast<int64_t>(eidx) ) {
                        auto etype      = ( vset.HybridElementTypeMesh() == true ) ? vset.ElementType(static_cast<size_t>(neighbor_idx)) : vset.ElementType(0);
                        auto nbor_etype = ( vset.HybridElementTypeMesh() == true ) ? vset.ElementType(static_cast<size_t>(neighbor_of_neighbor_idx)) : vset.ElementType(0);
                        if ( verbose_ )
                             cout << "\nNeighbor " << neighbor_idx << ": "
                                  << parseAbbreviated_FE_Type(etype)
                                  << " is connected to "
                                  << neighbor_of_neighbor_idx << ": "
                                  << parseAbbreviated_FE_Type(nbor_etype);
                                  
                        _test( parseFiniteElementDimension(static_cast<CSMP_FEM_TYPE>(etype)) ==
                               parseFiniteElementDimension(static_cast<CSMP_FEM_TYPE>(nbor_etype)) );
                        if ( parseFiniteElementDimension(static_cast<CSMP_FEM_TYPE>(etype)) !=
                             parseFiniteElementDimension(static_cast<CSMP_FEM_TYPE>(nbor_etype)) )
                          inconsistencies_found++;
                    }
                }
            }
        }
    }
      
    return (inconsistencies_found == 0);

} // end TestConsistencyWithCSMP_Conventions

template bool vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<3>& );
template bool vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<2>& );





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
                _test( neighbor < static_cast<int64_t>(n_elements) );
                if ( vset.HybridElementTypeMesh() ) {
                  _test( vset.PfvertsSize(static_cast<uint32_t>(neighbor)) == CSMP_ElementSpecifications::NeighborsPerElementOfType( vset.ElementType(static_cast<uint32_t>(neighbor)) ) );
                  _test( vset.PlistSize(static_cast<uint32_t>(neighbor)) == CSMP_ElementSpecifications::NodesPerElementOfType( vset.ElementType(static_cast<uint32_t>(neighbor)) ) );
                  }
             }
        }
 
 } // end TestConnectedElementNeighbors

template void vsetMakers_Test::TestConnectedElementNeighborNumbering(  const VSet<2>& );
template void vsetMakers_Test::TestConnectedElementNeighborNumbering(  const VSet<3>& );


} // end csmp
