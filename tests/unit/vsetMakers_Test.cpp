//
//  vsetMakers_Test.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 29/5/2024.
//  Copyright © 2024 The University of Melbourne. All rights reserved.
//

#include "vsetMakers_Test.h"

using namespace std;

namespace csmp {

// auxiliary functions
vector<set<size_t>> vsetMakers_Test::NodesOfElementFaces( const VSet<3U>& vset, size_t elmt )
 {
    assert( elmt < vset.Elements() );
    assert( vset.HybridElementTypeMesh() );

    const auto n_nbors = distance( vset.PfvertsBegin(elmt), vset.PfvertsEnd(elmt) );
    const auto e_type  = vset.ElementType( elmt );

    vector<set<size_t>> elmt_face_nodes( n_nbors );
    cout <<"\nprintNeighboursOfElement: element "<< elmt << endl;
    cout <<"\t"<<"faces and their nodes:"<< endl;
    for ( uint32_t face{0u}; face<n_nbors; ++face ) {
          const uint32_t nodes_per_face = CSMP_ElementSpecifications::NodesPerFaceForElementOfType( e_type, face );
          set<uint32_t> fnids( nodes_per_face );
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
 


vsetMakers_Test::run()
 {
    // 2D test cases for all element types except simplicies
    _test( testConsistencyWithCSMP_Conventions( const VSet<2U>& ) );
    
    
    // 3D test cases
    _test( testConsistencyWithCSMP_Conventions( const VSet<3U>& ) );

 } // end run


/**
     For single element meshes only.
*/
vsetMakers_Test::TestConsistencyWithCSMP_Conventions( const VSet<2U>& )
 {
    assert( vset.Elements() == 1 );
 
    std::vector<std::set<size_t>> NodesOfElementFaces( const VSet<3U>& vset, size_t elmt );
} // end TestConsistencyWithCSMP_Conventions


} // end csmp
