//
//  splitCornerTetrahedron.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 4/11/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "splitCornerTetrahedron.hpp"
#include "VSet.h"

namespace csmp {

/**
Assuming the following numbering:

First tetra
0 1 2 3 - ony one neighbor = nbr

Second tetra
0 3 2 4 -  4 neigbors: 0=, 1=, 2=, 3=cnr

The new elements are formed
0 1 2 4 - replacing cnr,             new neighbors: 0=new nbor, 1=nbor1, 2=new, 3=outside
0 3 1 4 - new,                           new neighbors: 0=new nbor, 1=cnr, 2=nbor2, 3=outside
1 3 2 4 - replacing neighbor,    new neighbors: 0=nbor0, 1=nbor1, 2=new, 3=outside

The corresponding neighbor elements are.

*/
void splitCornerTetrahedron( VSet<3>& vset, size_t cnr, size_t nbr )
 {
    assert( cnr < vset.Elements() );
    assert( nbr < vset.Elements() );
    
    // 1. creation of space for new element (assuming that there are no Faces of InterFaces)
    // -------------------------------------------------------------------------------------
    assert( vset.Faces() == 0 );
    assert( vset.InterFaces() == 0 );
    const size_t n_elements_new{ vset.Elements() + 1 };
    if ( vset.HybridElementTypeMesh() ) {
         vset.ResizeElementTypes( n_elements_new );
         vset.ElementType( n_elements_new-1U, vset.ElementType( cnr ) );
      }
    vset.ResizePlist( n_elements_new );
    vset.ResizePlist( n_elements_new, 4 ); // nodes of tetrahedron
    vset.ResizePfverts( n_elements_new );
    vset.ResizePfverts( n_elements_new, 4 ); // nbors of tetrahedron
    
    // 2. assignment of nodes
    // ----------------------
    // new element first
    vset.Plist( n_elements_new-1U, 0, vset.Plist( cnr, 0 ) );
    vset.Plist( n_elements_new-1U, 1, vset.Plist( cnr, 3 ) );
    vset.Plist( n_elements_new-1U, 2, vset.Plist( cnr, 1 ) );
    vset.Plist( n_elements_new-1U, 3, vset.Plist( nbr, 3 ) );
    // temp numbers
    const size_t cnr_n0{vset.Plist(cnr,0)}, cnr_n1{vset.Plist(cnr,3)},
                 cnr_n2{vset.Plist(cnr,1)}, cnr_n3{vset.Plist(nbr,3)},
                 nbr_n0{vset.Plist(cnr,1)}, nbr_n1{vset.Plist(cnr,3)},
                 nbr_n2{vset.Plist(cnr,2)}, nbr_n3{vset.Plist(nbr,3)};
    // reshaped corner element
    vset.Plist( cnr, 0, cnr_n0 );
    vset.Plist( cnr, 1, cnr_n1 );
    vset.Plist( cnr, 2, cnr_n2 );
    vset.Plist( cnr, 3, cnr_n3 );
    // reshaped neighbor element
    vset.Plist( nbr, 0, nbr_n0 );
    vset.Plist( nbr, 1, nbr_n1 );
    vset.Plist( nbr, 2, nbr_n2 );
    vset.Plist( nbr, 3, nbr_n3 );

    // 3. assignment of neigbors
    // -------------------------
    // new element first
    vset.Pfvert( n_elements_new-1U, 0, nbr );
    vset.Pfvert( n_elements_new-1U, 1, cnr );
    vset.Pfvert( n_elements_new-1U, 2, vset.Pfvert( nbr, 2 ) );
    vset.Pfvert( n_elements_new-1U, 3, vset.Pfvert( cnr, 3 ) ); // outside
    // temporaries
    const long64 cnr_p0{nbr}, cnr_p1{vset.Pfvert(nbr,1)},
                 cnr_p2{n_elements_new-1}, cnr_p3{vset.Pfvert(nbr,3)},
                 nbr_p0{vset.Pfvert(nbr,0)}, nbr_p1{vset.Pfvert(nbr,1)},
                 nbr_p2{n_elements_new-1}, nbr_p3{vset.Pfvert(cnr,3)};
    // reshaped corner element
    // new neighbors: 0=new nbor, 1=nbor1, 2=new, 3=outside
    vset.Pfvert( cnr, 0, cnr_p0 );
    vset.Pfvert( cnr, 1, cnr_p1 );
    vset.Pfvert( cnr, 2, cnr_p2 );
    vset.Pfvert( cnr, 3, cnr_p3 ); // outside
    // reshaped neighbor element
    // new neighbors: 0=nbor0, 1=nbor1, 2=new, 3=outside
    vset.Pfvert( nbr, 0, nbr_p0 );
    vset.Pfvert( nbr, 1, nbr_p1 );
    vset.Pfvert( nbr, 2, nbr_p2 );
    vset.Pfvert( nbr, 3, nbr_p3 ); // outside

 } // end splitCornerTetrahedron


} 
