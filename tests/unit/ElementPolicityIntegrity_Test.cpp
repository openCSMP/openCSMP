//
//  ElementPolicityIntegrity_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 9/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "ElementPolicyIntegrity_Test.h"

namespace csmp {

void ElementPolicyIntegrity_Test::run()
 {
   LocalVariables var_set( 1U, // size_t scalarsVars,
                           0U, // size_t vectorVars,
                           0U, // size_t tensorVars,
                           0U, // size_t array_count,
                           0U, // size_t array_length,
                           0U, // size_t flag_array_count,
                           0U, // size_t flag_array_length,
                           1U, // size_t total_data_depth,
                           1U ); // size_t total_flag_depth

   // testing the stand-alone element for a triangle
   ElementOnly<2U>  elmt1( 0, 3, 3 );
   // connecting it up
   Node<2U> nd1( 0, Point<2U>(0.,0.), var_set, IRREGULAR, MESH_VERTEX );
   Node<2U> nd2( 1, Point<2U>(3.,1.), var_set, IRREGULAR, MESH_VERTEX );
   Node<2U> nd3( 2, Point<2U>(1.,2.), var_set, IRREGULAR, MESH_VERTEX );
   // nd1.Assign( 0 /* parent number */, elmt1 );
   // nodes
   elmt1.Assign( 0, &nd1 );
   elmt1.Assign( 1, &nd2 );
   elmt1.Assign( 2, &nd3 );
   // neighbors, but there are none
   elmt1.Assign( 0, &elmt1 );
   
     
 } // end run

// explicit instantiations
template class ElementOnly<2U>;

} // end csmp
