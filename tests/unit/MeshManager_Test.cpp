//
//  MeshManager_Test.cpp
//  CSMP_GitHub_UnitTests-Intel
//
//  Created by Stephan Matthai on 23/08/2018.
//  Copyright © 2018 Stephan Matthai. All rights reserved.
//

#include "MeshManager_Test.h"

using namespace std;

namespace csmp {

void MeshManager_Test::run()
 {
    _test( TestEntityNumberingFunction() );
  
    _test( TestElementDeletionAndInsertion() );
  
    _test( TestFaceDeletionAndInsertion() );
  
 } // end run
 
 
/**
    Tests whether the elements are correctly numbered consecutively by
    the mesh manager.
*/
bool MeshManager_Test::TestEntityNumberingFunction()
  {
     return true;
  }

} // end csmp
