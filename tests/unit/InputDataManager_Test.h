/*
 *  InputDataManager_Test.h
 *  csmp_core
 *
 *  Created by Ali Tabatabaei on 11/4/10.
    refactored by SKM 2016.
 *
 */

#ifndef INPUT_DATA_MANAGER_TEST_H
#define INPUT_DATA_MANAGER_TEST_H

#include "Test.h"

namespace csmp {

class ANSYS_Model3D;

class InputDataManager_Test: public Test {
  public:
      explicit InputDataManager_Test( bool verbose=false );
      ~InputDataManager_Test();
      
      virtual void run();
  
  private:
    bool verbose_;

};

} // end csmp

#endif /* INPUT_DATA_MANAGER_TEST_H */
