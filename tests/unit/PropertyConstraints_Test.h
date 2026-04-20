/*
 *  PropertyConstraints_Test.h
 *  csmp_core
 *
 *  Created by SKM 2/2/2022.
 *
 */

#ifndef CSMP_PROPERTY_CONSTRAINTS_TEST_H
#define CSMP_PROPERTY_CONSTRAINTS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

/**
    Unit test for PropertyConstraints object.
    
    @note Dependent functionality 2: needs to be able o work with any kind of CSMP model
*/
class PropertyConstraints_Test: public Test {
  public:
      explicit PropertyConstraints_Test( bool verbose=false );
      ~PropertyConstraints_Test();
      
      virtual void run();
      
  private:
     bool TestBuildRegionsFromPropertyConstraints();
     
     /// tests that 'pointInVolumeElement()'  actually works using 'prism_test'
     bool PointInVolumeElementTest();
  
  private:
    bool verbose_;

};

} // end csmp

#endif /* CSMP_PROPERTY_CONSTRAINTS_TEST_H */
