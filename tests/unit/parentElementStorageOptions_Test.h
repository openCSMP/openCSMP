//
//  parentElementStorageOptions_Test.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 19/8/2024.
//

#ifndef CSMP_NODE_PARENT_ELEMENT_STORAGE_OPTIONS_TEST_H
#define CSMP_NODE_PARENT_ELEMENT_STORAGE_OPTIONS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

class parentElementStorageOptions_Test : public Test {
  public:
    virtual void run();
    
    const bool verbose_ = true;
    
  private:
    
    void Test_AlternativeContainersForSizeAndSpeed();
};

} // end csmp

#endif /* CSMP_NODE_PARENT_ELEMENT_STORAGE_OPTIONS_TEST_H */
