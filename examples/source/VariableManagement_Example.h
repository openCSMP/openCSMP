//
//  VariableManagement_Example.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 27/9/2022.
//

#include "Example.h"

namespace csmp {

/**
    Shows how to selectively load variables from CSMP binary files and store
    them back again.
*/
class VariableManagement_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
};

} // csmp
