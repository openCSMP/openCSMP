//
//  GenericFiniteVolumeTransport_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/01/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#ifndef GenericFiniteVolumeTransport_Test_h
#define GenericFiniteVolumeTransport_Test_h

#include "Test.h"

namespace csmp {

class GenericFiniteVolumeTransport_Test : public Test {
  public:
    virtual void run();
  
    void BenchmarkGlobalVersusParametricIntegration();
};



} // end csmp


#endif /* GenericFiniteVolumeTransport_Test_h */
