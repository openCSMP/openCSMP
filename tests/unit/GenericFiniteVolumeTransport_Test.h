//
//  GenericFiniteVolumeTransport_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/01/2017.
//

#ifndef CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_TEST_H
#define CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_TEST_H

#include "Test.h"

namespace csmp {

class GenericFiniteVolumeTransport_Test : public Test {
  public:
    virtual void run();
  
    void TestBasics();
    void BenchmarkGlobalVersusParametricIntegration();
};



} // end csmp


#endif /* CSMP_GENERIC_FINITE_VOLUME_TRANSPORT_TEST_H */
