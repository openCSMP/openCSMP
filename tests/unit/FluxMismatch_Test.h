#ifndef CSMP_FLUX_MISMATCH_TEST_H
#define CSMP_FLUX_MISMATCH_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

class FluxMismatch_Test : public Test
{
  public:
    void run() override final;
    
  private:
    int FluxMismatch( bool bPrescribedVelocity );
    
    constexpr static bool verbose_ = true;
};

} // csmp

#endif
