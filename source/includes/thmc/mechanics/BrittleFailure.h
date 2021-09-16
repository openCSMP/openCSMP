//
//  BrittleFailure.h
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 10/25/14.
//  Copyright (c) 2014 Stephan K. Matthai. All rights reserved.
//

#ifndef BRITTLE_FAILURE_H
#define BRITTLE_FAILURE_H

#include "StressInvariants.h"
#include "MechanicalProperties.h"

namespace csmp {

enum FAILURE : std::int8_t { COMPRESSIVE=-1,
               NONE=0,
               FRICTIONAL_SLIDING=1,
               SHEAR_FRACTURE=2,
               TENSILE_OPENING=3,
               MODE1_FRACTURE=4 };

/// converts enum to string to make it printable
std::string  parseFailure( FAILURE );

/// if the failure mode is not recognized, none is returned
FAILURE parseFailure( const std::string& failure );


class BrittleFailure  {
  public:
    /// evaluates failure criteria -1 to 4 without explicit consideration of pore pressure
    static FAILURE Evaluate( const MechanicalProperties&, const csmp::StressInvariants& );

    /// evaluates failure criteria -1 to 4 with consideration of pore pressure
    static FAILURE Evaluate( const MechanicalProperties&, const csmp::StressInvariants&, double64 pf );
};

} // csmp

#endif /* defined(BRITTLE_FAILURE_H) */
