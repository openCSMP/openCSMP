#ifndef CSMP_MAPLE_INTERFACE_H
#define CSMP_MAPLE_INTERFACE_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Model;

/**
@file MapleInterface.h
*/

/**
@addtogroup CSMPglobalFunctions
*/

/// creates 1D dataset for plotting; output needs to be pasted into a Maple worksheet
void writeVariableToMapleTextFile( const Model<1U>& sg, 
                                   const char* variable, uint32 timestep, double64 time );

/// creates 1D dataset from a target region for plotting; output needs to be pasted into a Maple worksheet
void writeVariableToMapleTextFile( const Model<1U>& sg, const char* group,
                                   const char* variable, uint32 timestep, double64 time );

/// as above, but for plotting pairs of variables
void writeVariablesToMapleTextFile( const Model<1U>& sg,
                                    const char* variable1, const char* variable2, 
                                    uint32 timestep, double64 time );
/**
  @}
*/

} // end csmp

#endif
