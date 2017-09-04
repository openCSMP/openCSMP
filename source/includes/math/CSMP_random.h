#ifndef CSMP_RANDOM_H
#define CSMP_RANDOM_H

#include "CSMP_definitions.h"

namespace csmp {

/*! file CSMP_random.h */

/// random number generator
/**
 * For now we are using the 64-bit Mersenne twister and just borrowing
 * its implementation. Eventually we want to be able to split streams,
 * which will require some additional operations in here.
 */
typedef std::mt19937_64 random_generator;
    

} // csmp

#endif
