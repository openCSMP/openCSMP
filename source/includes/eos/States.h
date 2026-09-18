// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STATES_H
#define STATES_H

/// enum to characterize possible phase states in the H2O-NaCl phase diagram
namespace csmp
{
  enum States
  {
    none,    // 0 i.e. state undefined!
    L,       // 1 liquid, for convenience only, use rather 'F'
    V,       // 2 vapor, for convenience only, use rather 'F'
    F,       // 3 single phase fluid
    VL,      // 4 vapor+liquid twophase 
    VH,      // 5 vapor+halite twophase
    LH,      // 6 liquid+halite twophase
    VLH,     // 7 vapor+liquid+halite twophase
    H,       // 8 halite
    M,       // 9 NaCl melt
    HM       // 10 Halite+NaCl melt
  };
}
#endif
