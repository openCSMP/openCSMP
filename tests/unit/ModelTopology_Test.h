// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_MODEL_TOPOLOGY_TEST_H
#define CSMP_MODEL_TOPOLOGY_TEST_H

#include "Test.h"
#include "ModelTopology.h"

namespace csmp {

/// PL Dec 2010
class ModelTopology_Test : public Test
{
public:
  virtual void run();
  
  static const bool verbose_ = true;
};

} // csmp

#endif // CSMP_MODELTOPOLOGY_TEST_H
