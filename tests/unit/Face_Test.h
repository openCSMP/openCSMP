// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_FACE_TEST_H
#define CSMP_FACE_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp{

  template<uint32_t> class Face;

  class Face_Test : public Test {
    enum{DIM=3};
    
  public:
    virtual void run();

    // PL2011
    static bool FaceUnitNormalPointsOutward( const Face<3>& f );

  };


} // csmp

#endif // CSMP_FACE_TEST_H
