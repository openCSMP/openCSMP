#ifndef FACE_TEST_H
#define FACE_TEST_H

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

#endif // FACE_TEST_H
