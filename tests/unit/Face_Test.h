#ifndef FACE_TEST_H
#define FACE_TEST_H

#include "Test.h"

namespace csmp{

  template<size_t> class Face;

  // PL2011
  class Face_Test : public Test
  {
    enum{DIM=3};
  public:
    virtual void run();

    static bool FaceUnitNormalPointsOutward( const Face<3>& f );

  };


} // csmp

#endif // FACE_TEST_H
