#ifndef ANSYS_MODEL3D_TEST_H
#define ANSYS_MODEL3D_TEST_H

#include "Test.h"

namespace csmp {

  // P. Lang 2011
  class ANSYS_Model3D_Test : public Test
    {
    public:
      virtual void run();
      
      void ModelRecoveryFromFileTest();
      
    private:
      static const bool verbose_ = false;
    };

   void create_ANSYS3D_Model( bool contiguous, bool reconstruct_from_file );


  } // csmp

#endif // ANSYS_MODEL3D_TEST_H
