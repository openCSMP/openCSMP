#ifndef SNEDDON_CRACK_COUPLED_3D_VVCASE_H
#define SNEDDON_CRACK_COUPLED_3D_VVCASE_H

#include "Test.h"

namespace csmp
  {

  class SneddonCrackCoupled3D_VVCase: public Test
    {
    public:
      SneddonCrackCoupled3D_VVCase(const char* prefix);

      void set_min_avg_error(double epsilon){ min_avg_error = epsilon; };

      virtual void run();

    private:

     double min_avg_error = 8.0;


  };




  } // csmp

#endif
