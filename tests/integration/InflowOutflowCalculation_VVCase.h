#ifndef NCFVTA_VVCASE_H
#define NCFVTA_VVCASE_H

#include <iostream>
#include <iomanip>
#include <limits>
#include "CSMP_number_types.h"
#include "Test.h"


namespace csmp{

//template<size_t> class Model;
template<size_t dim>
class InflowOutflowCalculation_VVCase : public Test
  {
    public:
        explicit InflowOutflowCalculation_VVCase();
        explicit InflowOutflowCalculation_VVCase( const char* prefix );
        ~InflowOutflowCalculation_VVCase();

        void run();

    private:
      const char* prefix_;
      string name_;
  };

} // csmp


#endif // NCFVTA_VVCASE_H
