#ifndef NCFVTA_VVCASE_H
#define NCFVTA_VVCASE_H

#include "Test.h"

namespace csmp {

//template<uint32_t> class Model;
template<uint32_t dim>
class InflowOutflowCalculation_VVCase : public Test
  {
    public:
        explicit InflowOutflowCalculation_VVCase();
        explicit InflowOutflowCalculation_VVCase( const char* prefix );
        ~InflowOutflowCalculation_VVCase();

        void run();

    private:
      const char* prefix_;
      std::string name_;
  };

} // csmp


#endif // NCFVTA_VVCASE_H
