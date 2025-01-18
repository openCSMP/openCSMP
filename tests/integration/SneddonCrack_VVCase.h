#ifndef SNEDDON_CRACK_VVCASE_H
#define SNEDDON_CRACK_VVCASE_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {
  
  template<uint32_t> class Point;

  class SneddonCrack_VVCase: public Test {
    public:
      SneddonCrack_VVCase(const char* prefix);
      virtual void run();

      void set_min_avg_error(double epsilon){ min_avg_error = epsilon; };


    private:

      double min_avg_error = 8.0;

      template<typename T>
      void WriteSolutionToFile( std::vector<std::vector<T>> data, std::string name);
      
      template<typename T, uint32_t dim>
      void Verify( std::map<Point<dim>,double>& PressureProfile, bool at_angle, T F , double tol );
    };




  } // csmp

#endif
