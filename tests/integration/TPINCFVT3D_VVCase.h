#ifndef TPINCFVT3D_VVCASE_H
#define TPINCFVT3D_VVCASE_H

#include <iostream>
#include <iomanip>
#include <limits>
#include "Test.h"


namespace csmp{

template<uint32_t> class Model;
struct TPINCFVT3D_TestData;

class TPINCFVT3D_VVCase : public Test
{
  public:
      explicit TPINCFVT3D_VVCase();
      explicit TPINCFVT3D_VVCase( const char* prefix );
      ~TPINCFVT3D_VVCase();

      void run();

      /// Set model
      void SetModel( Model<3>& model,
                     double mobilityRatio,
                     double velX,
                     double porosity );

      void OutputVTU( Model<3>& model,
                   const char* fileName,
                   long time );

      double TestFront( Model<3>& model,
                        double tolSat,
                        double satOil,
                        double tolFront,
                        double analyticFront
                      );

      void Test( TPINCFVT3D_TestData& data, double tolSat, double satOil, double tolFront );

  private:
      const char* prefix_;
      std::string name_;
};


} // csmp


#endif // TPINCFVT3D_VVCASE_H
