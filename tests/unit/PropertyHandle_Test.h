#ifndef PROPERTY_HANDLE_TEST_H
#define PROPERTY_HANDLE_TEST_H

#include "Test.h"
#include "PropertyHandle.h"
#include "Model.h"


namespace csmp
{

/// test revised by SKM, January 2015
class PropertyHandle_Test : public Test {

    public:
  
      // difference that may be resolved given that tes numbers range between -3000 and 3000
      explicit PropertyHandle_Test( double tolerance=5.0e-11, bool verbose=false );
      
      // delete dynamically allocated model
      ~PropertyHandle_Test();
      
      void run();

    private:
      csmp::Model<3>* model;

      const double TOLERANCE;
  
      // scalar
      PropertyHandle<3> elementVariable1;
      PropertyHandle<3> elementVariable2;

      PropertyHandle<3> nodeVariable1;
      PropertyHandle<3> nodeVariable2;

      PropertyHandle<3> IPVariable1;
      PropertyHandle<3> IPVariable2;

      // vector
      PropertyHandle<3> elementVariable3;
      PropertyHandle<3> elementVariable4;

      PropertyHandle<3> nodeVariable3;
      PropertyHandle<3> nodeVariable4;

      PropertyHandle<3> IPVariable3;
      PropertyHandle<3> IPVariable4;

      // tensor
      PropertyHandle<3> elementVariable5;
      PropertyHandle<3> elementVariable6;

      PropertyHandle<3> nodeVariable5;
      PropertyHandle<3> nodeVariable6;

      PropertyHandle<3> IPVariable5;
      PropertyHandle<3> IPVariable6;
  
      const bool verbose_;
  };
}

#endif // PROPERTY_HANDLE_TEST_H
