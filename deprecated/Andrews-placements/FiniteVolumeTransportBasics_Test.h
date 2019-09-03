#ifndef CSMP_FINITEVOLUMEBASICS_TEST_H
#define CSMP_FINITEVOLUMEBASICS_TEST_H

#include "Test.h"

namespace csmp 
{
    
    template<size_t dim> class Model;

class FiniteVolumeTransportBasics_Test : public Test 
{

  public:
	FiniteVolumeTransportBasics_Test();
	~FiniteVolumeTransportBasics_Test();
	void run(); // runs all the tests for the class (register other methods)

  private:
    void test_constant_velocity_field(Model<3U>& model);
    template<size_t dim> void test_placements(Model<dim>& model);

    void test_b25();
  }; //end class

} //end csmp

#endif
