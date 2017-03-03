#ifndef CSMP_REGION_MONITOR_TEST_H
#define CSMP_REGION_MONITOR_TEST_H

#include "RegionMonitor.h"
#include "Test.h"
#include "CSMP_number_types.h"

namespace csmp 
{

class RegionMonitorTest : public Test 
{

  public:
	RegionMonitorTest();
	~RegionMonitorTest();
	void run(); // runs all the tests for the class (register other methods)
	void RegionMonitorScalarPropertyIntegrals2D();
	void RegionMonitorScalarPropertyIntegrals3D();
	
  private:
  double64 fTolerance;
  
}; //end class

} //end csmp

#endif