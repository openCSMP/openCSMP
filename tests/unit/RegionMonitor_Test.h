#ifndef CSMP_REGION_MONITOR_TEST_H
#define CSMP_REGION_MONITOR_TEST_H

#include "RegionMonitor.h"
#include "Test.h"

namespace csmp 
{

class RegionMonitor_Test : public Test 
{

  public:
	RegionMonitor_Test();
	~RegionMonitor_Test();
	void run(); // runs all the tests for the class (register other methods)
	void RegionMonitorScalarPropertyIntegrals2D();
	void RegionMonitorScalarPropertyIntegrals3D();
	
  private:
  double fTolerance;
  
}; //end class

} //end csmp

#endif
