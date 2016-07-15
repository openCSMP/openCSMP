#ifndef CSP_ELEMENT_TEST_H
#define CSP_ELEMENT_TEST_H

#include "Element.h"
#include "Test.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "Node.h"

namespace csmp
{

template <> class Element<3U>;

class Element_Test : public Test 
{

  public:
	Element_Test();
	~Element_Test();
	void run(); // runs all the tests for the class (register other methods)
	void lengthMeasurement2D_Test();
	//void lengthMeasurement3D_Test();
	
  
  private:
  double64 fTolerance;
  
}; //end class

} //end csp

#endif
