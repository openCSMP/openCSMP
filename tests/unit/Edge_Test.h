#ifndef CSMP_EDGE_TEST_H
#define CSMP_EDGE_TEST_H

#include "Edge.h"
#include "Test.h"

namespace csmp
{

//template <> class Edge<3U>;

class Edge_Test : public Test 
{

  public:
	Edge_Test();
	~Edge_Test();
	void run(); // runs all the tests for the class (register other methods)
  
}; //end class

} //end csmp

#endif
