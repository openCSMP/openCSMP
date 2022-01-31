/** 
@class LinearCuboid_Test
@author Hani Akbari (Jan. 2017)
This class tests implementation of LinearCuboid element.

1- In constructor, LinearCuboid_Test(bool) we define 8 nodes of a cuboid element 
in standard oredering and each node has 3 coordinates. 

2- Run() tests integral of NN and dNdN i.e. integrates multilplication of shape functions and 
gradient of shape functions over element and prints the result.

3- Run() tests values of shape functions at nodes : TestSumShapesAtBaryCenter(const Element<3U>&). 
It also tests sum of values of shape functions at barycenter that must be one :
TestSumShapesAtBaryCenter(const Element<3U>& )

4- _equalTest(a,b,tol) returns false if |a-b| > tol.

*/
#ifndef LINEARCUBOID_TEST_H
#define LINEARCUBOID_TEST_H

#include "Test.h"
#include "LinearCuboid.h"
#include "Element.h"

namespace csmp {

	class LinearCuboid_Test : public Test {
	public:
		explicit LinearCuboid_Test(bool verbose);
		~LinearCuboid_Test();
		virtual void run(); // It tests integrals of NN and dNdN as well.
		void TestInterpolationFunctionValues(const Element<3U>& e);
		void TestSumShapesAtBaryCenter(const Element<3U>& e);
		bool _equalTest(double a, double b, double tol) const;

	private:
		LinearCuboid    lcuboid_;
		Element<3U>     element_;
		Node<3U>        n0, n1, n2, n3, n4, n5, n6, n7;
		const double  tolerance_factor_;
		bool            verbose_;
	};

}
#endif

