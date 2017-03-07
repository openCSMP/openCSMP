#ifndef CSMP_ELEMENT_TEST_H
#define CSMP_ELEMENT_TEST_H

#include "Element.h"
#include "Test.h"

namespace csmp
{

class Element_Test : public Test
{

  public:
	explicit Element_Test( bool verbose=false );
	~Element_Test();
  
	virtual void run();
  
  /// using a quadrilateral element for testing 
	void ElementLengthTest2D();
	void ElementLengthTest3D();
  
  /// is it correctly located
  void FaceBaryCenterTest();

  /// is it correctly located
  void BaryCenterTest();
  
  /// for all element types tests wether the face normals are outward pointing
  void UnitNormalTest();
	
  
  private:
    double64 fTolerance;
    const bool verbose_;
};

} //end csmp

#endif
