#ifndef CSMP_ELEMENT_TEST_H
#define CSMP_ELEMENT_TEST_H

#include "Element.h"
#include "Test.h"

namespace csmp
{

class Element_Test : public Test
{

  public:
	  Element_Test();
  
    virtual void run();
    
    /// using a quadrilateral element for testing 
    void ElementLengthTest2D();
    void ElementLengthTest3D();
    
    /// is it correctly located
    void FaceBaryCenterTest();

    /// is it correctly located
    void BaryCenterTest();
    
    void MoveSemanticsTest();
    
    /// tests that all members of the supplied elements are the same
    template<uint32_t dim>
    void CompareElements( const Element<dim>&, const Element<dim>& );
    
    void VariableAccessAndIterators();
  
  private:
    double fTolerance   = 1.0e-3;
    const bool verbose_ = true;
};

} //end csmp

#endif
