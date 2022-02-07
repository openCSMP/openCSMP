#ifndef CSMP_COLONY_WITH_ELEMENT_TEST_H
#define CSMP_COLONY_WITH_ELEMENT_TEST_H

#include "plf_colony.h"
#include "Element.h"
#include "Test.h"

namespace csmp
{

/**
    Tests how plf::colony manages csmp Element objects.
*/
class Colony_Test : public Test
{

  public:
    Colony_Test();
    ~Colony_Test();
    
    virtual void run();
    
    /// using a quadrilateral element for testing 
    bool TestColonyWith_int();
    bool TestColonyWith_Element2();
  
  private:
    const bool verbose_ = true;
};

} //end csmp

#endif /* CSMP_COLONY_WITH_ELEMENT_TEST_H */
