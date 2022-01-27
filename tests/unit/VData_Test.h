

#ifndef VDATA_TEST_H
#define VDATA_TEST_H

#include "Test.h"

namespace csmp {

template<size_t> class VSet;

// MISSING
// TODO: test read / write of VSet/VData with faces and interfaces
// TODO: test extraction of manifolds

/// PL Nov 2010 & SKM 23/12022
class VData_Test : public Test {
  public:
    virtual void run();
    
    // checks whether mesh is still intact after corner elements were split
    // test model FracBox
    bool TestReplacementOfCornerTetrahedra();
};


// builds fracbox model with 1872 elements and several tetrahedra spanning the corners
void create_FracBoxModel( VSet<3U>& );
    


} // csmp

#endif // VDATA_TEST_H
