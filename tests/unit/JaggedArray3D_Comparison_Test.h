//
//  JaggedArray3D_Comparison_Test.hpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 25/11/2024.
//

#ifndef CSMP_JAGGED_ARRAY_3D_COMPARISON_TEST_H
#define CSMP_JAGGED_ARRAY_3D_COMPARISON_TEST_H

#include "Test.h"

namespace csmp {

class JaggedArray3D_Comparison_Test : public Test {
   public:
     virtual void run();
     
     void ComparePerformanceWithVectorOfVectors();
     
     static const bool verbose_ = true;
};

} // end csmp

#endif /* CSMP_JAGGED_ARRAY_3D_COMPARISON_TEST_H */
