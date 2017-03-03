#ifndef POINT_TEST_H
#define POINT_TEST_H
/*
 *  Point_Test.h
 *  Created by Stephan Matthai on 10/26/10.
 *
 */
#include "Test.h"

namespace csmp {

class Point_Test : public Test {
  public:
    Point_Test();
    virtual ~Point_Test(){};
    
    virtual void run();
    
  private:
    void Test_1D_Point();
    void Test_2D_Point();
    void Test_3D_Point();
};

} // end 

#endif

