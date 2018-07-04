#ifndef Matrix_Test_H
#define Matrix_Test_H

#include "Test.h"
#include "Matrix.h"
#include "TensorVariable.h"
#include "Point.h"

namespace csmp {

class Matrix_Test : public Test {
  public:
    explicit Matrix_Test( bool verbose=false );
    ~Matrix_Test();
  
    virtual void run();

    Matrix A;
    Matrix ResA;
    Matrix B;
    Matrix ResAB;
    Matrix C;
    Matrix ResC;
    Matrix D;
    Matrix ResD;
    Matrix E;
    Matrix F;
    Matrix ResF;
    Matrix G;
    Matrix ResG;
    Matrix H;
    Matrix ResH;
    Matrix Matrix2x2;
    Matrix Matrix3x3;

    TensorVariable<1> TensorVariable1U;
    TensorVariable<2> TensorVariable2U;
    TensorVariable<3> TensorVariable3U;
    VectorVariable<1> VectorVariable1U;
    VectorVariable<2> VectorVariable2U;
    VectorVariable<3> VectorVariable3U;

    Point<1> Point1U;
    Point<2> Point2U;
    Point<3> Point3U;

    std::vector<double64> x, y, sol_y, solB, solD;
    const bool verbose_;
};

}

#endif // Matrix_Test_H

