#ifndef DenseMatrix_Test_H
#define DenseMatrix_Test_H

#include "Test.h"
#include "DenseMatrix.h"
#include "vector"
#include "TensorVariable.h"
#include "Point.h"
#include <iomanip>

namespace csmp
{
    class DenseMatrix_Test : public Test
    {
      public:
        explicit DenseMatrix_Test( bool verbose=false );
        ~DenseMatrix_Test();
        virtual void run();
        DenseMatrix<4> A;
        DenseMatrix<4> ResA;
        DenseMatrix<4> B;
        DenseMatrix<4> ResAB;
        DenseMatrix<4> C;
        DenseMatrix<4> ResC;
        DenseMatrix<4> D;
        DenseMatrix<4> ResD;
        DenseMatrix<4> E;
        DenseMatrix<4> F;
        DenseMatrix<4> ResF;
        DenseMatrix<4> G;
        DenseMatrix<4> ResG;
        DenseMatrix<4> H;

        TensorVariable<1> TensorVariable1U;
        TensorVariable<2> TensorVariable2U;
        TensorVariable<3> TensorVariable3U;
        VectorVariable<1> VectorVariable1U;
        VectorVariable<2> VectorVariable2U;
        VectorVariable<3> VectorVariable3U;

        DenseMatrix<3> DenseMatrix3x3;

        Point<1> Point1U;
        Point<2> Point2U;
        Point<3> Point3U;

        std::vector<double> x, y, sol_y, solB, solD;
      
        const bool verbose_;
    };
}

#endif // DenseMatrix_Test_H

