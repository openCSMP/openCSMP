
#ifndef SparseMatrix_Test_H
#define SparseMatrix_Test_H

#include "Test.h"
#include "CSMP_number_types.h"
#include "SparseMatrix.h"
#include "vector"
#include "TensorVariable.h"
#include "Point.h"
#include "Matrix.h"
#include <iomanip>

namespace csmp
{
    class SparseMatrix_Test : public Test
    {
      public:
        SparseMatrix_Test();
        ~SparseMatrix_Test();
        void run();
        SparseMatrix A;
        SparseMatrix ResA;
        SparseMatrix B;
        SparseMatrix BB;
        SparseMatrix BBB;
        SparseMatrix ResAB;
        SparseMatrix C;
        SparseMatrix ResC;
        SparseMatrix D;
        SparseMatrix ResD;
        SparseMatrix E;
        SparseMatrix F;
        SparseMatrix ResF;
        SparseMatrix G;
        SparseMatrix ResG;
        SparseMatrix H;
        SparseMatrix SparseMatrix2x2;
        SparseMatrix SparseMatrix3x3;

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
        std::vector<size_t> sizetVector;

        double64 fTolerance;




    };
}

#endif // SparseMatrix_Test_H

