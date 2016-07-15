#ifndef CUBICSPLINE_TEST_H
#define CUBICSPLINE_TEST_H

#include "Test.h"
#include "CSMP_number_types.h"
#include "vector"
#include "CubicSpline.h"

namespace csmp
{
    class CubicSpline_Test : public Test
    {
      public:
        CubicSpline_Test();
		~CubicSpline_Test();
        virtual void run();

      private:
        void InitializeTest( const char* );
        void InitializeTest( const std::vector<double64>&,
                             const std::vector<double64>&,
                             const double64,
                             const double64 );
        void ValueTest();
        void DerivativeTest();
        void MaxDerivativeTest();
        void Range_xTest();
        void Range_fxTest();
        const double64 tolerance_;
        void OutTest();

        CubicSpline cspline1;
        CubicSpline cspline2;
        CubicSpline cspline1_copy;
        CubicSpline* cspline2_copy;
        string datafile;
        std::vector<double64> x;
        std::vector<double64> fx;
        std::vector<double64> dfx;
        double64 df1;
        double64 dfn;

    };
}
#endif // CUBICSPLINE_TEST_H
