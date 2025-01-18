#ifndef CUBICSPLINE_TEST_H
#define CUBICSPLINE_TEST_H

#include "Test.h"
#include "CubicSpline.h"

namespace csmp
{
    class CubicSpline_Test : public Test
    {
      public:
        explicit CubicSpline_Test( bool verbose=false );
		    ~CubicSpline_Test();
        
        virtual void run();

      private:
        void InitializeTest( const char* );
        void InitializeTest( const std::vector<double>&,
                             const std::vector<double>&,
                             const double,
                             const double );
        void ValueTest();
        void DerivativeTest();
        void MaxDerivativeTest();
        void Range_xTest();
        void Range_fxTest();
        const double tolerance_;
        void OutTest();

        CubicSpline cspline1;
        CubicSpline cspline2;
        CubicSpline cspline1_copy;
        CubicSpline* cspline2_copy;
        std::string datafile;
        std::vector<double> x;
        std::vector<double> fx;
        std::vector<double> dfx;
        double df1;
        double dfn;
      
        const bool verbose_;
    };
}
#endif // CUBICSPLINE_TEST_H
