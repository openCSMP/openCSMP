#ifndef PARAMETER_TEST_H
#define PARAMETER_TEST_H

#include "Test.h"

namespace csmp
{
    class Parameter;
    class Parameter_Test :public Test
    {
    public:
        explicit Parameter_Test();
        ~Parameter_Test();
        void IsWithinRangeTest();
        void RangeTest();
        void run();
    private:
        Parameter* t_param;

    };
};

#endif // PARAMETER_TEST_H
