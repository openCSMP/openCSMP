// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_PARAMETER_TEST_H
#define CSMP_PARAMETER_TEST_H

#include "Test.h"

namespace csmp
{
    struct Parameter;
    
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

#endif // CSMP_PARAMETER_TEST_H
