// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef INDEX_TEST_H
#define INDEX_TEST_H

#include "Test.h"

namespace csmp
{
	struct Index;
    class Index_Test: public Test
    {
    public:
        Index_Test();
        ~Index_Test();
        void IndexOperatorsTest(Index&);
        void run();
    private:
    };
}
#endif // INDEX_TEST_H
