// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CREATE_UNIQUE_COMBINATIONS_TEST_H
#define CREATE_UNIQUE_COMBINATIONS_TEST_H

#include "Test.h"
#include "Point.h"

namespace csmp {

/**
@author SKM
@date 25/5/22

*/
	class createUniqueCombinations_Test : public Test {
	public:
		explicit createUniqueCombinations_Test();
		~createUniqueCombinations_Test();
    
		virtual void run();

	private:
		static const bool verbose_ = true;
	};

}
#endif /* CREATE_UNIQUE_COMBINATIONS_TEST_H */

