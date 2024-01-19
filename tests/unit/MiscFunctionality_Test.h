#ifndef MISC_FUNCTIONALITY_TEST_H
#define MISC_FUNCTIONALITY_TEST_H

#include "Test.h"
#include "Point.h"

namespace csmp {

/**
@author SKM
@date 25/5/22

*/
	class MiscFunctionality_Test : public Test {
	public:
		explicit MiscFunctionality_Test();
		~MiscFunctionality_Test();
    
		virtual void run();

	private:
		static const bool verbose_ = true;
	};

}
#endif /* MISC_FUNCTIONALITY_TEST_H */

