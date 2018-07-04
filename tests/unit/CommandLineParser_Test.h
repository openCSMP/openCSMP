#ifndef CommandLineParser_Test_h__
#define CommandLineParser_Test_h__

#include <cstdlib>
#include <sstream>

#include "Test.h"
#include "CommandLineParser.h"

namespace csmp {

class CommandLineParser_Test : public Test {
public:
    CommandLineParser_Test();
    
    void run();
private:
    // Basic tests
    void intArgumentTest();
    void doubleArgumentTest();

	#ifndef _MSC_VER
		void stringArgumentTest();
	#endif

    void boolArgumentTest();
    void boolWithArgumentTest();
    void boolWithEndMarkerTest();
    void unusedDefaultTest();

    // Multiple argument test
    void sameArgumentTest();
    void longOptTest();
    void combinedTest();
    void sequenceTest();
    
    // Error tests
    void missingArgumentTest();
    void missingArgumentDefaultTest();
    void combinedMissingArgumentTest();
    void combinedMissingArgumentDefaultTest();
    void wrongOptionTest();
    void combinedWrongOptionTest();
    void endOfArgumentMarkerOptionTest();
    
    // Other functionality
    void nonSeparatedTest();
    void combinedNonSeparatedTest();
    void nonArgumentTest();
    void endOfOptionTest();
    void noNonArgumentsTest();
    void notUsedArgumentsTest();
    void outputTest();

};

} // end csp

#endif
