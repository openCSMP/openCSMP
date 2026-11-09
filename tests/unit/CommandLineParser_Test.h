// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COMMAND_LINE_PARSER_TEST_H
#define COMMAND_LINE_PARSER_TEST_H

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
