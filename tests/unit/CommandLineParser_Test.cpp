#include "CommandLineParser_Test.h"

#include <cstring>

using namespace std;

namespace csmp {

CommandLineParser_Test::CommandLineParser_Test() {}

void CommandLineParser_Test::run() {
    // Basic tests
    intArgumentTest();
    doubleArgumentTest();
#ifndef _MSC_VER
    stringArgumentTest();
#endif
    boolArgumentTest();
    boolWithArgumentTest();
    boolWithEndMarkerTest();
    unusedDefaultTest();

    // Multiple argument test
    sameArgumentTest();
    longOptTest();
    combinedTest();
    sequenceTest();
    
    // Error tests
    missingArgumentTest();
    missingArgumentDefaultTest();
    combinedMissingArgumentTest();
    combinedMissingArgumentDefaultTest();
    wrongOptionTest();
    combinedWrongOptionTest();
    endOfArgumentMarkerOptionTest();
    
    // Other functionality
    nonSeparatedTest();
    combinedNonSeparatedTest();
    nonArgumentTest();
    endOfOptionTest();
    noNonArgumentsTest();
    notUsedArgumentsTest();
    outputTest();
}

void CommandLineParser_Test::intArgumentTest()
{
  const char **argv = new const char*[3];
  argv[0] = "program";
  argv[1] = "-i";
  argv[2] = "10";

  CommandLineParser parser(3, argv);
  int i = parser.get<int>("-i");
  int h = parser.get<int>("-h", 11);

  _test(i == 10);
  _test(h == 11);
 
  delete[] argv;
}

void CommandLineParser_Test::doubleArgumentTest()
{
  const char **argv = new const char*[3];
  argv[0] = "program";
  argv[1] = "-f";
  argv[2] = "2.0";
  
  CommandLineParser parser(3, argv);
  double f = parser.get<double>("-f");
  double g = parser.get<double>("-g", 3.0);
  
  _test(f == 2.0);
  _test(g == 3.0);
  
  delete[] argv;
}

#ifndef _MSC_VER
void CommandLineParser_Test::stringArgumentTest()
{
  const char **argv = new const char*[5];
  argv[0] = "program";
  argv[1] = "-s";
  argv[2] = "Hello";
  argv[3] = "-c";
  argv[4] = "World";
  
  csmp::CommandLineParser parser(5, argv);
 
  string s = parser.get<const char*>("-s");
  const char* c = parser.get<const char*>("-c");
  string t = parser.get<const char*>("-t", "Gotcha");
   
  _test(!std::strcmp(c, "World"));
  _test(s == "Hello");
  _test(t == "Gotcha");
  
  delete[] argv;
}
#endif

void CommandLineParser_Test::boolArgumentTest()
{
  const char **argv = new const char*[3];
  argv[0] = "program";
  argv[1] = "-b";
  
  CommandLineParser parser(2, argv);
  bool f = parser.get<bool>("-b");
  bool g = parser.get<bool>("-g");
  
  _test(f);
  _test(!g);

  delete[] argv;
}

void CommandLineParser_Test::boolWithArgumentTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-b";
    argv[2] = "Hallo";
    
    bool f;
    try
	{
        CommandLineParser parser(3, argv);
        f = parser.get<bool>("-b");
    }
	catch (std::runtime_error& e)
	{
      _succeed();
	  delete[] argv;
      return;
    }
    _fail("boolWithArgumentTest");
    
	delete[] argv;
}

void CommandLineParser_Test::boolWithEndMarkerTest()
{
  const char **argv = new const char*[4];
  argv[0] = "program";
  argv[1] = "-b";
  argv[2] = "--";
  argv[3] = "Hallo";
  
  bool f;
  CommandLineParser parser(4, argv);
  try
  {
    f = parser.get<bool>("-b");
  }
  catch (std::runtime_error& e)
  {
    _fail("boolWithEndMarkerTest: exception");
  }
  _test(f);
  _test(parser.getOperands() == "Hallo");
  
  delete[] argv;
}

void CommandLineParser_Test::unusedDefaultTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0e-2";
    
    CommandLineParser parser(3, argv);
    double d = parser.get<double>("-d", 3.0);
    
    _test(d == 2.0e-2);
    
	delete[] argv;
}

void CommandLineParser_Test::sameArgumentTest()
{
    const char **argv = new const char*[5];
    argv[0] = "program";
    argv[1] = "-s";
    argv[2] = "Hello";
    argv[3] = "-s";
    argv[4] = "World";
    
    try 
	{
        CommandLineParser parser(5, argv);
        string str = parser.get<const char*>("-s");
    }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }

    _fail("sameArgumentTest: argument twice accepted");
    
	delete[] argv;
}

void CommandLineParser_Test::longOptTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "--double";
    argv[2] = "2.0";
    
    CommandLineParser parser(3, argv);
    double d = parser.get<double>("--double");
    
    _test(d == 2.0);
    
	delete[] argv;
}

void CommandLineParser_Test::combinedTest()
{
    const char **argv = new const char*[5];
    argv[0] = "program";
    argv[1] = "--double";
    argv[2] = "1.0e3";
    argv[3] = "-i";
    argv[4] = "10";
    
    CommandLineParser parser(5, argv);
    double d = parser.get<double>("-d", "--double");
    int i = parser.get<int>("-i", "--integer");
    
    _test(d == 1.0e3);
    _test(i == 10);
    
	delete[] argv;
}

void CommandLineParser_Test::sequenceTest()
{
    const char **argv = new const char*[5];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "1.0e3";
    argv[3] = "-i";
    argv[4] = "10";
    
    CommandLineParser parser(5, argv);
    parser.get<int>("-i");
    parser.get<double>("-d");
    
	delete[] argv;
}

void CommandLineParser_Test::missingArgumentTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    try	{
         CommandLineParser parser(3, argv);
         parser.get<double>("-e");
      }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }
    _fail("missingArgumentTest");
    
	delete[] argv;
}

void CommandLineParser_Test::missingArgumentDefaultTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("-e", 3.0);
    }
	catch (std::runtime_error& e)
	{
        _fail("missingArgumentDefaultTest");
		delete[] argv;
        return;
    }
    _test(d == 3.0);
    
	delete[] argv;
}

void CommandLineParser_Test::combinedMissingArgumentTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("-w", "--wrong");
    }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }
    _fail("combinedMissingArgumentTest()");
    
	delete[] argv;
}

void CommandLineParser_Test::combinedMissingArgumentDefaultTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("-w", "--wrong", 3.0);
    }
	catch (std::runtime_error& e)
	{
        _fail("combinedMissingArgumentDefaultTest");
		delete[] argv;
        return;
    }
    _test(d == 3.0);
    
	delete[] argv;
}


void CommandLineParser_Test::wrongOptionTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("d", 3.0);
    }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }
    _fail("wrongOptionTest");
    
	delete[] argv;
}

void CommandLineParser_Test::combinedWrongOptionTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("-d", "-wrong", 3.0);
    }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }

    _fail("combinedWrongOptionTest");
	
	delete[] argv;
}

void CommandLineParser_Test::endOfArgumentMarkerOptionTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";
    
    double d;
    try
	{
        CommandLineParser parser(3, argv);
        d = parser.get<double>("-d", "--", 3.0);
    }
	catch (std::runtime_error& e)
	{
        _succeed();
		delete[] argv;
        return;
    }

    _fail("endOfArgumentMarkerOptionTest");
    
	delete[] argv;
}

void CommandLineParser_Test::nonSeparatedTest()
{
    const char **argv = new const char*[2];
    argv[0] = "program";
    argv[1] = "-d2.0";
    
    double d;
    try
	{
        CommandLineParser parser(2, argv);
        d = parser.get<double>("-d");
    }
	catch (std::runtime_error& e)
	{
        _fail("nonSeparatedTest: exception");
    }

    _test(d == 2.0);
    
	delete[] argv;
}

void CommandLineParser_Test::combinedNonSeparatedTest()
{
    const char **argv = new const char*[2];
    argv[0] = "program";
    argv[1] = "-d2.0";
    
    double d;
    try
	{
        CommandLineParser parser(2, argv);
        d = parser.get<double>("-d", "--double");
    }
	catch (std::runtime_error& e)
	{
        _fail("nonSeparatedTest: exception");
        std::cout << e.what() << std::endl;
    }

    _test(d == 2.0);
    
	delete[] argv;
}

void CommandLineParser_Test::nonArgumentTest()
{
    const char **argv = new const char*[5];
    argv[0] = "program";
    argv[1] = "-g";
    argv[2] = "Neue";
    argv[3] = "Zuercher";
    argv[4] = "Zeitung";
    
    CommandLineParser parser(5, argv);
    string nonOpt = parser.getOperands();
    
    _test(parser.operands() == true);
    _test(nonOpt == "ZuercherZeitung");
    
	delete[] argv;
}

void CommandLineParser_Test::endOfOptionTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "--";
    argv[2] = "Hallo";
    
    CommandLineParser parser(3, argv);
    string nonOpt = parser.getOperands();
    
    _test(parser.operands() == true);
    _test(nonOpt == "Hallo");
    
	delete[] argv;
}

void CommandLineParser_Test::noNonArgumentsTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-s";
    argv[2] = "Hallo";
    
    CommandLineParser parser(3, argv);
    string nonOpt = parser.getOperands();
    
    _test(parser.operands() == false);
    _test(nonOpt.empty() == true);
    
	delete[] argv;
}

void CommandLineParser_Test::notUsedArgumentsTest()
{
    const char **argv = new const char*[3];
    argv[0] = "program";
    argv[1] = "-d";
    argv[2] = "2.0";

    CommandLineParser parser(3, argv);
    _test(parser.unused());
    _test(parser.getUnused() == "-d 2.0");
    
    parser.get<double>("-d");
            
    _test(!parser.unused());
    _test(parser.getUnused().empty());
    
	delete[] argv;
}

void CommandLineParser_Test::outputTest()
{
    const char **argv = new const char*[5];
    argv[0] = "program";
    argv[1] = "-g";
    argv[2] = "Neue";
    argv[3] = "Zuercher";
    argv[4] = "Zeitung";
    
    CommandLineParser parser(5, argv);
    std::ostringstream os;
    os << parser;

    _test(os.str() == "program -g Neue Zuercher Zeitung");
    
	delete[] argv;
}

} // end csmp

