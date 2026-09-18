Unit Test Framework
===================
Author: Adrian Burri (burriad@student.ethz.ch)
Date: 20.10.2003


Motivation
----------
Only a small fraction of a developer's work is actually coding an application. 
More time is spent debugging, which is often an unproductive and tedious work. 
Unit testing offers a way to alleviate the problem by providing a means to catch 
many bugs while programming.

The core idea is to write a collection of tests which test the behaviour of a 
class (or a cluster of classes). In the CSP framework, these tests are written 
in a C++ framework that makes it easy to add new tests and offers means to 
automate them (i.e. the computer can test the correctness of the tests' execution 
- no manual checking required) so that they can be run after every change in the code.

But testing is more than a way to reduce bugs in code. Unit tests are a key ingredient
 to a programming paradigm called Extreme Programming (XP) [1].In XP, tests simplify
 changing written code  because the developer sees whether the assumptions of the code
 - expressed as tests - still hold or if the behaviour is altered. This makes it easier
 to refactor code [2]. Moreover, tests assist the developer to understand the code, 
both when writing and when changing it.


Framework
---------
The core of the framework is taken from an article in the C/C++ User Journal [3]. 
The classes explained therein are stored in the CVS repository in the 'unit_test' folder.

Test:
The basic idea is that personal test classes should be derived from the framework's Test class. 
In the personal test class, one method for each test is created an registered as a test function 
in the run() method defined in Test. The test class creates and manages the classes to be tested 
(called 'fixture' in most testing frameworks)

Methods of class Test:
_test(bool condition): Evaluates a boolean expression. If the condition is true, the test passes, 
otherwise it fails. If a test fails, a detailed summary of the failure is created.
_equal(double expr, double value, double tol): Checks an equality condition and allows for a user 
defined tolerance.
_fail(string comment): Can be regarded as a condition that always fails.
_succeed(): Can be regarded as a condition that is always true.

run(): Runs the test cases. All test methods must be registered here.
(Remark: _test, _equal and _fail are macros, that's why an underscore is used)

Suite:
Test classes are gathered together in test suites. A test suite executes its tests and creates 
a report about a test run.

Methods of class Suite:
setStream(ostream* stream): Sets output stream (cout will do for most purposes).
addTest(Test* test): Adds a test class to the suite.
addSuite(Suite* suite): Adds a whole test suite to the suite. Allows for test hierarchies.
run(): Runs all tests in the suite. Calls the run() method of all suites and test classes it contains.

ToDo: Extending the framework to do functional tests

Guidelines
----------
The following are guidelines about how to use the testing framework as well as naming conventions.

Tests:
* Test classes are derived from the Test class in the framework or subclasses thereof
* For each actual class in the CSP framework, a related test class is created
* The name of the class is <CSP_class>_Test
* The test class is stored in the unit_test folder

Test suites:
* For each library (folder), one test suite collecting all test cases from the actual folder is created
* For each library (folder) a target which executes the test suite is created

Additional guidelines about what and when should be tested is indicated next:

Usage:
* Tests can be used by all developers to verify the correctness of a class or a collection of classes. 
Therefore, the test classes are put under CVS as well (in the same folder as the class it is testing)
* Tests should cover basic functionality of a class
* Tests should include simple relations between classes (e.g. null references, if allowed)
* Tests should capture normal AND exceptional behaviour (i.e. tests that fail to indicate wrong usage 
and limitations of the code)
* When catching a bug, write a test that exposes the bug before correcting it. (This test can be 
subsequently used to make sure that one works with a version of the code in which the bug is actually fixed)
* Tests are kept and run regularly even if the component is regarded as stable
* The tests are run after each change in the source code (i.e. after compiling; this could 
be automated as well)
* All tests are run at least weekly
* Additional tools used to create the test data (i.e. scripts, tables etc.) are mentioned 
in the code of the test classes and are stored under CVS as well

Example
-------
main_test_suite.cpp contains an example of a test program. It only tests the CommandLineParser 
class now, but more could (and should) be added to this file.

References
----------
[1] Extreme Programming Explained, Embrace Change; Beck, K.; Addison-Wesley
[2] Refactoring, Improving the Design of Existing Code; Fowler, M.; Addison-Wesley
[3] The Simplest Automated Unit Test Framework That Could Possibly Work; Chuck 
	Allison; C/C++ User Journal (http://www.cuj.com/documents/s=8035/cuj0009allison1/)
