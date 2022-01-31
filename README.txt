README.txt - Welcome to Open Source version of the Complex Systems Modelling Platform (CSMP++)
(last updated SKM 31/1/2022)

Organisation of this repository:
- source/include        header files from the API library
- source/modules        source files
- tests/unit            unit tests illutrating the functionality of classes and functions
- tests/integration     tests demonstrating the functionality of classes working together
- examples              START HERE - then use ExperimentalExamnple.cpp to write your own first code therein
- doc                   documentation : CSMP User Guide, FEM_conventions etc.

Remarks:
- In most cases there is a compilation unit (*.h and *.cpp files) for each class 

Issues:
- have a look at the trackers in Trello: https://trello.com/invite/b/Tcpg8IpD/1fb97a2da7c35fa3349f388915b20098/csmp-development
- https://trello.com/invite/b/WvR3rTvV/8bd188de9876be7443a7751053b54750/csmp-current-issues

Please start with the source code of the examples which required data from the respective /data
directory.

Compile the libraries in versions for debugging and running the code. 
For the former define DEBUG and NDEBUG for the release version. 
Note that the debug version contains many assert
statements that will safe-guard you when developing new code, highlighting when assumptions 
made in the design of CSMP are violated.
These are not contained in the release version, which merely error handles the contained 
and throws csmp::Exception and standard exception objects.

Compilation with the commercial Algebraic Multigrid Solver for Systems (SAMG), Fraunhofer Gesellschaft, Germany,
Provide the compiler with the flags

-DCSMP_WITH_SAMG_SOLVER -DSAMG_MULTIPLE_INSTANCES -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE

Compilation on different platforms:

Apple:
======

 -DPYRAMID_TRIANGULAR_FACETS
-march=core2 \
-std=c++14
 -stdlib=libc++ -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-float-equal

 XCode - see csmp-api-library/development/XCode-creating-CSMP-examples-project.pdf

for debug add: -O0
for release add: -O3

Linux:
======
-DCSMP_WITH_SAMG_SOLVER
 -DSAMG_MULTIPLE_INSTANCES -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS


Clang:
-std=c++14
 -stdlib=libc++ -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-float-equal  -march=corei7-avx
GCC:
-Winline
 -Wall -fomit-frame-pointer -march=sandybridge -std=c++14 

ICC:
-fomit-frame-pointer
 -xCORE-AVX-I -std=c++14 -O3 -ipo-separate

for debug add: -O0
for release add: -O3


Windows: (MS Visual Studio, Community Edition)
==============================================

-DCSMP_WITH_SAMG_SOLVER
 -DSAMG_MULTIPLE_INSTANCES -DPYRAMID_TRIANGULAR_FACETS -D_CRT_SECURE_NO_WARNINGS \
/bigobj
 /EHsc

for debug add: /O0 -DEBUG
for release add: /O2 -NDEBUG
 /Oi /Ot
