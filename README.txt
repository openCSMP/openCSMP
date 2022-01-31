README.txt - Welcome to Open Source version of the Complex Systems Modelling Platform (CSMP++)

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

Compilation for different platforms;

Apple:
======

 -DPYRAMID_TRIANGULAR_FACETS
-march=core2 \
-std=c++14
 -stdlib=libc++ -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-float-equal

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
