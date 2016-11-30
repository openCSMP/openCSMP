rem g++ -c ../../source/CSMP_number_types_main.cpp
rem g++ -o CSMP_number_types_main CSMP_number_types_main.o
rem ./CSMP_number_types_main
rem copy ./CSMP_number_types.h ../../source/includes/model/

bjam --toolset=msvc-14.0 --build-type=complete --layout=system optimization=speed address-model=64 link=static threading=multi 

set LM_LICENSE_FILE="../../bin/license.dat"
set LD_LIBRARY_PATH="../../bin/"

rem copy ../../bin/windows/csmptest.exe ../../bin/data/csmptest.exe
rem cd ../../bin/data
rem csmptest.exe
