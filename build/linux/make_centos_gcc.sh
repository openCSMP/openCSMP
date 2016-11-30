g++ -c ../../source/CSMP_number_types_main.cpp
g++ -o CSMP_number_types_main CSMP_number_types_main.o
./CSMP_number_types_main
cp ./CSMP_number_types.h ../../source/includes/model/

bjam -j32 --toolset=clang --build-type=complete --layout=system optimization=speed link=shared threading=multi cxxflags="-std=c++11"

env LM_LICENSE_FILE="../../bin/license.dat"
env LD_LIBRARY_PATH="../../bin/data"

cp ../../bin/linux/csmptest ../../bin/data/
cd ../../bin/data
./csmptest
