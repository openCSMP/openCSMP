#export LM_LICENSE_FILE="/home/junchulk/samg/license.dat"
#export LD_LIBRARY_PATH="/home/junchulk/csmp/csmp-api-library/bin/linux/"
export LD_LIBRARY_PATH="../../bin/linux/"

if [ -f "CSMP_number_type.h" ]
  then
    g++ -c ../../source/CSMP_number_types_main.cpp
    g++ -o CSMP_number_types_main CSMP_number_types_main.o
    ./CSMP_number_types_main
    cp ./CSMP_number_types.h ../../source/includes/model/
fi

bjam -j2 --toolset=gcc --build-type=complete --layout=system optimization=speed link=shared threading=multi cxxflags="-std=c++11"

cp ../../bin/linux/csmptest ../../bin/data/
cd ../../bin/data
./csmptest