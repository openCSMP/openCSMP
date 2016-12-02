#export LM_LICENSE_FILE="/usr/local/home/junchulk/samg/license.dat"
#export LD_LIBRARY_PATH="/usr/local/home/junchulk/csmp/csmp-api-library/bin/linux/:/opt/intel/lib/intel64/"

if [ -f "CSMP_number_type.h" ]
  then
    icc -c ../../source/CSMP_number_types_main.cpp
    icc -o CSMP_number_types_main CSMP_number_types_main.o
    ./CSMP_number_types_main
    cp ./CSMP_number_types.h ../../source/includes/model/
fi

bjam -j32 --toolset=intel --build-type=complete --layout=system optimization=speed link=shared threading=multi cxxflags="-std=c++11"

cp ../../bin/linux/csmptest ../../bin/data/
cd ../../bin/data
./csmptest