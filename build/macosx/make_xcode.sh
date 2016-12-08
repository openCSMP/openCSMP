export LM_LICENSE_FILE="/Users/gocdslave/samg/license.dat"
export LD_LIBRARY_PATH="/Users/gocdslave/csmp/csmp-api-library/bin/macosx/"

if [ -f "CSMP_number_type.h" ]
  then
    xcrun g++ -c ../../source/CSMP_number_types_main.cpp
    xcrun g++ -o CSMP_number_types_main CSMP_number_types_main.o
    ./CSMP_number_types_main
    cp ./CSMP_number_types.h ../../source/includes/model/
fi

bjam -j32 --toolset=darwin --build-type=complete --layout=system optimization=speed link=static threading=multi cxxflags="-std=c++11"

cp ../../bin/macosx/csmptest ../../bin/data/
cd ../../bin/data
./csmptest