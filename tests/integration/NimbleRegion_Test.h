#ifndef NIMBLE_REGION_TEST_H
#define NIMBLE_REGION_TEST_H

#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

/**
     Models the construction and modification of a NimbleRegion object.
 
     @author SKM
     @date 25/7/2019
*/
class NimbleRegion_Test : public Test	{
	public:
    NimbleRegion_Test( const std::string& model="tutorial1_input",
                       const std::string& variables_file="NimbleRegion_Test-variables.txt" );
                       
		virtual void run();
  
  private:
    Model<2U>*  model2D_;
};

} // csmp

#endif // NIMBLE_REGION_TEST_H
