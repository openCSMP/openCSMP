//
//  ModelSubDomain_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef ModelSubDomain_Test_hpp
#define ModelSubDomain_Test_hpp

#include "Test.h"
#include "ModelSubDomain.h"
#include "CSMP_mathUtilities.h"
#include "Box.h"

namespace csmp {

class ModelSubDomain_Test : public Test {
  public:
  
    // TODO: does not test any of the public interface methods of ModelSubDomain
    virtual void run();

    /// tests method that creates neighbor connectivity inside of CSMP
    bool Test_EstablishNeighborConnectivity();

    /// compares node locations and connectivity
    template<uint32_t dim,template<uint32_t> class simplicial_complex>
    bool CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>&,
                                 const ModelSubDomain<dim,simplicial_complex>&,
                                 bool verbose );
};

} // end csmp

#endif /* ModelSubDomain_Test_hpp */
