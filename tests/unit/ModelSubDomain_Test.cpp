//
//  ModelSubDomain_Test.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "ModelSubDomain_Test.hpp"
#include "Model.h"
#include "CSMP_highLevelUtilities.h"
#include "vset_makers.h"

using namespace std;

namespace csmp {


void ModelSubDomain_Test::run()
  {
     VSet<3U>   vset;
     const bool skewed(false), isoparametric(true);
     test_Create_Prism_Hexa_VSet( vset, skewed );
     Model<3U>  model( vset, isoparametric );
    
     cerr <<"\nModelSubDomain_Test::run: original model.";
     model.Out();
    
     model.OutputToBinaryFile("ModelSubDomain_Test");
    
     Model<3U>  model2( string("ModelSubDomain_Test") );
     //model2.InputFromBinaryFile("ModelSubDomain_Test");
    
     cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
     model2.Out();
  }

} // end csmp
