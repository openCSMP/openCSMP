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
#include "VTK_Interface.h"
#include "vset_makers.h"
//#include "ModelComparator.h"
#include "Element.h"
#include "Box.h"

using namespace std;

namespace csmp {


/**
     Questions
     - are the flags correctly assigned
     - establish clear responsibility sharing between MeshManager and AllElements region !

     Observations
     - AllElements in single domain model is non-unique, but should be unique
 
*/
void ModelSubDomain_Test::run()
  {
     VSet<3U>   vset;
     const bool skewed(false), isoparametric(true), verbose(true);
     test_Create_Prism_Hexa_VSet( vset, skewed );
     Model<3U>          model1( vset, isoparametric );
    
     //model1.CreateProperty( "box flag", "none", SCALAR, NODE );
     //model1.CreateProperty( "box flag element", "none", SCALAR, ELEMENT );
     //boxFlagsToVariable( model1, "box flag", "box flag element" );
     //VTK_Interface<3U>  vtk_output;
     //vtk_output.OutputDataToVTK( model1, "node_flag", "box flag", 0 );
     //vtk_output.OutputDataToVTK( model1, "element_flag", "box flag element", 0 );
    
     cerr <<"\nModelSubDomain_Test::run: original model.";
     // TODO: numbering of boundary nodes does not seem to be correct
     vset.Out();
     //if ( verbose ) model1.Out();
     model1.Region("All Elements").Out();
    
     model1.OutputToBinaryFile("ModelSubDomain_Test");
    
     // testing: model2.InputFromBinaryFile("ModelSubDomain_Test");
     Model<3U>  model2( string("ModelSubDomain_Test") );
     //vtk_output.OutputDataToVTK( model2, "node_flag", "box flag", 1 );
     //vtk_output.OutputDataToVTK( model2, "element_flag", "box flag element", 1 );
    
     cerr <<"\nModelSubDomain_Test::run: model reconstructed from disk.";
     //if ( verbose ) model2.Out();
    
     _test( CompareModelSubdomains( model1.Region("All Elements"), model2.Region("All Elements"), verbose ) );
     _test( CompareModelSubdomains( model1.Region("Model"), model2.Region("Model"), verbose ) );
  }



//template bool compareModelSubdomains( const ModelSubDomain<3,Element>&, const ModelSubDomain<3,Element>&, false );



} // end csmp
