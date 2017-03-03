#include "Visitor_TestSuite.h"
#include "CopyReplaceVisitor_Test.h"
#include "ANSYS_Model3D.h"
#include "InputDataManager.h"

namespace csmp{

void Visitor_TestSuite::run()
{
  
  // Establishing test Model
  model_ = new ANSYS_Model3D( "FracBox", "CSMP-2phase-variables.txt", true );
                                            
  // Adding Visitor Tests
  suite_.addTest( new CopyReplaceVisitor_Test( model_ ) );
  //suite_.addTest( new ConstraintPointToNodePropertyVisitor_Test( *model_ ) ); DEACTIVATED UNTIL VISITOR IS FIXED
  
} // run

Visitor_TestSuite::~Visitor_TestSuite()
{
    if( model_ != 0 )
        delete model_;
}

} // csmp
