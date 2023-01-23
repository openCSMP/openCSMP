//
//  Experimental_Example.cpp
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//  Copyright (c) 2014 Stephan Matthai. All rights reserved.
//

#include "Experimental_Example.h"
#include "compareFloats.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#endif
#ifdef CSMP_WITH_MESCHACH
#include "Gauss_Solver.h"
#endif

#include <tuple>
//#include "BoundaryInterface_Test.h"
#include "Boundary.h"
#include "Region.h"
#include "ErrorHandler.h"
#include "ANSYS_Model3D.h"
#include "VTK_Interface.h"
#include "VTU_Interface.h"
#include "FaceConstructionData.h"


using namespace std;

namespace csmp {

void Experimental_Example::Specifications()
  {
     SetTitle( "Experimental_Example" );
     SetDifficulty( 1 );
     SetCategory( "Software Functionality" );
     AddAuthor( "You!" );
     AddDescription( "source in: Experimental_Example.cpp" );
     AddDescription( "Empty example for the user to experiment with" );
     AddRequirement( "none" );
     AddRequirement( "no predefined model or variables file" );
  }





/**
     Put CSMP code that you would like to test here and run it as part of the 
     example suite.

*/
void Experimental_Example::Run()
  {
    // Making a material matrix for illustrating the FE integral accumulation pattern for vector variable
    const uint32_t dim{ 2U }; // spatial dimension of model
    const uint32_t nodes{ 3U }; // a triangle element
    const uint32_t dof{ dim * nodes }; // degrees of freedom associated with element
    const double   value_label_for_entry{ 1. };
    
    DenseMatrix<DM_MIN>  LHS( dof, dof );
    LHS = 0.;
        
    // for a VECTOR node variable the components are expanded out to achieve an ordering comp1, comp2, comp3... in the diagonal and RH vector
    // tagging the different components with integers, meaning that if the scalar is 3, the component becomes 13 and so forth
    {
        for ( uint32_t i{0U}; i < nodes; ++i )
          {
             for ( uint32_t j{0U}; j < nodes; ++j )
               {
                  // for positive diagonal elements
                  if ( i == j ) {
                      for ( uint32_t k{0U}; k<dim; ++k )
                        LHS(i*dim+k,i*dim+k) = value_label_for_entry + k + 10;
                    }
                  // for negative off-diagonal elements
                  else {
                       for ( uint32_t k{0U}; k<dim; ++k )
                         LHS(i*dim+k, j*dim+k) = -(value_label_for_entry + k + 10);
                    }
               }
          }
        LHS.Out(0);
    }
    
} // end






} // csmp
