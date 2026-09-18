// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  DFM_PermeabilityAndFluxRatioAnalyzer_Test.cpp
//  CSMP_DFM_Analyzer
//
//

#include "DFM_PermeabilityAndFluxRatioAnalyzer_Test.h"
#include "TestSuite.h"

using namespace std;

namespace csmp {

DFM_PermeabilityAndFluxRatioAnalyzer_Test::DFM_PermeabilityAndFluxRatioAnalyzer_Test()
 : DFM_analyzer_( "aperture_20_frac", true /* with VTK output */ )
 {
 }


/** Unit Testing
   
    1. correct identification of model boundaries 
    
    2. integrity of test model: no mesh points out etc.
    
    3. conservative and accurate measurement of flux through the model
    
    4. correct determination of flow volumes
*/
void DFM_PermeabilityAndFluxRatioAnalyzer_Test::run()
{
   /// checking whether all perimeter nodes of the region "Model" are identified
   _test( DFM_analyzer_.TestBoundaryIntegrity() );

   /// checking the BOX boundary flagging of the test model
   /*
      - the interior lower-dim elements must not be flagged as boundary
   */
   DFM_analyzer_.OutputBoxBoundaryFlagsAsNumbers("fluid pressure");

   // checking that correctness of the flux calculation using the model 'box_hrz_fault'
   // (for the correct results see the Maple spreadsheet."
   const bool with_vtk_output(true);
   DFM_analyzer_.DiagonalTensorAnalysis( "box_hrz_fault", with_vtk_output );

} // Run()





} // end csmp
