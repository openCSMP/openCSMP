#include "Averaging_Example.h"

#include "ANSYS_Model2D.h"
#include "Region.h"
#include "VTU_Interface.h"
#include "PL_Utilities.h"

namespace csmp {

void Averaging_Example::Specifications()
  {
    SetTitle( "Extrapolation & Averaging" );
    SetDifficulty( 2 );
    SetCategory( "Software Functionality" );
    AddAuthor( "P. Lang" );
    AddDescription( "Comparison of element to node extrapolation approaches" );
    AddDescription( "Illustration of VTU Interface" );
    AddDescription( "source in: Averaging_Example.cpp" );
    AddRequirement( "file set: 'LeftRight'");
    AddRequirement( "CSMP-1phase-variables.txt");
  }



/**
    Visualisation of element-to-node extrapolation options within csmp.
*/
void Averaging_Example::Run()
{
  // initializing model and properties
  ANSYS_Model2D model( "LeftRight", "CSMP-1phase-variables.txt", true );
  PropertyHandle<2U> sourceSink( model, "fluid volume source", SCALAR, ELEMENT );
  PropertyHandle<2U> nodalSourceSinkByDistance( model, "nodal fluid volume source distance", SCALAR, NODE );
  PropertyHandle<2U> nodalSourceSinkByVolume( model, "nodal fluid volume source volume", SCALAR, NODE );
  PropertyHandle<2U> nodalSourceSinkByCount( model, "nodal fluid volume source count", SCALAR, NODE );

  // resetting target variables
  nodalSourceSinkByDistance = 0.;
  nodalSourceSinkByVolume = 0.;
  nodalSourceSinkByCount = 0.;

  // setting source variables
  sourceSink = 1.;

  // setting up VTU interface for visualization
  VTU_Interface<2U> vtu( model, "Extrapolation Comparison" );

  // creating a reference to the model region
  Region<2U>& region( model.Region( "Model" ) );

  // using the ModelSubDomain methods to extrapolate single-valued field
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source distance" );
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source volume", false );

  // output
  vtu.OutputDataToVTU( "ExtrapolationByDistance", "nodal fluid volume source distance", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByVolume", "nodal fluid volume source volume", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByCount", "nodal fluid volume source count", "Model", static_cast<int>(0) );

  // creating a perturbed field and extrapolating this
  randomPerturb( model, "fluid volume source", 20. );

  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source distance" );
  region.ExtrapolateCellToNodeProperty( "fluid volume source", "nodal fluid volume source volume", false );

  vtu.OutputDataToVTU( "PerturbedFluidVolumeSource", "fluid volume source", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByDistance", "nodal fluid volume source distance", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByVolume", "nodal fluid volume source volume", "Model", static_cast<int>(1) );
  vtu.OutputDataToVTU( "ExtrapolationByCount", "nodal fluid volume source count", "Model", static_cast<int>(1) );

}

} // csmp
