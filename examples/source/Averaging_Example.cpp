#include "Averaging_Example.h"

#include "ANSYS_Model3D.h"
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
  ANSYS_Model3D model( "LeftRight", "CSMP-1phase-variables.txt", true );
  PropertyHandle<3U> sourceSink( model, "fluid volume source", SCALAR, ELEMENT );
  PropertyHandle<3U> nodalSourceSinkByDistance( model, "nodal fluid volume source distance", SCALAR, NODE );
  PropertyHandle<3U> nodalSourceSinkByVolume( model, "nodal fluid volume source volume", SCALAR, NODE );
  PropertyHandle<3U> nodalSourceSinkByCount( model, "nodal fluid volume source count", SCALAR, NODE );

  // resetting target variables
  nodalSourceSinkByDistance = 0.;
  nodalSourceSinkByVolume = 0.;
  nodalSourceSinkByCount = 0.;

  // setting source variables
  sourceSink = 1;

  // setting up VTU interface for visualization
  VTU_Interface<3U> vtu( model, "Extrapolation Comparison" );

  // creating a reference to the model region
  Region<3>& region( model.Region( "Model" ) );

  // using the ModelSubDomain methods
  region.ExtrapolateElementToNodeProperty( "fluid volume source", "nodal fluid volume source distance" );
  region.ExtrapolateElementToNodeProperty( "fluid volume source", "nodal fluid volume source volume", false );

  // using shahos node count function
  extrapolateElementToNodalVariable( model, "Model", "fluid volume source", "nodal fluid volume source count" );

  // output
  vtu.OutputDataToVTU( "ExtrapolationByDistance", "nodal fluid volume source distance", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByVolume", "nodal fluid volume source volume", "Model", static_cast<int>(0) );
  vtu.OutputDataToVTU( "ExtrapolationByCount", "nodal fluid volume source count", "Model", static_cast<int>(0) );

}

} // csmp
