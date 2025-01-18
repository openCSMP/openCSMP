#ifndef TWOPHASEMODEL_TESTS_H
#define TWOPHASEMODEL_TESTS_H

#include "Test.h"
#include "TestSuite.h"

#include "Model1D.h"
#include "TwoPhaseModel_Test.h"
#include "PropertyHandle.h"
#include "VectorVariable.h"
#include "TestSuite.h"
#include "BrooksCorey.h"
#include "Corey.h"
#include "FourarLenormand.h"
#include "VanGenuchten.h"
#include "TwoPhaseFileBased.h"
#include "Experimental2PhaseModel.h"
#include "ExperimentalRT.h"
#include "FractureMatrixUpscaled.h"
#include "BrooksCoreyWithHysteresis.h"
#include "LinearTwoPhaseModel.h"

namespace csmp{

/// PL Nov 2010
class TwoPhaseModel_TestSuite : public Test
{
public:

  explicit TwoPhaseModel_TestSuite( TestSuite& suite ) : suite_(suite), free_( true ){}
  virtual ~TwoPhaseModel_TestSuite() { free(); }
  virtual void run();
  void free();

private:

  void AssignSaturationValues( Model<1U>* );

  TestSuite&  suite_;
  bool        free_;

  Model<1U>*  rock_model_;
  Model<1U>*  fracture_rock_model_;

};

/// resets model properties before test
void TwoPhaseModel_TestSuite::free()
{
	if( !free_  )
	{
		if( rock_model_ != NULL )
			delete rock_model_;
		if( fracture_rock_model_ != NULL )
			delete fracture_rock_model_;
		free_ = true;
	}

}


/// assign saturation values
void TwoPhaseModel_TestSuite::AssignSaturationValues( Model<1U>* model )
{

    // Saturation values
    Region<1>&  sg(model->Region("Model"));
    static csmp::Index   satw_key = model->Database().StorageKey("saturation water");
    static csmp::Index   satn_key = model->Database().StorageKey("saturation oil");
    const double       sat_incr(1./model->Mesh().Nodes());
    ScalarVariable       saturation;
    // generating a range of saturation values for water and oil
    for ( vector<Node<1U>*>::iterator
          it=sg.NodesBegin(); it!=sg.NodesEnd(); it++ )
      {
         saturation() = 0. + sat_incr * (*it)->Idx();
         (*it)->Store( satw_key, saturation );
         saturation = 1. - saturation;
         (*it)->Store( satn_key, saturation );
      }

}

void TwoPhaseModel_TestSuite::run()
{

  free_ = false;

  const uint32_t  N_ELEMENTS(100);
  const double length(1.);

  // creating 1D TestModel with a matrix and fracture region

  fracture_rock_model_  = new Model1D<1U>( "FractureRockModel1D", "CSMP-2phase-variables.txt", length, N_ELEMENTS );

  vector<uint32_t>   elms;
  elms.reserve( N_ELEMENTS );
  for ( uint32_t i = 0; i < 40; ++i )
    elms.push_back(i);
  for ( uint32_t i = 60; i < N_ELEMENTS; ++i )
    elms.push_back(i);
  fracture_rock_model_->FormRegionFrom( "ROCK", elms );
  elms.erase( elms.begin(), elms.end() );
  for ( uint32_t i=40; i<60; i++ )
    elms.push_back(i);
  fracture_rock_model_->FormRegionFrom( "FRACTURE", elms );


  // creating 1D TestModel with just a matrix region

  rock_model_           = new Model1D<1U>( "RockModel1D", "CSMP-2phase-variables.txt", length, N_ELEMENTS );

  elms.erase( elms.begin(), elms.end() );
  for ( uint32_t i = 0; i<N_ELEMENTS; i++ )
    elms.push_back(i);
  rock_model_->FormRegionFrom( "ROCK", elms );


  // Assign properties

  // Adding properties that are not in the variables file
  if ( !fracture_rock_model_->Database().IsDefined( "fracture matrix interface area" ) )
      fracture_rock_model_->CreateProperty( "fracture matrix interface area",  "SI", SCALAR, ELEMENT );
  if ( !fracture_rock_model_->Database().IsDefined( "fracture porosity" ) )
      fracture_rock_model_->CreateProperty( "fracture porosity",  "SI", SCALAR, ELEMENT );
  if ( !fracture_rock_model_->Database().IsDefined( "fracture matrix flux ratio" ) )
      fracture_rock_model_->CreateProperty( "fracture matrix flux ratio",  "SI", SCALAR, ELEMENT );
  if ( !fracture_rock_model_->Database().IsDefined( "block radius" ) )
      fracture_rock_model_->CreateProperty( "block radius",  "SI", SCALAR, ELEMENT );
  if ( !fracture_rock_model_->Database().IsDefined( "fracture aperture" ) )
      fracture_rock_model_->CreateProperty( "fracture aperture",  "SI", SCALAR, ELEMENT );

  if ( !fracture_rock_model_->Database().IsDefined( "previous water saturation barycenter" ) )
      fracture_rock_model_->CreateProperty( "previous water saturation barycenter",  "SI", SCALAR, ELEMENT );
  fracture_rock_model_->InputPropertyValue( "previous water saturation barycenter", makeScalar(PLAIN,0.2) );
  if ( !fracture_rock_model_->Database().IsDefined( "inflection water saturation barycenter" ) )
      fracture_rock_model_->CreateProperty( "inflection water saturation barycenter",  "SI", SCALAR, ELEMENT );
  fracture_rock_model_->InputPropertyValue( "inflection water saturation barycenter", makeScalar(PLAIN,0.2) );

  if ( !fracture_rock_model_->Database().IsDefined( "initial saturation water" ) )
      fracture_rock_model_->CreateProperty( "initial saturation water",  "SI", SCALAR, NODE );
  if ( !fracture_rock_model_->Database().IsDefined( "maximum residual oil saturation" ) )
      fracture_rock_model_->CreateProperty( "maximum residual oil saturation",  "SI", SCALAR, NODE );
  fracture_rock_model_->InputPropertyValue( "maximum residual oil saturation", makeScalar(PLAIN,0.4) );

  // Matrix properties

  VectorVariable<1U>  Vd(DIRICH,  1.0e-10 );   // Prescribed velocity

  const double fluid_pressure ( 5.0e7 );
  const double visc_water     ( 1.0e-3 );
  const double visc_oil       ( 3.0e-3 );
  const double dens_water     ( 1000. );
  const double dens_oil       ( 800. );

  vector<double> porosity       ( 1, 0.5 );
  vector<double> permeability   ( 1, 1.0e-12 );
  vector<double> swr            ( 1, 0.2 );
  vector<double> snr            ( 1, 0.25 );

  vector<double> lambda         ( 1, 2.0 );     // Brooks Corey, Corey
  vector<double> entry_pressure ( 1, 2.0e3 );   // Brooks Corey, Corey

  vector<double> exp_w          ( 1, 2.0 );     // Corey exponent for wetting phase
  vector<double> exp_n          ( 1, 1.5 );     // Corey exponent for non-wetting phase
  vector<double> krwr           ( 1, 1.0 );     // Corey end-point wetting phase relative permeability
  vector<double> krnr           ( 1, 0.98 );    // Corey end-point non-wetting phase relative permeability

  vector<double> VG_n           ( 1, 3.0 );     // Van Genuchten rel.perm. parameter
  vector<double> VG_alpha       ( 1, 0.37 );    // Van Genuchten cap. press. parameter

  vector<double> fracture_aperture( 1, 1.0-3 ); // Fourar Lenormand

  // Fracture Matrix interaction for rock model
  rock_model_->InputPropertyValue( "fracture aperture",                     makeScalar(PLAIN,fracture_aperture[0]) );

  // Matrix and global properties for rock model
  rock_model_->InputPropertyValue( "fluid pressure",                        makeScalar(PLAIN,fluid_pressure) );
  rock_model_->InputPropertyValue( "velocity",                              Vd );
  rock_model_->InputPropertyValue( "viscosity water",                       makeScalar(PLAIN,visc_water) );
  rock_model_->InputPropertyValue( "viscosity oil",                         makeScalar(PLAIN,visc_oil) );
  rock_model_->InputPropertyValue( "density water",                         makeScalar(PLAIN,dens_water) );
  rock_model_->InputPropertyValue( "density oil",                           makeScalar(PLAIN,dens_oil) );
  rock_model_->InputPropertyValue( "permeability",                          makeScalar(PLAIN,permeability[0]) );
  rock_model_->InputPropertyValue( "porosity",                              makeScalar(PLAIN,porosity[0]) );
  rock_model_->InputPropertyValue( "residual saturation wetting phase",     makeScalar(PLAIN,swr[0]) );
  rock_model_->InputPropertyValue( "residual saturation non-wetting phase", makeScalar(PLAIN,snr[0]) );
  rock_model_->InputPropertyValue( "van genuchten parameter",               makeScalar(PLAIN,VG_n[0]) );
  rock_model_->InputPropertyValue( "van genuchten alpha",                   makeScalar(PLAIN,VG_alpha[0]) );
  rock_model_->InputPropertyValue( "brooks corey parameter",                makeScalar(PLAIN,lambda[0]) );
  rock_model_->InputPropertyValue( "entry pressure",                        makeScalar(PLAIN,entry_pressure[0]) );
  rock_model_->InputPropertyValue( "corey exponent water",                  makeScalar(PLAIN,exp_w[0]) );
  rock_model_->InputPropertyValue( "corey exponent oil",                    makeScalar(PLAIN,exp_n[0]) );
  rock_model_->InputPropertyValue( "relperm endpoint water",                makeScalar(PLAIN,krwr[0]) );
  rock_model_->InputPropertyValue( "relperm endpoint oil",                  makeScalar(PLAIN,krnr[0]) );

  AssignSaturationValues( rock_model_);

  // Fracture Matrix interaction for fracture-rock model
  fracture_rock_model_->InputPropertyValue( "fracture aperture",                     makeScalar(PLAIN,fracture_aperture[0]) );
  fracture_rock_model_->InputPropertyValue( "fracture matrix interface area",        makeScalar(PLAIN,2.0) );
  fracture_rock_model_->InputPropertyValue( "fracture porosity",                     makeScalar(PLAIN,0.1) );
  fracture_rock_model_->InputPropertyValue( "fracture matrix flux ratio",            makeScalar(PLAIN,5.0) );
  fracture_rock_model_->InputPropertyValue( "initial saturation water",              makeScalar(PLAIN,0.01) ); // as needed for transfer function
  fracture_rock_model_->InputPropertyValue( "volume flux",                           makeScalar(PLAIN,1.0e-9) ); // for scaling of transfer contribution
  fracture_rock_model_->InputPropertyValue( "block radius",                          makeScalar(PLAIN,0.7) );

  // Matrix and global properties for fracture-rock model
  fracture_rock_model_->InputPropertyValue( "fluid pressure",                        makeScalar(PLAIN,fluid_pressure) );
  fracture_rock_model_->InputPropertyValue( "velocity",                              Vd );
  fracture_rock_model_->InputPropertyValue( "viscosity water",                       makeScalar(PLAIN,visc_water) );
  fracture_rock_model_->InputPropertyValue( "viscosity oil",                         makeScalar(PLAIN,visc_oil) );
  fracture_rock_model_->InputPropertyValue( "density water",                         makeScalar(PLAIN,dens_water) );
  fracture_rock_model_->InputPropertyValue( "density oil",                           makeScalar(PLAIN,dens_oil) );
  fracture_rock_model_->InputPropertyValue( "permeability",                          makeScalar(PLAIN,permeability[0]) );
  fracture_rock_model_->InputPropertyValue( "porosity",                              makeScalar(PLAIN,porosity[0]) );
  fracture_rock_model_->InputPropertyValue( "residual saturation wetting phase",     makeScalar(PLAIN,swr[0]) );
  fracture_rock_model_->InputPropertyValue( "residual saturation non-wetting phase", makeScalar(PLAIN,snr[0]) );
  fracture_rock_model_->InputPropertyValue( "van genuchten parameter",               makeScalar(PLAIN,VG_n[0]) );
  fracture_rock_model_->InputPropertyValue( "van genuchten alpha",                   makeScalar(PLAIN,VG_alpha[0]) );
  fracture_rock_model_->InputPropertyValue( "brooks corey parameter",                makeScalar(PLAIN,lambda[0]) );
  fracture_rock_model_->InputPropertyValue( "entry pressure",                        makeScalar(PLAIN,entry_pressure[0]) );
  fracture_rock_model_->InputPropertyValue( "corey exponent water",                  makeScalar(PLAIN,exp_w[0]) );
  fracture_rock_model_->InputPropertyValue( "corey exponent oil",                    makeScalar(PLAIN,exp_n[0]) );
  fracture_rock_model_->InputPropertyValue( "relperm endpoint water",                makeScalar(PLAIN,krwr[0]) );
  fracture_rock_model_->InputPropertyValue( "relperm endpoint oil",                  makeScalar(PLAIN,krnr[0]) );

  // Fracture properties
  Region<1>& fractures = fracture_rock_model_->Region( "FRACTURE" );
  fractures.InputPropertyValue( "permeability", makeScalar(PLAIN,1.0E-9) );
  fractures.InputPropertyValue( "porosity", makeScalar(PLAIN,1.0) );
  fractures.InputPropertyValue( "residual saturation wetting phase", makeScalar(PLAIN,0.0) );
  fractures.InputPropertyValue( "residual saturation non-wetting phase", makeScalar(PLAIN,0.0) );
  fractures.InputPropertyValue( "entry pressure", makeScalar(PLAIN,0.0) );
  fractures.InputPropertyValue( "brooks corey parameter", makeScalar(PLAIN,0.0) );

  AssignSaturationValues( fracture_rock_model_);

  // TwoPhaseModel tests

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new BrooksCorey<1U>( rock_model_->Database(),
                                                "permeability", "viscosity oil", "viscosity water",
                                                "density oil", "density water",
                                                "brooks corey parameter",
                                                "entry pressure",
                                                NULL,
                                                "saturation water",
                                                "residual saturation non-wetting phase", "residual saturation wetting phase"),
                                          "csmp::TwoPhaseModel_BrooksCorey_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure") );

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new Corey<1U>( rock_model_->Database(),
                                                "permeability", "viscosity oil", "viscosity water",
                                                "density oil", "density water",
                                                "brooks corey parameter",
                                                "entry pressure",
                                                NULL,
                                                "corey exponent oil","corey exponent water",
                                                "relperm endpoint oil","relperm endpoint water",
                                                "saturation water",
                                                "residual saturation non-wetting phase", "residual saturation wetting phase"),
                                          "csmp::TwoPhaseModel_Corey_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new LinearTwoPhaseModel<1U>( rock_model_->Database(),
                                                 "permeability", "viscosity oil", "viscosity water",
                                                 "density oil", "density water",
                                                 "entry pressure",
                                                  NULL,
                                                 "saturation water",
                                                 "residual saturation non-wetting phase", "residual saturation wetting phase"),
                                          "csmp::TwoPhaseModel_LinearTwoPhaseModel_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new VanGenuchten<1U>( rock_model_->Database(),
                                                 "permeability", "viscosity oil", "viscosity water",
                                                 "density oil", "density water",
                                                 "van genuchten parameter",
                                                 "van genuchten alpha",
                                                 "saturation water",
                                                 "residual saturation non-wetting phase", "residual saturation wetting phase"),
                                          "csmp::TwoPhaseModel_VanGenuchten_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new FourarLenormand<1U>( rock_model_->Database(),
                                                 "fracture aperture",
                                                 "permeability", "viscosity oil", "viscosity water",
                                                 "density oil", "density water",
                                                 "saturation water",
                                                 "residual saturation non-wetting phase", "residual saturation wetting phase"),
                                          "csmp::TwoPhaseModel_FourarLenormand_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );

  suite_.addTest( new TwoPhaseModel_Test( rock_model_,
                                          new TwoPhaseFileBased<1U>( rock_model_->Database(), "TwoPhaseFileBasedIn.txt" ),
                                          "csmp::TwoPhaseModel_TwoPhaseFileBased_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );


  /*
  suite_.addTest( new TwoPhaseModel_Test( model_,
                                          new BrooksCoreyWithHysteresis<1U>( model_->Database(),
                                                          "brooks corey parameter",
                                                          "entry pressure" ),
                                          "csmp::TwoPhaseModel_BrooksCoreyWithHysteresis_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure" ) );
                                          
  suite_.addTest( new TwoPhaseModel_Test( fracture_rock_model_,
                                          new FractureMatrixUpscaled<1U>( model_->Database(),
                                                          "permeability", "viscosity oil", "viscosity water",
                                                           "density oil", "density water",
                                                           "brooks corey parameter", "entry pressure",
                                                           "fracture matrix interface area",
                                                           "fracture porosity", "porosity",
                                                           "fracture matrix flux ratio",
                                                           "volume flux", "initial saturation water",
                                                           "block radius" ),
                                          "csmp::TwoPhaseModel_FractureMatrixUpscaled_Test",
                                          "nodal relative permeability oil",
                                          "nodal relative permeability water",
                                          "capillary pressure",
                                          20 ) );
  */

} // run



} // csmp

#endif // TWOPHASEMODEL_TESTS_H
