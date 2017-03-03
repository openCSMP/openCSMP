#include "RelativePermeabilityModel_Test.h"

#include "ErrorHandler.h" 
#include "ModelTopology.h"
#include "InputDataManager.h"
#include "Region.h"
#include "Model.h"
#include "Model1D.h"
#include "LineElementMesher.h"

// new polymorphic relative permeability classes
#include "LinearTwoPhaseModel.h"
#include "BrooksCorey.h"
#include "VanGenuchten.h"
#include "TwoPhaseFileBased.h"
//#include "ToddQuadraticModel.h"
#include "Experimental2PhaseModel.h"
#include "GenericTransferFunction.h"
#include "FractureMatrixUpscaled.h"
#include "FourarLenormand.h"

// auxiliary functions called by this main program
#include "PropertyHandle.h"
#include "CSMP_highLevelUtilities.h"


using namespace std;

#define DIM 1U

namespace csmp {

RelativePermeabilityModel_Test::RelativePermeabilityModel_Test()
 : model_ptr_(0)
 {
    const uint32  N_ELEMENTS(100);
    model_ptr_ = new Model1D<1U>("Model1D", "CSMP-2phase-variables_upscaled.txt", 50., N_ELEMENTS );
    cerr << "Number of elemnts: " << model_ptr_ -> Mesh().Elements() << endl;
    // making some groups: rock (elements 1-40, 61-100) and fracture (elements 41-60)
    vector<size_t>   elms;  elms.reserve( N_ELEMENTS );
    for ( uint32 i=0; i<40; i++ ) elms.push_back(i);
    for ( uint32 i=60; i<N_ELEMENTS; i++ ) elms.push_back(i);
    model_ptr_->FormRegionFrom( "ROCK", elms );
    elms.erase( elms.begin(), elms.end() );
    for ( uint32 i=40; i<60; i++ ) elms.push_back(i);
    model_ptr_->FormRegionFrom( "FRACTURE", elms );
 }


RelativePermeabilityModel_Test::~RelativePermeabilityModel_Test()
 {
    delete model_ptr_;
 }





void RelativePermeabilityModel_Test::run()
 {
   InitializeFlowProperties();
   
   TwoPhaseModel<DIM>*  relperm_model(0);

   // Linear model
   // ------------
   try {
        // value of zero invokes linear model
        model_ptr_->InputPropertyValue( "brooks corey parameter", makeScalar(PLAIN,0.) );
        relperm_model = new BrooksCorey<DIM>( model_ptr_->Database(), 
                                             "brooks corey parameter", "entry pressure" );
        Test( *relperm_model, true );
        delete relperm_model;
        model_ptr_->InputPropertyValue( "brooks corey parameter", makeScalar(PLAIN,2.) );
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in LinearTwoPhaseModel");
     }
   _succeed();

   // vanGenuchten Model
   // -------------------

   try {
        relperm_model = new VanGenuchten<DIM>( model_ptr_->Database(),
                                              "van genuchten parameter", "van genuchten alpha" );
        Test( *relperm_model, true );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in vanGenuchten");
     }
   _succeed();


   // Todd quadratic model
   // --------------------
/*
   try {
        relperm_model = new ToddQuadraticModel<DIM>( model_ptr_->Database() );
        Test( *relperm_model, true );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in ToddQuadraticModel");
     }
   _succeed();
*/
   // Brooks-Corey model
   // ------------------

   try {
        relperm_model = new BrooksCorey<DIM>( model_ptr_->Database(),
                                             "brooks corey parameter", "entry pressure" );
        Test( *relperm_model, true );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in BrooksCorey");
     }
   _succeed();


   // Experimental2PhaseModel model
   // -----------------------------
   try {
        // write some experimental relative permeability curves
        /*
            two-column datafile with the format:
            header line - just filename and explanation what's tabulated
            blank line
            number of x,y value pairs (the tabulated function)
            derivate1 and derivative_n at the function origin and endpoint
            x-y value pairs (one per line) 1..n 
        */
        
        ofstream  krw_data("krw_data.txt");
        assert( krw_data.is_open() );
        krw_data <<"krw_data for testing of Experimental2PhaseModel"<< endl << endl;
        krw_data << 5U << endl;
        krw_data <<"0.5 1." << endl;
        krw_data <<"0. 0.\n0.2 0.1\n0.4 0.3\n0.6 0.6\n0.8 0.9\n1. 1."<< endl;
        krw_data.close();

        ofstream  kro_data("kro_data.txt");
        assert( kro_data.is_open() );
        kro_data <<"kro_data for testing of Experimental2PhaseModel"<< endl << endl;
        kro_data << 5U << endl;
        kro_data <<"-1. -1." << endl;
        kro_data <<"0. 1.\n0.2 0.8\n0.4 0.6\n0.6 0.4\n0.8 0.2\n1. 0."<< endl;
        kro_data.close();

        ofstream  pc_data("pc_data.txt");
        assert( pc_data.is_open() );
        pc_data <<"pc_data for testing of Experimental2PhaseModel"<< endl << endl;
        pc_data << 5U << endl;
        pc_data <<"-1000. -10000." << endl;
        pc_data <<"0. 5.0e7\n0.2 1.0e6\n0.4 1.0e5\n0.6 1.0e4\n0.8 0.5e4\n1. 0."<< endl;
        pc_data.close();

        relperm_model = new Experimental2PhaseModel<DIM>( model_ptr_->Database(),
                                                                  "krw_data", "kro_data", "pc_data" );
        Test( *relperm_model, true );
        delete relperm_model;

     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in Experimental2PhaseModel");
     }
   _succeed();


   // FractureMatrixUpscaled model
   // ----------------------------
   // parameters employed in the upscaling
   PropertyHandle<1U>  Af( *model_ptr_, "fracture matrix interface area", SCALAR, ELEMENT );
   PropertyHandle<1U>  phi_f( *model_ptr_, "fracture porosity", SCALAR, ELEMENT );
   PropertyHandle<1U>  qfqm( *model_ptr_, "fracture matrix flux ratio", SCALAR, ELEMENT );
   PropertyHandle<1U>  R( *model_ptr_, "block radius", SCALAR, NODE );
   PropertyHandle<1U>  swi( *model_ptr_, "initial saturation water", SCALAR, NODE );
   
   model_ptr_->InputPropertyValue( "fracture matrix interface area",  makeScalar(PLAIN,2.0) );
   model_ptr_->InputPropertyValue( "fracture porosity",               makeScalar(PLAIN,0.1) );
   model_ptr_->InputPropertyValue( "fracture matrix flux ratio",      makeScalar(PLAIN,5.0) );
   model_ptr_->InputPropertyValue( "initial saturation water",        makeScalar(PLAIN,0.01) ); // as needed for transfer function
   model_ptr_->InputPropertyValue( "volume flux",                     makeScalar(PLAIN,1.0e-9) ); // for scaling of transfer contribution
   model_ptr_->InputPropertyValue( "block radius",                    makeScalar(PLAIN,0.7) );
   try {
      relperm_model = new FractureMatrixUpscaled<DIM>( model_ptr_->Database(),
                                                                "permeability", "viscosity oil", "viscosity water",
                                                                "density oil", "density water", 
                                                                "brooks corey parameter", "entry pressure",
                                                                "fracture matrix interface area",
                                                                "fracture porosity", "porosity",
                                                                "fracture matrix flux ratio",
                                                                "volume flux", "initial saturation water",
                                                                "block radius" );

        Test( *relperm_model, false );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in FractureMatrixUpscaled");
     }
   _succeed();


// FOURAR LENORMAND
// ----------------------------
   try {
        PropertyHandle<1U>  a( *model_ptr_, "fracture aperture", SCALAR, ELEMENT );
        model_ptr_->InputPropertyValue( "fracture aperture",  makeScalar(PLAIN,1.0e-3) );
        randomPerturb( *model_ptr_, "fracture aperture", 20. ); // by 20 percent
        relperm_model = new FourarLenormand<DIM>( model_ptr_->Database(), "fracture aperture" );
        Test( *relperm_model, true );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in FourarLenormand");
     }
   _succeed();


// TWO PAHSE FILE BASED
// ----------------------------
   try {
        relperm_model = new TwoPhaseFileBased<DIM>( model_ptr_->Database(), "RelPerms.txt" );
        Test( *relperm_model, true );
        delete relperm_model;
     }
   catch(...) {
        _fail("RelativePermeabilityModel_Test::run: Error in TwoPhaseFileBased");
     }
   _succeed();




   
 } // end run





void RelativePermeabilityModel_Test::InitializeFlowProperties()
 {
    VectorVariable<DIM>   Vd(DIRICH,1.0e-10); // Vd(DIRICH,0.135); 
    model_ptr_->InputPropertyValue( "velocity", Vd );
    model_ptr_->InputPropertyValue( "permeability",                          makeScalar(PLAIN,1.0e-12) );
    model_ptr_->InputPropertyValue( "porosity",                              makeScalar(PLAIN,0.2) );
    model_ptr_->InputPropertyValue( "fluid pressure",                        makeScalar(PLAIN,5.0e7) );
    model_ptr_->InputPropertyValue( "viscosity water",                       makeScalar(PLAIN,1.0e-03) );
    model_ptr_->InputPropertyValue( "viscosity oil",                         makeScalar(PLAIN,3.0e-03) );
    model_ptr_->InputPropertyValue( "density water",                         makeScalar(PLAIN,1000.) );
    model_ptr_->InputPropertyValue( "density oil",                           makeScalar(PLAIN,800.) );
    model_ptr_->InputPropertyValue( "van genuchten parameter",               makeScalar(PLAIN,3.) );
    model_ptr_->InputPropertyValue( "van genuchten alpha",                   makeScalar(PLAIN,0.37) );
    model_ptr_->InputPropertyValue( "brooks corey parameter",                makeScalar(PLAIN,2.) );
    model_ptr_->InputPropertyValue( "residual saturation wetting phase",     makeScalar(PLAIN,0.2) );
    model_ptr_->InputPropertyValue( "residual saturation non-wetting phase", makeScalar(PLAIN,0.25) );
    model_ptr_->InputPropertyValue( "entry pressure",                        makeScalar(PLAIN,2.0e3) );
 }



/**
     To run test manually with a user-specified relative permeability model.
*/
void  RelativePermeabilityModel_Test::Test( TwoPhaseModel<1U>& relperm,
                                            bool extended_property_set )                      
 { 

    const csmp::Index   satw_key = model_ptr_->Database().StorageKey("saturation water");
    const csmp::Index   satn_key = model_ptr_->Database().StorageKey("saturation oil");
    const double64      sat_incr(1./model_ptr_->Mesh().Nodes());
    ScalarVariable  saturation;
    
    // memorizing the original saturation values

    model_ptr_->CopyReplace( "saturation oil", "previous saturation oil" );
    Region<1U>& sg(model_ptr_->Region("Model"));
    // generating a range of saturation values for water and oil
    for ( vector<Node<1U>*>::iterator
          it=sg.NodesBegin(); it!=sg.NodesEnd(); it++ )
      {
         saturation() = 0. + sat_incr * (*it)->Idx();
         (*it)->Store( satw_key, saturation );
         saturation = 1. - saturation;
         (*it)->Store( satn_key, saturation );
      }
        
    /// @todo (2-C) Rm rtti (entire file, lots of refs)
    // creating an output file name
    string  relperm_model_name( typeid(relperm).name() );
    // exception for LinearTwoPhaseModel
    double64 lmin, lmax;
    model_ptr_->MinMaxOf("brooks corey parameter", lmin, lmax );
    if ( lmax <= numeric_limits<double64>::epsilon() ) relperm_model_name = typeid(LinearTwoPhaseModel<1U>).name();
    // removing not permitted characters 
    for ( string::iterator it=relperm_model_name.begin(); it!=relperm_model_name.end(); it++ )
         if ( !isalpha(*it) and !isdigit(*it) ) (*it) = '_'; 
    relperm_model_name += ".txt";
       
    // setting up the output file   
    ofstream  ofs( relperm_model_name.c_str() );
    bool      is_nan_test_failed(false), value_out_of_range(false);
    ofs <<"'"<< relperm_model_name <<"' written by testRelativePermeabilityModel(";
    ofs << typeid( relperm ).name() <<") for display in Excel."<< endl;
    if ( extended_property_set ) {
        ofs <<"sw\tseff\tkrw\tkrn\tmob_t\tfn\tdfn/dSn ";
        ofs <<"\tG\tdGdS\tpc\tdpcdSn\tadvection-mult\tdiffusion-mult\tgravity-G-mult\tgravity-dGds-mult\tpc-diff-mult\tshock speed"<< endl;
      }
    else ofs <<"sw\tseff\tkrw\tkrn\tmob_t\tfn\tG\tpc\tdiffusion-mult\tgravity-G-mult"<< endl;

    // computing multiphase flow properties and writing these to file
    for ( vector<Element<1U>*>::iterator 
          it=sg.ElementsBegin(); it!=sg.ElementsEnd(); it++ )
      {

         // setting up the relative permeability model
         // ---------------------------------------------
         relperm.Initialize( *(*it) );
         relperm.InitializeForNode( *(*it), 1U );
         relperm.EffectiveSaturation();
        
         /// @todo (2-C) Rm rtti
         // water saturation & effective water saturation
         // ---------------------------------------------
         const double64 sw(relperm.Saturation(1U));
         if ( sw < 0. or sw > 1. ) {
              cerr <<"\nRelativePermeabilityModel_Test::Test:(";
              cerr << typeid( relperm ).name() <<") sw out of range. Testing cannot be performed."<< endl;
              return;
           }
         ofs << relperm.Saturation(1U) <<"\t"; 

         if ( relperm.EffectiveSaturation() < 0. or 1. < relperm.EffectiveSaturation() ) {
              cerr <<"\n"<< typeid( relperm ).name() <<" seff("<< sw <<") out of range."<< endl;
              value_out_of_range = true;
              _test(false);
           }
         ofs << relperm.EffectiveSaturation() <<"\t"; 
        
         // 1. test relative permeability
         // -------------------------------------------------------------------------------------------
         if ( relperm.krw_Phase() < 0. or relperm.krw_Phase() > 1. ) {
              cerr <<"\n"<< typeid( relperm ).name() <<"::krw("<< sw <<") out of range: "<< relperm.krw_Phase() << endl;
              value_out_of_range = true;
              _test(false);
           }
         if ( relperm.krn_Phase() < 0. or relperm.krn_Phase() > 1. ) {
              cerr <<"\n"<< typeid( relperm ).name() <<"::krn("<< sw <<") out of range: "<< relperm.krn_Phase() << endl;
              value_out_of_range = true;
              _test(false);
           }
         ofs << relperm.krw_Phase() <<"\t"<< relperm.krn_Phase() <<"\t";
         
         // 2. total mobility
         // -------------------------------------------------------------------------------------------
         ofs << relperm.TotalMobility() <<"\t";

         // 3. fractional flow and its saturation derivative
         // -------------------------------------------------------------------------------------------
         if ( relperm.f_Phase(2U) < 0. or relperm.f_Phase(2U) > 1. ) {
              cerr <<"\n"<< typeid( relperm ).name() <<"::f_Phase("<< sw <<") out of range: "<< relperm.f_Phase(2U) << endl;
              value_out_of_range = true;
              _test(false);
           }
         ofs << relperm.f_Phase( 2U ) <<"\t";
         
         if ( extended_property_set ) { 
              if ( (is_nan_test_failed=isnan(relperm.dfds())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::dfds("<< sw <<") is NaN (erratic value)."<< endl;
                   _test(false);
                }
              ofs << relperm.dfds() <<"\t";
           }
           
         // 4. G and derivative of G-function
         // -------------------------------------------------------------------------------------------
         if ( (is_nan_test_failed=isnan(relperm.G())) ) {
              cerr <<"\n"<< typeid( relperm ).name() <<"::G("<< sw <<") is NaN (erratic value)."<< endl;
              _test(false);
           }
         ofs << relperm.G() <<"\t";

         if ( extended_property_set ) { 
              if ( (is_nan_test_failed=isnan(relperm.dGds())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::dGds("<< sw <<") is NaN (erratic value)."<< endl;
                   _test(false);
                }
              ofs << relperm.dGds() <<"\t";
           }

         // 5. capillary pressure and its saturation derivative
         // -------------------------------------------------------------------------------------------
         if ( (is_nan_test_failed=isnan(relperm.pc_Phase())) or relperm.pc_Phase() > 1e9 ) { // 1GPa
              cerr <<"\n"<< typeid( relperm ).name() <<"::pc_Phase("<< sw <<") is NaN (erratic value) or greater than 1 GPa."<< endl;
              _test(false);
           }
         ofs << relperm.pc_Phase() <<"\t";

         if ( extended_property_set ) { 
              if ( (is_nan_test_failed=isnan(relperm.dpcds_Phase())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::dpcds_Phase("<< sw <<") is NaN (erratic value)."<< endl;
                   _test(false);
                }
              ofs << relperm.dpcds_Phase() <<"\t";
           }
         
         // 6. the 3 characteristic multipliers: advection, diffusion, gravity
         // -------------------------------------------------------------------------------------------
         if ( extended_property_set ) {
              // uses fractional flow derivative  
              if ( (is_nan_test_failed=isnan(relperm.AdvectionMultiplier())) ) {
                    cerr <<"\n"<< typeid( relperm ).name() <<"::AdvectionMultiplier("<< sw <<") is NaN (erratic value)."<< endl;
                    _test(false);
                }
           }
         if ( (is_nan_test_failed=isnan(relperm.DiffusionMultiplier( 2U )) ) ) {
              cerr <<"\n"<< typeid( relperm ).name() <<"::DiffusionMultiplier("<< sw <<") is NaN (erratic value)."<< endl;
              _test(false);
           }
         if ( (is_nan_test_failed=isnan(relperm.GravityMultiplier_G())) ) {
               cerr <<"\n"<< typeid( relperm ).name() <<"::GravityMultiplier_G("<< sw <<") is NaN (erratic value)."<< endl;
               _test(false);
           }
           
         if ( extended_property_set ) ofs << relperm.AdvectionMultiplier() <<"\t";
         ofs << relperm.DiffusionMultiplier( 2U ) <<"\t";
         ofs << relperm.GravityMultiplier_G() <<"\t";

         if ( extended_property_set ) {  
              if ( (is_nan_test_failed=isnan(relperm.GravityMultiplier_dGds())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::GravityMultiplier_dGds("<< sw <<") is NaN (erratic value)."<< endl;
                   _test(false);
                }
              if ( (is_nan_test_failed=isnan(relperm.CapillaryDiffusionMultiplier())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::CapillaryDiffusionMultiplier("<< sw <<") is NaN (erratic value)."<< endl;
                   _test(false);
                }
              if ( (is_nan_test_failed=isnan(relperm.ShockSpeed())) ) {
                   cerr <<"\n"<< typeid( relperm ).name() <<"::ShockSpeed() is NaN (erratic value)."<< endl;
                   _test(false);
                }

              ofs << relperm.GravityMultiplier_dGds() <<"\t";
              ofs << relperm.CapillaryDiffusionMultiplier() <<"\t";
              ofs << relperm.ShockSpeed(); 
           }
         ofs << endl;
      }
      
    ofs.close();

    // restoring the original saturation values

    model_ptr_->CopyReplace( "previous saturation oil", "saturation oil" );
    for ( vector<Node<1U>*>::iterator 
          it=sg.NodesBegin(); it!=sg.NodesEnd(); it++ ) {
         saturation() = 1. - (*it)->Read( satn_key );
         (*it)->Store( satw_key, saturation );
      }


    if ( value_out_of_range ) {
         cerr <<"\nrelative permeability model: '"<< typeid(relperm).name() <<"' Error."<< endl;
         throw out_of_range("RelativePermeabilityModel_Test::Test: relperm model produces out of range values.");
      }
    if ( is_nan_test_failed ) {
         cerr <<"\nrelative permeability model: '"<< typeid(relperm).name() <<"' Error."<< endl;
         throw range_error("RelativePermeabilityModel_Test::Test: relperm model produces NaN output.");
      }
    cout <<"\nRelativePermeabilityModel_Test::Test: file '"<< relperm_model_name <<"' written successfully."<< endl; 
 
 } // end Test



const PropertyDatabase<1>& RelativePermeabilityModel_Test::Database() const
 { return model_ptr_->Database(); }



} // end csmp

