#include "CSMP_definitions.h"
#include "PropertyDatabase.h"
#include "Model1D.h"
#include "Region.h"
#include "InputDataManager.h"
#include "FlowFunctionsModule.h"
#include "TwoPhaseModelWithHysteresis_Test.h"
#include "BrooksCoreyWithHysteresis.h"


using namespace std;

namespace csmp {
  
/**
     1. the process test.. loop over saturation and check the process is picking the correct path
*/
vector <pair<double,double> >  TwoPhaseModelwithHysteresis_Test::Extract_data( Model<1U>& mdl, TestCases subfunctions)
 {
    
    vector< pair<double,double> > SwPc;
    
    FlowFunctionsModule2<1U> SatFunctions( mdl.Database(), 9.8 ) ;
    Region<1U>& mdl_domain{mdl.Region("Model")};
  
    for (auto it = mdl_domain.ElementsBegin(); it != mdl_domain.ElementsEnd(); ++it ) // loop over elements
    {
        double Sw = SatFunctions.Sw((*it)) ;
        double extract_data(0) ;
        
        switch (subfunctions) {
            
          case CapillaryPressure :
            extract_data = SatFunctions.pc((*it)) ;
            break;
            
          case DerivativeOfCapillaryPressure :
            extract_data = SatFunctions.dpcds((*it)) ;
            break;
            
          case krw :
            extract_data = SatFunctions.krw((*it)) ;
            break;
            
          case krn :
            extract_data = SatFunctions.krn((*it)) ;
            break;
            
          default:
            cout << "Error: None of Cases has been set" ;
            break;
        }
        
        SwPc.push_back(make_pair(Sw,extract_data)) ;
    }
    
    return SwPc;
    
} // end Extract_data
  
  
  


  
  
/**
 10. run over saturation from primary drainage to primary imbibitions..
 
 This is an Example how to have these parameters:
 Here are these parameters fitted to lab data from the paper  by:  Behrooz Raeesi, Norman R. Morrow, and Geoffrey Mason (2014)
 
 a) For Primary Drainage from the paper Raeesi Morrow Mason 2014 capillary hysteresis of Berea (see fig 3a)
 [1.1522295  0.40043515]
 [-0.01509935  0.87575084]
 
 b) For Drainage from the paper Raeesi Morrow Mason 2014 capillary hysteresis of Berea (see fig 3a)
 [0.6061638  0.63360703]
 [-1.6446064e-13  9.5839138e+00]
 
 c) For Imbib from the paper Raeesi Morrow Mason 2014 capillary hysteresis of Berea (see fig 3a)
 [0.31500375 0.8576149 ]
 [-5.1522971e-07  4.4824004e+00]
 
 array<array<double,2>,2> a1_ =  {{0.634,0.858},{9.58e+00,4.48}} ;
 array<array<double,2>,2> c1_ {{0.606,0.315},{-1.6e-13,-5.e-07}} ;
 
 array<array<double,2>,2> a1_ =  {{0.4,0.858},{0.9,4.48}} ;
 array<array<double,2>,2> c1_ {{1.15,0.315},{-0.015,-5.e-07}} ;
 
 array<array<double,2>,2> a1_ = {{{{0.5,0.25}},{{0.25,0.5}}}} ;  // these are just for test
 array<array<double,2>,2> c1_ = {{{{3.0,3.0}},{{-2.0,-2.0}}}} ;  // these are just for test
 */
 
void  TwoPhaseModelwithHysteresis_Test::runOverSaturationRange( TestCases local_case ) {
     //10.1 first create the Main Drainage and Imbibition curves...
  
  
    double length(1);
    size_t   elements(3);
  
    double DS(0.01);
  
    double srH2O(0.3) ;
    double srCO2(0.1) ;
  
    Point<1U> origin(0) ;
    Point<1U> destination(1);
  
    Model1D<1U> mdl( "Model1D", "CO2-geo-sequestration-variables.txt", length, elements );
  
    Region<1U>& model_domain = mdl.Region( "Model" );

    FlowFunctionsModule2<1U> SatFunctions( mdl.Database() , 9.8 ) ;
  
  
  
  /**
   
   10.2) Assign the storage keys and properties for the test.
   
   */
  
  // Build model and initialise with sensible saturation and parameter values
  
  // initialisation of model
  mdl.InputPropertyValue("porosity", makeScalar(PLAIN, 0.20)) ;
  mdl.InputPropertyValue("residual saturation aqueous phase", makeScalar(PLAIN, srH2O)) ;
  mdl.InputPropertyValue("residual saturation carbonic phase", makeScalar(PLAIN, srCO2)) ;
  mdl.InputPropertyValue("change of saturation carbonic phase", makeScalar(PLAIN, 0));
  
  mdl.InputPropertyValue("saturation aqueous phase", makeScalar(PLAIN, 0)) ;
  
  mdl.InputPropertyValue("fluid pressure", makeScalar(PLAIN, 50000000)) ;
  mdl.InputPropertyValue("temperature", makeScalar(PLAIN, 60)) ;
  mdl.InputPropertyValue("salinity", makeScalar(PLAIN, 0)) ;
  
  
  /**
   
   12.3.2) define the Hysteresis relative permeability model Parameters on each elements
   
   Important
   
   AWD, AOD, CWD, COD, AWI, AOI, CWI, COI  based on the paper by Skjaeveland et al. 2000
   Note: In "*-variable.txt" file, this has to be defined as:
   Hysteresis relative permeability model Parameters  HisBCparam  none  8  0.00E+00  5.00E+07  ELEMENT
   ^
   This number counts as a "DataDepth" for the array
   
   array<array<double,2>,2> a_ = {{{{0.5,0.25}},{{0.25,0.5}}}} ;
   array<array<double,2>,2> c_ = {{{{3.0,3.0}},{{-2.0,-2.0}}}} ;
   
   array<array<double,2>,2> a_ =  {{{{0.4,0.858}},{{0.9,4.48}}}} ;
   array<array<double,2>,2> c_ {{{{1.15,0.315}},{{-0.015,-5.e-07}}}} ;
  
  array<array<double,2>,2> a_ =  {{{{0.634,0.858}},{{9.58e+00,4.48}}}} ;
  array<array<double,2>,2> c_ =  {{{{0.606,0.315}},{{-1.6e-13,-5.e-07}}}} ;
  
  */
  
  array<array<double,2>,2> a_ = {{{{0.52,0.26}},{{0.24,0.4}}}} ;
  array<array<double,2>,2> c_ = {{{{5.0,3.0}},{{-2.0,-2.0}}}} ;
  
  
  ArrayVariable arrayvariable1_(8, PLAIN);
  
  typedef BrooksCoreySaturationFunctionsWithHysteresis<1,FlowFunctionsModule2> bcp;
  
  arrayvariable1_(bcp::AWD) = a_[bcp::H2O][bcp::DRAINAGE];
  arrayvariable1_(bcp::AOD) = a_[bcp::CO2][bcp::DRAINAGE];
  arrayvariable1_(bcp::CWD) = c_[bcp::H2O][bcp::DRAINAGE];
  arrayvariable1_(bcp::COD) = c_[bcp::CO2][bcp::DRAINAGE];
  
  arrayvariable1_(bcp::AWI) = a_[bcp::H2O][bcp::IMBIBITION];
  arrayvariable1_(bcp::AOI) = a_[bcp::CO2][bcp::IMBIBITION];
  arrayvariable1_(bcp::CWI) = c_[bcp::H2O][bcp::IMBIBITION];
  arrayvariable1_(bcp::COI) = c_[bcp::CO2][bcp::IMBIBITION];
  
  for ( auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it )
    (*it)->Store( SatFunctions.key_kri_param, arrayvariable1_);
  
  cout << "\n 1D Initialised value of the Hysteresis relative permeability model Parameters in the model." << endl;
  
  
  // Mahyar: this is all parameters we need for Drainage and Imbibition curves.. it is an array with 8 double number as ordered as:
  // AWD, AOD, CWD, COD, AWI, AOI, CWI, COI  based on the paper by Skjaeveland et al. 2000
  
  /**
   
   10.3) Start using the Saturation functions functionallity in the 1D model
   
   */
  
  
  /**
   
   12.3.1) accumulate the nodes with water satuaration from 0 to 1 ...
   
   */
  
  string filename(parseTestCases(local_case)) ;
  filename.append("_MainDrainageandImbibitions.csv") ;
  cout << filename << "\n" ;

  ofstream outputfile;
  outputfile.open (filename);
  outputfile<<"Sw,"<<parseTestCases(local_case)<<endl ;
  
  cerr <<"\n" ;
  cerr <<"The Main Drainage test....\n" ;
  cerr <<"\n" ;
  
  for ( double S = 0; S<= 1; S+=DS) {   // Drainage loop
    
    for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    {
      ScalarVariable Sw_ ((*it)->Status(SatFunctions.key_sH2O), S) ;
      ScalarVariable SCO2_ ((*it)->Status(SatFunctions.key_sCO2), 1-S) ;
      ScalarVariable pSCO2_ (SCO2_*.9) ;
      
      (*it)->Store( SatFunctions.key_sH2O,   Sw_   ) ;
      (*it)->Store( SatFunctions.key_sCO2,   SCO2_ ) ;
      (*it)->Store( SatFunctions.key_sCO2_0, pSCO2_) ;
    }
    
    auto SwPc1_st_Drainage = Extract_data(mdl,local_case) ;
    for (auto data : SwPc1_st_Drainage) outputfile << data.first <<" , "<< data.second <<endl  ;
    
  }
  
  for ( double S = 0; S<= 1; S+=DS) { //Imbibitaions loop
      
    for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    {
      ScalarVariable Sw_ ((*it)->Status(SatFunctions.key_sH2O), S) ;
      ScalarVariable SCO2_ ((*it)->Status(SatFunctions.key_sCO2), 1-S) ;
      ScalarVariable pSCO2_ (SCO2_*1.1) ;
      
      if (pSCO2_.Component(1) > 1) pSCO2_ = SCO2_ ;
      
      (*it)->Store( SatFunctions.key_sH2O,   Sw_   ) ;
      (*it)->Store( SatFunctions.key_sCO2,   SCO2_ ) ;
      (*it)->Store( SatFunctions.key_sCO2_0, pSCO2_) ;
    }
    
    auto SwPc1_st_Imbib = Extract_data(mdl,local_case) ;
    for (auto data : SwPc1_st_Imbib) outputfile << data.first <<" , "<< data.second <<endl  ;
    
  }
  
  outputfile.close() ;
  
  
  mdl.InputPropertyValue("change of saturation carbonic phase", makeScalar(PLAIN, 0));

  
  string filename2(parseTestCases(local_case)) ;
  filename2.append("_Drainage2Imbibitions.csv") ;
  cout << filename2 << "\n" ;

  ofstream outputfile2;
  outputfile2.open (filename2);
  outputfile2<<"Sw,"<<parseTestCases(local_case)<<endl ;

  cerr <<"\n" ;
  cerr <<"The Drainage to Imbibition test....\n" ;
  cerr <<"\n" ;
  
  
  ScalarVariable Sw0_  (ANY, 1.-srCO2) ;
  ScalarVariable PSCO2 (ANY, 1) ;

  DS = 0.001;
  
  bool index_im(false) ;
  
  for ( size_t time = 0; time < 1000 ; time++) {   // Drainage loop
    
    if (index_im) {
      Sw0_  = Sw0_       + DS; // make up imbibition
      PSCO2 = (1.- Sw0_) + DS;
    }
    else
    {
      Sw0_  = Sw0_       - DS; // make up drainage
      PSCO2 = (1.- Sw0_) - DS;
    }
    
    if ( abs( Sw0_.Component(1)-0.40 ) <= DS ) index_im = true ;  // imbib one
    if ( abs( Sw0_.Component(1)-0.5 ) <= DS ) index_im = false ;  // drainge one

    
    if ( Sw0_.Component(1)  < srH2O        ) break ;
    if ( Sw0_.Component(1)  > 1.-srCO2     ) break ;

    
    
    for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    {
      (*it)->Store( SatFunctions.key_sH2O,   Sw0_   ) ;
      (*it)->Store( SatFunctions.key_sCO2,   1.-Sw0_ ) ; // assume SCO2 = 1 - SH2O
      (*it)->Store( SatFunctions.key_sCO2_0, PSCO2  ) ;
    }
    
    auto SwPc1_st_Drainage2Imbib = Extract_data(mdl,local_case) ;
    for (auto data : SwPc1_st_Drainage2Imbib) outputfile2 << data.first <<" , "<< data.second <<endl  ;
    
  }
  
  outputfile2.close() ;
    
} // end runOverSaturationRange

  
  
  
  
  

/**
 
 12) This is main file to run...
 
 */
void csmp::TwoPhaseModelwithHysteresis_Test::run() {
  /**
   
   12.1) Create an 1D model without any spliting point in the model.
   
   */
  double length   = 10;
  size_t   elements = 20;
  
  double srH2O(0.) ;
  double srCO2(0.) ;
  
  Point<1U> origin = 0 ;
  Point<1U> destination = 10;
  
  Model1D<1U> model( "Model1D", "CO2-geo-sequestration-variables.txt", length, elements );
  printModelDimensions( model, true );
  cout << "\n 1D Model created" << endl;
  
  Region<1U>&              model_domain = model.Region( "Model" );
  FlowFunctionsModule2<1U> SatFunctions( model.Database() , 9.8 ) ;

  /**
   
   12.2) Assign the storage keys and properties for the test.
   
   */
  
  // Build model and initialise with sensible saturation and parameter values
  
  // initialisation of model
  model.InputPropertyValue("porosity", makeScalar(PLAIN, 0.20)) ;
  model.InputPropertyValue("residual saturation aqueous phase", makeScalar(PLAIN, srH2O)) ;
  model.InputPropertyValue("residual saturation carbonic phase", makeScalar(PLAIN, srCO2)) ;
  
  model.InputPropertyValue("saturation aqueous phase", makeScalar(PLAIN, 0)) ;
  
  model.InputPropertyValue("fluid pressure", makeScalar(PLAIN, 50000000)) ;
  model.InputPropertyValue("temperature", makeScalar(PLAIN, 60)) ;
  model.InputPropertyValue("salinity", makeScalar(PLAIN, 0)) ;
  
  model.InputPropertyValue("model time", makeScalar(PLAIN, 0)) ;
    
    

    
  // Mahyar: this is all parameters we need for Drainage and Imbibition curves.. it is an array with 8 double number as ordered as:
  // AWD, AOD, CWD, COD, AWI, AOI, CWI, COI  based on the paper by Skjaeveland et al. 2000
  
  /**
   
   12.3) Start using the Saturation functions functionallity in the 1D model
   
   */
  
  
  
  
  /**
   
   12.3.1) accumulate the nodes with water satuaration from 0 to 1 ...
   
   */
  double i(0.0) ;
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it )
    {
      ScalarVariable variable0_ ((*it)->Status(SatFunctions.key_sH2O), i/elements) ;
      (*it)->Store( SatFunctions.key_sH2O, variable0_);
      (*it)->Store( SatFunctions.key_sCO2, 1.0-variable0_) ; // assume SCO2 = 1 - SH2O
      (*it)->Store( SatFunctions.key_sCO2_0,1.0-variable0_) ;
      i++ ;
    }
    
  /**
   
   12.3.2) define the Hysteresis relative permeability model Parameters on each elements
   
   
   
   Important
   
   AWD, AOD, CWD, COD, AWI, AOI, CWI, COI  based on the paper by Skjaeveland et al. 2000
   Note: In "*-variable.txt" file, this has to be defined as:
   Hysteresis relative permeability model Parameters  HisBCparam  none  8  0.00E+00  5.00E+07  ELEMENT
   ^
   This number counts as a "DataDepth" for the array
   
   array<array<double,2>,2> a_ = {{{{0.5,0.25}},{{0.25,0.5}}}} ;
   array<array<double,2>,2> c_ = {{{{3.0,3.0}},{{-2.0,-2.0}}}} ;
   
   array<array<double,2>,2> a_ =  {{{{0.4,0.858}},{{0.9,4.48}}}} ;
   array<array<double,2>,2> c_ {{{{1.15,0.315}},{{-0.015,-5.e-07}}}} ;
   
   
   */
  
  array<array<double,2>,2> a_ =  {{{{0.634,0.858}},{{9.58e+00,4.48}}}} ;
  array<array<double,2>,2> c_ =  {{{{0.606,0.315}},{{-1.6e-13,-5.e-07}}}} ;
  
  
  ArrayVariable arrayvariable1_(8, PLAIN);
  
    typedef BrooksCoreySaturationFunctionsWithHysteresis<1,FlowFunctionsModule2> bcp;
  
  arrayvariable1_(bcp::AWD) = a_[bcp::H2O][bcp::DRAINAGE];
  arrayvariable1_(bcp::AOD) = a_[bcp::CO2][bcp::DRAINAGE];
  arrayvariable1_(bcp::CWD) = c_[bcp::H2O][bcp::DRAINAGE];
  arrayvariable1_(bcp::COD) = c_[bcp::CO2][bcp::DRAINAGE];
  
  arrayvariable1_(bcp::AWI) = a_[bcp::H2O][bcp::IMBIBITION];
  arrayvariable1_(bcp::AOI) = a_[bcp::CO2][bcp::IMBIBITION];
  arrayvariable1_(bcp::CWI) = c_[bcp::H2O][bcp::IMBIBITION];
  arrayvariable1_(bcp::COI) = c_[bcp::CO2][bcp::IMBIBITION];
  
  for ( auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it )
    (*it)->Store( SatFunctions.key_kri_param, arrayvariable1_);
  
  cout << "\n 1D Initialised value of the Hysteresis relative permeability model Parameters in the model." << endl;
  
  
  /**
   
   12.3.3) printing out water satuaration of each nodes
   
   */
  
  for ( auto it = model_domain.NodesBegin(); it != model_domain.NodesEnd(); ++it ) {
    auto Sw = (*it)->Read(SatFunctions.key_sH2O) ;
    cerr <<"Water Saturation defined at Nodes: "<< Sw <<"\n" ;
    
  }
  
  // outputing into csv file of the Sw , Pc , dPcdS , krw , krn , dkrwds , dkrnds , fw , Sh , Sv
  
  ofstream outputfile;
  outputfile.open ("Capillary_Pressure_correlations_for_Mixed_wet_Reservoirs.csv");
  outputfile<<"Sw,Pc,dPcdS,krw,krn,dkrwds,dkrnds"<<endl ;
  
  
  /**
   
   12.3.4) loop over each element and assign the process and constants to calculate the
   Pc , dPcdS , krw , krn , dkrwds , dkrnds , fw , Sh , Sv and etc.
   
   */
  
  for (auto it = model_domain.ElementsBegin(); it != model_domain.ElementsEnd(); ++it ) // loop over elements
  {
    auto idx_ = (*it)->Idx() ;
    
    /**
     Important:
     
     Mahyar: here I define contants for imbibitions and drainage
     for example for primary dranage aw = 0.5, ao = 0.25, cw =3.0, co=0.0.
     for the second dranage aw = 0.5, ao = 0.25, cw =3.0, co=-2.0.
     for the second imbibitions aw = 0.25, ao = 0.5, cw = 3.0, co=-2.0.
     
     */
    
    
    // Water Saturation of model at the BaryCenter.. this is not the same as water saturation of Nodes.. its is interploated from the nodes at BaryCenter
    
  
    double Sw = SatFunctions.Sw((*it)) ;
    
    cerr <<"\n" ;
    cerr <<"Entering into Element: " << idx_ << "\n" ;
    cerr <<"Water Saturation at the BaryCenter: "<< Sw <<"\n" ;
    
    // water Saturation
    // ---------------------------------------------
    _test( !( Sw < 0. || Sw > 1. ) );
    
    
    //  Check the water saturation is in the range of
   
    double pc = SatFunctions.pc((*it)) ;
    
    cerr <<"Capillary pressure: "<< pc <<"\n" ;
    
    // checking capilary pressure
    // ---------------------------------------------
    //_test( !isnan(pc) );
    
    double dpcds = SatFunctions.dpcds((*it)) ;
    
    cerr <<"Derivative of capillary pressure: "<< dpcds <<"\n" ;
    
    // checking 1st derivative of capilary pressure
    // ---------------------------------------------
    _test( !isnan(dpcds) );
    
    double krw = SatFunctions.krw((*it)) ;
    
    cerr <<"Relative permability of Water: "<< krw <<"\n" ;
    
    // checking relative permeability of water
    // ---------------------------------------------
    _test( !isnan(krw) );
    _test(!( krw < 0. || krw > 1. )) ;
    
    double dkrwds = SatFunctions.dkrwds((*it)) ;
    
    cerr <<"Derivative of Relative permability of Water: "<< dkrwds <<"\n" ;
    
    // checking 1st derivative of relative permeability of water
    // ---------------------------------------------
    _test( !isnan(dkrwds) );
    
    double krn = SatFunctions.krn((*it)) ;
    
    cerr <<"Relative permability of Co2: "<< krn <<"\n" ;
    
    // checking relative permeability of CO2
    // ---------------------------------------------
    _test( !isnan(krn) );
    _test(!( krn < 0. || krn > 1. )) ;
    
    double dkrnds = SatFunctions.dkrnds((*it)) ;
    
    cerr <<"Derivative of Relative permability of Co2: "<< dkrnds <<"\n" ;
    
    
 
    
    /// Test has been done in all functions for Flow Functions
    
    
    
    
    outputfile<<Sw<<","<<pc<<","<<dpcds<<","<<krw<<","<<krn<<","<<dkrwds<<","<<dkrnds<<endl ; // otherwise we can not visulise with python
    
    cerr <<"\n" ;
    
    
    cerr << "============================================================================\n" ;
    cerr << "    Entering into Prescribed saturation loop, to test the functionality  \n" ;
    cerr << "============================================================================\n" ;
    
    
    // Functionality for the prescribed saturation
    // ---------------------------------------------
    for (double S(0.0); S <= 1.0; S+=0.1){
      
      cerr <<"\n" ;
      cerr <<"Prescribed saturation:  "<< S <<"\n" ;
      
      double S_at = SatFunctions.EffectiveSaturation_at((*it), S) ;
      
      // water Saturation
      // ---------------------------------------------
      _test(!( S_at < 0. || S_at > 1.  ));
      
      krw = SatFunctions.krw_at((*it), S) ;
      
      cerr <<"Relative permability of Water: "<< krw <<"\n" ;
      
      // checking relative permeability of water
      // ---------------------------------------------
      _test( !isnan(krw) );
      _test(!( krw < 0. || krw > 1. )) ;
      
      dkrwds = SatFunctions.dkrwds_at((*it), S) ;
      
      cerr <<"Derivative of Relative permability of Water: "<< dkrwds <<"\n" ;
      
      // checking 1st derivative of relative permeability of water
      // ---------------------------------------------
      _test( !isnan(dkrwds) );
      
      krn = SatFunctions.krn_at((*it), S) ;
      
      cerr <<"Relative permability of Co2: "<< krn <<"\n" ;
      
      // checking relative permeability of CO2
      // ---------------------------------------------
      _test( !isnan(krn) );
      _test(!( krn < 0. || krn > 1. )) ;
      
      dkrnds = SatFunctions.dkrnds_at((*it), S) ;
      
      cerr <<"Derivative of Relative permability of Co2: "<< dkrnds <<"\n" ;
      
      
      // checking 1st derivative of relative permeability of CO2
      // ---------------------------------------------
      _test( !isnan(dkrnds) );
      
      
    }
    
    cerr << "============================================================================\n" ;
    cerr << "    Finishing Prescribed saturation loop  \n" ;
    cerr << "============================================================================\n" ;
    
    
    
    cerr <<"\n" ;
    cerr << Sw<<" , "<<pc<<endl ;
    
    
  } // end loop over elements
  outputfile.close() ;
  
  
  cerr << "\n The Unit Test has been done ....\n" ;
  
  //
  //  The unit test has been done.. the following of code is going to run over different saturation
  //  and create an example of turing point in saturations curve....
  //
  
  cerr << "============================================================================\n" ;
  cerr << "        Start run over saturations for dranage and imbibition loops         \n" ;
  cerr << "============================================================================\n" ;
  
  runOverSaturationRange(CapillaryPressure) ;     // Pass
  
  runOverSaturationRange(krw) ;                    // Pass
  
  runOverSaturationRange(krn) ;                   // Pass
  
} // run()

} // csmp
