#include "CSMP_definitions.h"

#include "Model1D.h"
#include "InputDataManager.h"

#include "FlowFunctionsBC_Hysteretic.h"
#include "BrooksCoreySaturationFunctionswithHysteresis.h"

#include "TwoPhaseModelWithHysteresis_Test.h"


using namespace std;

namespace csmp {
  
  // default constructor
  TwoPhaseModelwithHysteresis_Test::TwoPhaseModelwithHysteresis_Test()
  {
  }
 
  
  
  
  
  // simple destructore
  TwoPhaseModelwithHysteresis_Test::~TwoPhaseModelwithHysteresis_Test()
  {
  }
  
  
  
  
  
  
  
  
/**
 1. the process test.. loop over saturation and check the process is picking the correct path
*/
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::Process_Test(Model<1U>& mdl) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      auto E1 = (*it)->AtBarycenter() ;
      SatFunctions.InitialiseBrooksCoreyParameters(E1) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      bool flag_  = (Sw > Sw_min)&&(Sw < Sw_max) ;
      
      if (flag_) {
        
        double64 pc = SatFunctions.pc(E1) ;
        SwPc.push_back(make_pair(Sw,pc)) ;
        
      }
    }
    
    return SwPc ;
}
  
  
  
  
  
  
  
  
  
/**
2.  Imbibition test... loop over saturation and calculate Pc
*/
  
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::Imbibition_Test(Model<1U>& mdl, const array<array<double64,2>,2>&  a_, const array<array<double64,2>,2>&  c_) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    TWO_PHASE_FLOW_PROCESS Process_path = IMBIBITION ;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      auto E1 = (*it)->AtBarycenter() ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ , Process_path) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      bool flag_  = (Sw > Sw_min)&&(Sw < Sw_max) ;
      
      if (flag_) {
        
        double64 pc = SatFunctions.pc(E1) ;
        SwPc.push_back(make_pair(Sw,pc)) ;
        
      }
    }
    
    return SwPc ;
}

  
  
  
  

  
  
/**
   3.  Drainage test... loop over saturation and calculate Pc
*/
  
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::Drainage_Test(Model<1U>& mdl, const array<array<double64,2>,2>&  a_, const array<array<double64,2>,2>&  c_) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    TWO_PHASE_FLOW_PROCESS Process_path = DRAINAGE ;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      
      auto E1 = (*it)->AtBarycenter() ;
      
      SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ , Process_path) ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      bool flag_  = (Sw > Sw_min)&&(Sw < Sw_max) ;
      
      if (flag_) {
        
        double64 pc = SatFunctions.pc(E1) ;
        SwPc.push_back(make_pair(Sw,pc)) ;
        
      }
    }
    
    return SwPc ;
}

  
  
  
  
  
  
  
  
  
/**
 4. Update water and oil residual saturations
*/
void TwoPhaseModelwithHysteresis_Test::UpdatePseduoResidualSaturations(Model<1U>& mdl , double64 Srw , double64 Sro ){
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    
    for ( auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ){
      
      ScalarVariable variable1_ ((*it)->Status(var.key_psrH2O), Srw) ;
      (*it)->Store( var.key_psrH2O, variable1_);
      
      ScalarVariable variable2_ ((*it)->Status(var.key_psrCO2), Sro) ;
      (*it)->Store( var.key_psrCO2, variable2_);
      
    }
}
  
  
  
  
  
  
  
  
  
  
  
/**
 5. The transition from Drainage to Imbibition and reverse test... loop over saturation and calculate Pc
*/
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::TransitionProcess_Test(Model<1U>& mdl) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    
    double64  tol_(0.01) ;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    double64 newSro(0) ;
    double64 newSrw(0) ;
    
    TWO_PHASE_FLOW_PROCESS ProcessPath = DRAINAGE ;
    
    
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      
      auto E1 = (*it)->AtBarycenter() ;
      
      double64 psrCO2 = (*it)->Read(var.key_psrCO2) ;
      double64 psrH2O = (*it)->Read(var.key_psrH2O) ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      SatFunctions.InitialiseBrooksCoreyParameters(E1) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
       //, Sw_min, Sw_max
      
      TWO_PHASE_FLOW_PROCESS ProcessPath = SatFunctions.FlowProcess( E1 );
      
      bool flag_(false) ;
      
      switch ( ProcessPath ) {
        case DRAINAGE:   // that means current is Dranaige and Sw_min will change to Imbibitions
          flag_ = (abs(Sw-Sw_min)/Sw_min < tol_) ;
          
          if ( flag_ ) {
            
            newSro = SatFunctions.OilResidualSaturation(E1) ;
            cerr <<"At saturation: "<< Sw <<"\n" ;
            cerr <<"old oil residual pseudo saturation: "<< psrCO2 <<"\n" ;
            cerr <<"new oil residual pseudo saturation: "<< newSro <<"\n" ;
            cerr <<"\n" ;
            
            newSrw = SatFunctions.WaterResidualSaturation(E1, newSro) ;
            
            cerr <<"old water residual pseudo saturation: "<< psrH2O <<"\n" ;
            cerr <<"new water residual pseudo saturation: "<< newSrw <<"\n" ;
            cerr <<"\n" ;
            
          }
          break;
          
        case IMBIBITION:
           flag_ = (abs(Sw-Sw_max)/Sw_max < tol_) ;
          
          if ( flag_ ) {
            
            newSrw = psrH2O ;
            newSro = psrCO2 ;
            
            SatFunctions.WaterAndOilResidualSaturationImbibitionToDrainage(E1, newSrw, newSro) ;  // Mahyar: Here the oil saturation will change...
            
            cerr <<"At saturation: "<< Sw <<"\n" ;
            cerr <<"old oil residual pseudo saturation: "<< psrCO2 <<"\n" ;
            cerr <<"new oil residual pseudo saturation: "<< newSro <<"\n" ;
            cerr <<"old water residual pseudo saturation: "<< psrH2O <<"\n" ;
            cerr <<"new water residual pseudo saturation: "<< newSrw <<"\n" ;
            cerr <<"\n" ;
            
          }
          break;
          
          default:
          cout<<"Error: The Transition neither is Drainage nor Imbibitions.."<<endl;
          break;
          
      }
    }
    
    UpdatePseduoResidualSaturations(mdl, newSrw , newSro);
    
    switch ( ProcessPath ) {
        
      case DRAINAGE:
         CreateArtifitialImbibitionProcess(mdl);
        break;
        
      case IMBIBITION:
        CreateArtifitialDrainageProcess(mdl);
        break;
        
    }
   
    SwPc = Process_Test(mdl) ;
    

    return SwPc ;
}
  
  
  
  
  
  
  
  
  


  
/**
 6. The Drainage To Imbibition test... loop over saturation and calculate Pc
*/
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::DrainageToImbibition_Test(Model<1U>& mdl, const array<array<double64,2>,2>&  a_, const array<array<double64,2>,2>&  c_) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    
    double64  tol_(0.01) ;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    double64 newSro(0) ;
    double64 newSrw(0) ;
    
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      
      auto E1 = (*it)->AtBarycenter() ;
      
      double64 psrCO2 = (*it)->Read(var.key_psrCO2) ;
      double64 psrH2O = (*it)->Read(var.key_psrH2O) ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      TWO_PHASE_FLOW_PROCESS Process_path = DRAINAGE ;
      
      SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ , Process_path) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      bool flag_ = (abs(Sw-Sw_min)/Sw_min < tol_) ;
      
      if ( flag_ ) {
        
        newSro = SatFunctions.OilResidualSaturation(E1) ;
        cerr <<"At saturation: "<< Sw <<"\n" ;
        cerr <<"old oil residual Pseduo saturation: "<< psrCO2 <<"\n" ;
        cerr <<"new oil residual Pseduo saturation: "<< newSro <<"\n" ;
        cerr <<"\n" ;
        
        newSrw = SatFunctions.WaterResidualSaturation(E1, newSro) ;
        
        cerr <<"old water residual Pseduo saturation: "<< psrH2O <<"\n" ;
        cerr <<"new water residual Pseduo saturation: "<< newSrw <<"\n" ;
        cerr <<"\n" ;
        
      }
    }
    
    
    UpdatePseduoResidualSaturations(mdl , newSrw, newSro);
    
    SwPc = Imbibition_Test(mdl, a_, c_) ;
    
    return SwPc ;
}

  
  
  
  
  
  
  
  
  
  
  
/**
   7. The Imbibition to Drainage test... loop over saturation and calculate Pc
*/
vector <pair<double64,double64> >  TwoPhaseModelwithHysteresis_Test::ImbibitionToDrainage_Test(Model<1U>& mdl, const array<array<double64,2>,2>&  a_, const array<array<double64,2>,2>&  c_) {
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    vector< pair<double64,double64> > SwPc;
    
    double64  tol_(0.01) ;
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    double64 newSro(0) ;
    double64 newSrw(0) ;
    
    
    for (auto it = mdl.Region("Model").ElementsBegin(); it != mdl.Region("Model").ElementsEnd(); ++it ) // loop over elements
      
    {
      
      auto E1 = (*it)->AtBarycenter() ;
      
      double64 psrCO2 = (*it)->Read(var.key_psrCO2) ;
      double64 psrH2O = (*it)->Read(var.key_psrH2O) ;
      
      TWO_PHASE_FLOW_PROCESS Process_path = IMBIBITION ;
      SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ , Process_path) ;
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      double64 Sw_min = E1.Obtain(var.key_SwDrToImb) ;
      double64 Sw_max = E1.Obtain(var.key_SwImbToDr) ;
      
      
      bool flag_ = (abs(Sw-Sw_max)/Sw_max < tol_) ;
      
      if ( flag_ ) {
        
        newSrw = psrH2O ;
        newSro = psrCO2 ;
        
        SatFunctions.WaterAndOilResidualSaturationImbibitionToDrainage(E1, newSrw, newSro) ;  // Mahyar: Here the oil saturation will change...
        
        cerr <<"old oil residual Pseduo saturation: "<< psrCO2 <<"\n" ;
        cerr <<"new oil residual Pseduo saturation: "<< newSro <<"\n" ;
        cerr <<"old water residual Pseduo saturation: "<< psrH2O <<"\n" ;
        cerr <<"new water residual Pseduo saturation: "<< newSrw <<"\n" ;
        cerr <<"\n" ;
        
      }
    }
    
    
    UpdatePseduoResidualSaturations(mdl , newSrw, newSro);
    
    SwPc = Drainage_Test(mdl, a_, c_) ;
    
    return SwPc ;
}
 
  
  
  
  
  
  
  
  
  
  
  
  
/**
 8. This is a function make a parameter files for Primary Drainage and pure Imbibitions
*/
void   Primary(double64 c_[2][2], TWO_PHASE_FLOW_PROCESS Process_path) {
    
    switch (Process_path) {
      case DRAINAGE:
        c_[CO2][DRAINAGE]   = 0.0 ;
        break;
        
      case IMBIBITION:
        c_[H2O][IMBIBITION] = 0.0 ;
        break;
    }
    
}

  

  
  
  
  
  
  
  
  
  
  
  
/**
   9. This is a function check the saturation range between Srw, and Sro
*/
bool TwoPhaseModelwithHysteresis_Test::Out_of_Bond(double64 Sw, double64 srH2O, double64 srCO2)
  {
    return ( (Sw < srH2O || Sw > 1.-srCO2 )) ;
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
/**
 10. run over saturation from primary drainage to primary imbibitions..
*/
void  TwoPhaseModelwithHysteresis_Test::runOverSaturationRange(Model<1U>& mdl) {
    
    
    /**
     
     10.1 first create the Main Drainage and Imbibition curves...
     
     */
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions=FlowFunctionsBC_Hysteretic<1U>(mdl.Database()) ;
    
    ofstream outputfile;
    outputfile.open ("MainDrainageandImbibitions.csv");
    outputfile<<"Sw , Pc "<<endl ;
    
/**
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
    
       array<array<double64,2>,2> a1_ =  {{0.634,0.858},{9.58e+00,4.48}} ;
       array<array<double64,2>,2> c1_ {{0.606,0.315},{-1.6e-13,-5.e-07}} ;
  
       array<array<double64,2>,2> a1_ =  {{0.4,0.858},{0.9,4.48}} ;
       array<array<double64,2>,2> c1_ {{1.15,0.315},{-0.015,-5.e-07}} ;
 
       array<array<double64,2>,2> a1_ = {{{{0.5,0.25}},{{0.25,0.5}}}} ;  // these are just for test
       array<array<double64,2>,2> c1_ = {{{{3.0,3.0}},{{-2.0,-2.0}}}} ;  // these are just for test
 */
  
    array<array<double64,2>,2> a1_ =  {{{{0.634,0.858}},{{9.58e+00,4.48}}}} ;
    array<array<double64,2>,2> c1_ =  {{{{0.606,0.315}},{{-1.6e-13,-5.e-07}}}} ; // in psi
    
    cerr <<"\n" ;
    cerr <<"The Main Drainage test....\n" ;
    cerr <<"\n" ;
    
    mdl.InputPropertyValue("previous drainage saturation endpoint of aqueous phase", makeScalar(PLAIN, 0)) ;  //update the starting Imbibitions point
    mdl.InputPropertyValue("previous imbibition saturation endpoint of aqueous phase", makeScalar(PLAIN, 1)) ; //update the ending  Imbibitions point
    
    CreateArtifitialDrainageProcess(mdl);
    auto SwPc1_st_Drainage = Process_Test(mdl) ;
    for (auto data : SwPc1_st_Drainage) outputfile << data.first <<" , "<< data.second <<endl  ;
    
    cerr <<"\n" ;
    cerr <<"Then Main Imbibition test....\n" ;
    cerr <<"\n" ;
    
    
    CreateArtifitialImbibitionProcess(mdl);
    auto SwPc1_st_Imbibition = Process_Test(mdl) ;
    for (auto data : SwPc1_st_Imbibition) outputfile << data.first <<" , "<< data.second <<endl  ;
    
    outputfile.close() ;
    
    /**
     
     10.2 Second do transition from Drainage to Imbibitions
     
     */
    
    ofstream outputfile2;
    outputfile2.open ("FirstPrimaryDrainageandImbibitions.csv");
    outputfile2<<"Sw , Pc "<<endl ;
    
    
    cerr <<"\n" ;
    cerr <<"Then first Dranaige test....\n" ;
    cerr <<"\n" ;
    
    mdl.InputPropertyValue("previous drainage saturation endpoint of aqueous phase", makeScalar(PLAIN, 0.45)) ;  //update the starting Imbibitions point
    mdl.InputPropertyValue("previous imbibition saturation endpoint of aqueous phase", makeScalar(PLAIN, 0.6)) ; //update the ending  Imbibitions point
    
    
    CreateArtifitialDrainageProcess(mdl);
    
    auto SwPc0 = Process_Test(mdl) ;
    for (auto data : SwPc0) outputfile2 << data.first <<" , "<< data.second <<endl  ;
    
    cerr <<"\n" ;
    cerr <<"Then transition from Dranaige to Imbibition test....\n" ;
    cerr <<"\n" ;
    
    outputfile2.close() ;
    ofstream outputfile3;
    outputfile3.open ("TransitionsfromPrimaryDrainageandImbibitions.csv");
    outputfile3<<"Sw , Pc "<<endl ;
    
    CreateArtifitialDrainageProcess(mdl);
  
    //auto SwPc1 = TransitionProcess_Test(mdl) ;
    
    auto SwPc1 = DrainageToImbibition_Test(mdl, a1_, c1_) ;
    for (auto data : SwPc1) outputfile3 << data.first <<" , "<< data.second <<endl  ;
    
    cerr <<"\n" ;
    cerr <<"Then transition from Imbibition to Drainage test....\n" ;
    cerr <<"\n" ;
    
    CreateArtifitialImbibitionProcess(mdl);
  
    //auto SwPc2 = TransitionProcess_Test(mdl) ;
    
    auto SwPc2 = ImbibitionToDrainage_Test(mdl, a1_ , c1_) ;
    for (auto data : SwPc2) outputfile3 << data.first <<" , "<< data.second <<endl  ;
    
    outputfile3.close() ;
    
}
  
  
  
  
  
  
  
  
  
/**
   11) create drainage process
 */
void TwoPhaseModelwithHysteresis_Test::CreateArtifitialDrainageProcess(Model<1U>& mdl){
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    for ( auto it = mdl.Region("Model").NodesBegin(); it != mdl.Region("Model").NodesEnd(); ++it )
    {
      double64 SCO2 = (*it)->Read(var.key_sCO2) ;
      double64 S = (*it)->Read(var.key_SwDrToImb) ;
      
      if ( SCO2 < 1.-S)
      {
        ScalarVariable variable1_((*it)->Status(var.key_sCO2_1), SCO2*1.1) ;
        (*it)->Store(var.key_sCO2_1, variable1_) ;
      }
    }
}
  
  
  
  
  
  
  
  
  
  
  
/**
 11) create imbibition process
*/
void TwoPhaseModelwithHysteresis_Test::CreateArtifitialImbibitionProcess(Model<1U>& mdl){
    
    const PropertyDatabase<1>&   p_ref = mdl.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    for ( auto it = mdl.Region("Model").NodesBegin(); it != mdl.Region("Model").NodesEnd(); ++it )
    {
      double64 SCO2 = (*it)->Read(var.key_sCO2) ;
      double64 S = (*it)->Read(var.key_SwImbToDr) ;
      if ( SCO2 > 1.-S)
      {
        ScalarVariable variable1_((*it)->Status(var.key_sCO2_1), SCO2*0.9) ;
        (*it)->Store(var.key_sCO2_1, variable1_) ;
      }
    }
}
  
  
  
  
  
  
  
  
  
  

  
  
  
/**
 
 12) This is main file to run...
 
*/
void csmp::TwoPhaseModelwithHysteresis_Test::run() {
    /**
     
     12.1) Create an 1D model without any spliting point in the model.
     
     */
    double64 length   = 10;
    size_t   elements = 1000;
    
    double64 srH2O(0.3) ;
    double64 srCO2(0.35) ;
    
    double64 water_mu(0.001) ;
    double64 CO2_mu(0.00002);
  
    Point<1U> origin = 0 ;
    Point<1U> destination = 10;
    
    Model1D<1U> model( "Model1D", "TwoPhaseModelwithHysteresis_Test-variables.txt", length, elements );
    Region<1U>& rref = model.Region( "Model" );
    const PropertyDatabase<1>&   p_ref = model.Database();
    variables::VariableSet_CO2GeoSequestration var( p_ref );
    
    printModelDimensions( model, true );
    cout << "\n 1D Model created" << endl;
    
    /**
     
     12.2) Assign the storage keys and properties for the test.
     
     */
    
    // Build model and initialise with sensible saturation and parameter values
    
    // initialisation of model
    model.InputPropertyValue("porosity", makeScalar(PLAIN, 0.20)) ;
    model.InputPropertyValue("viscosity carbonic phase", makeScalar(PLAIN, CO2_mu)) ;
    model.InputPropertyValue("viscosity aqueous phase", makeScalar(PLAIN, water_mu)) ;
    
    model.InputPropertyValue("residual saturation aqueous phase", makeScalar(PLAIN, srH2O)) ;
    model.InputPropertyValue("residual saturation carbonic phase", makeScalar(PLAIN, srCO2)) ;
    
    
    //  Brooks Corey model
    model.InputPropertyValue("entry pressure", makeScalar(PLAIN, 10000)) ;
    model.InputPropertyValue("brooks corey parameter", makeScalar(PLAIN, 2.)) ;
    
    // Brooks Corey Hysteresis model
    
    model.InputPropertyValue("previous drainage saturation endpoint of aqueous phase", makeScalar(PLAIN, srH2O)) ;
    model.InputPropertyValue("previous imbibition saturation endpoint of aqueous phase", makeScalar(PLAIN, 1.-srCO2)) ;
    
    model.InputPropertyValue("pseudo residual saturation aqueous phase", makeScalar(PLAIN, srH2O)) ;
    model.InputPropertyValue("pseudo residual saturation carbonic phase", makeScalar(PLAIN, srCO2)) ;
    
    model.InputPropertyValue("saturation aqueous phase", makeScalar(PLAIN, 0)) ;
  
    // Mahyar: this is all parameters we need for Drainage and Imbibition curves.. it is an array with 8 double number as ordered as:
    // AWD, AOD, CWD, COD, AWI, AOI, CWI, COI  based on the paper by Skjaeveland et al. 2000
    
    /**
     
     12.3) Start using the Saturation functions functionallity in the 1D model
     
     */
    
    FlowFunctionsBC_Hysteretic<1U> SatFunctions = FlowFunctionsBC_Hysteretic<1U>(model.Database()) ;
    
    
    /**
     
     12.3.1) accumulate the nodes with water satuaration from 0 to 1 ...
     
     */
    double64 i(0.0) ;
    for ( auto it = rref.NodesBegin(); it != rref.NodesEnd(); ++it )
    {
      ScalarVariable variable0_ ((*it)->Status(var.key_sH2O), i/elements) ;
      (*it)->Store( var.key_sH2O, variable0_);
      (*it)->Store( var.key_sCO2, 1.0-variable0_) ; // assume SCO2 = 1 - SH2O
      (*it)->Store( var.key_sCO2_1,1.0-variable0_) ;
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
    
      array<array<double64,2>,2> a_ = {{{{0.5,0.25}},{{0.25,0.5}}}} ;
      array<array<double64,2>,2> c_ = {{{{3.0,3.0}},{{-2.0,-2.0}}}} ;
  
      array<array<double64,2>,2> a_ =  {{{{0.4,0.858}},{{0.9,4.48}}}} ;
      array<array<double64,2>,2> c_ {{{{1.15,0.315}},{{-0.015,-5.e-07}}}} ;
     
     
     */
  
       array<array<double64,2>,2> a_ =  {{{{0.634,0.858}},{{9.58e+00,4.48}}}} ;
       array<array<double64,2>,2> c_ =  {{{{0.606,0.315}},{{-1.6e-13,-5.e-07}}}} ;
  
    
    ArrayVariable arrayvariable1_(8, PLAIN);
    
    arrayvariable1_(AWD) = a_[H2O][DRAINAGE];
    arrayvariable1_(AOD) = a_[CO2][DRAINAGE];
    arrayvariable1_(CWD) = c_[H2O][DRAINAGE];
    arrayvariable1_(COD) = c_[CO2][DRAINAGE];
    
    arrayvariable1_(AWI) = a_[H2O][IMBIBITION];
    arrayvariable1_(AOI) = a_[CO2][IMBIBITION];
    arrayvariable1_(CWI) = c_[H2O][IMBIBITION];
    arrayvariable1_(COI) = c_[CO2][IMBIBITION];
    
   
    
    for ( auto it = rref.ElementsBegin(); it != rref.ElementsEnd(); ++it )
    {
      (*it)->Store( var.key_HisBCparam, arrayvariable1_);
    }
    
    cout << "\n 1D Initialised value of the Hysteresis relative permeability model Parameters in the model." << endl;
    
    
    /**
     
     12.3.3) printing out water satuaration of each nodes
     
     */
    
    for ( auto it = rref.NodesBegin(); it != rref.NodesEnd(); ++it )
    {
      
      auto Sw = (*it)->Read(var.key_sH2O) ;
      cerr <<"Water Saturation defined at Nodes: "<< Sw <<"\n" ;
      
    }
    
    // outputing into csv file of the Sw , Pc , dPcdS , krw , krn , dkrwds , dkrnds , fw , Sh , Sv
    
    cout << "\n 1D Checking values of Water Saturation, Mobility, " << endl;
    
    ofstream outputfile;
    outputfile.open ("Capillary_Pressure_correlations_for_Mixed_wet_Reservoirs.csv");
    outputfile<<"Sw , Pc , dPcdS , krw , krn , dkrwds , dkrnds , fw , Sh , Sv"<<endl ;
    
    
    /**
     
     12.3.4) loop over each element and assign the process and constants to calculate the
     Pc , dPcdS , krw , krn , dkrwds , dkrnds , fw , Sh , Sv and etc.
     
     */
    
    for (auto it = rref.ElementsBegin(); it != rref.ElementsEnd(); ++it ) // loop over elements
    {
      auto idx_ = (*it)->Idx() ;
      auto E1 = (*it)->AtBarycenter() ;  // Mahyar: Never, ever delete this line... The Finite Element Placement of E1  is important...
      
      
      /**
       Important:
       
        Mahyar: here I define contants for imbibitions and drainage
        for example for primary dranage aw = 0.5, ao = 0.25, cw =3.0, co=0.0.
        for the second dranage aw = 0.5, ao = 0.25, cw =3.0, co=-2.0.
        for the second imbibitions aw = 0.25, ao = 0.5, cw = 3.0, co=-2.0.
       
       
        There is other possibilities to initialise parameters.. check the header file of BrooksCoreySaturationFunctionsWithHysteresis.h
       
        SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ , Drainage) ;
        SatFunctions.InitialiseBrooksCoreyParameters(E1, a_ , c_ ) ; // Mahyar: to Initialaise Drainage process .. otherwise check the new CO2 saturation to to allocate the correct process
       
       */
      
      SatFunctions.InitialiseBrooksCoreyParameters(E1) ;
      
      
      // Water Saturation of model at the BaryCenter.. this is not the same as water saturation of Nodes.. its is interploated from the nodes at BaryCenter
      
      double64 Sw = SatFunctions.Sw(E1) ;
      
      cerr <<"\n" ;
      cerr <<"Entering into Element: " << idx_ << "\n" ;
      cerr <<"Water Saturation at the BaryCenter: "<< Sw <<"\n" ;
      
      // water Saturation
      // ---------------------------------------------
      _test( !( Sw < 0. || Sw > 1. ) );
      
      
      //  Check the water saturation is in the range of
      
      double64 pc = SatFunctions.pc(E1) ;
      
      cerr <<"Capillary pressure: "<< pc <<"\n" ;
      
      // checking capilary pressure
      // ---------------------------------------------
      //_test( !isnan(pc) );
      
      double64 dpcds = SatFunctions.dpcds(E1) ;
      
      cerr <<"Derivative of capillary pressure: "<< dpcds <<"\n" ;
      
      // checking 1st derivative of capilary pressure
      // ---------------------------------------------
      _test( !isnan(dpcds) );
      
      double64 krw = SatFunctions.krw(E1) ;
      
      cerr <<"Relative permability of Water: "<< krw <<"\n" ;
      
      // checking relative permeability of water
      // ---------------------------------------------
      _test( !isnan(krw) );
      _test(!( krw < 0. || krw > 1. )) ;
      
      double64 dkrwds = SatFunctions.dkrwds(E1) ;
      
      cerr <<"Derivative of Relative permability of Water: "<< dkrwds <<"\n" ;
      
      // checking 1st derivative of relative permeability of water
      // ---------------------------------------------
      _test( !isnan(dkrwds) );
      
      double64 krn = SatFunctions.krn(E1) ;
      
      cerr <<"Relative permability of Co2: "<< krn <<"\n" ;
      
      // checking relative permeability of CO2
      // ---------------------------------------------
      _test( !isnan(krn) );
      _test(!( krn < 0. || krn > 1. )) ;
      
      double64 dkrnds = SatFunctions.dkrnds(E1) ;
      
      cerr <<"Derivative of Relative permability of Co2: "<< dkrnds <<"\n" ;
      
      
      // checking 1st derivative of relative permeability of CO2
      // ---------------------------------------------
      _test( !isnan(dkrnds) );
      
      
      double64 lambda_w = SatFunctions.Mobility(E1, 0) ;
      
      cerr <<"Mobility of Water: "<< lambda_w <<"\n" ;
      
      // water Mobility
      // ---------------------------------------------
      _test( !isnan(lambda_w) );
      _test( !( lambda_w < 0. || lambda_w > (1./water_mu)) );   // Mahyar: upper limits of water mobility is 1/mu
      
      double64 dlambda_w = SatFunctions.MobilityDerivative(E1, 0) ;
      
      cerr <<"1st derivative of Mobility of Water: "<< dlambda_w <<"\n" ;
      
      // checking 1st derivative of water Mobility
      // ---------------------------------------------
      _test( !isnan(dlambda_w) );
      
      
      double64 tlambda = SatFunctions.TotalMobility(E1) ;
      
      cerr <<"Total Mobility of Water and Co2: "<< tlambda <<"\n" ;
      
      // checking 1st derivative of total Mobility
      // ---------------------------------------------
      _test( !isnan(tlambda) );
      _test( !( tlambda < 0. || tlambda > (1./water_mu + 1./CO2_mu) ) ); // Mahyar: upper limits of total mobility is sum of 1/mu
      
      
      double64 Plambda = SatFunctions.MobilityProduct(E1) ;
      
      cerr <<"Mobility product of Water and Co2: "<< Plambda <<"\n" ;
      
      // checking 1st product of Mobilities
      // ---------------------------------------------
      _test( !isnan(Plambda) );
      _test( !( Plambda < 0. || Plambda > (1./water_mu * 1./CO2_mu) ) ); // Mahyar: upper limits of Product of mobility is product of 1/mu
      
      
      double64 dPlambda =  SatFunctions.MobilityProductDerivative(E1) ;
      
      cerr <<"1st derivative of Mobility product of Water: "<< dPlambda <<"\n" ;
      
      double64 Flowfraction = SatFunctions.f(E1, 0) ;
      
      cerr <<"frcational flow function of water: "<< Flowfraction <<"\n" ;
      
      // checking Flow fraction
      // ---------------------------------------------
      _test( !isnan(Flowfraction) );
      _test( !( Flowfraction < 0. || Flowfraction > 1. ) );
      
      double64 FirstDerivative_FlowFunctions = SatFunctions.dfds(E1, 0) ;
       _test( !isnan(FirstDerivative_FlowFunctions) );
      
      cerr <<"dfds of Water: "<< FirstDerivative_FlowFunctions <<"\n" ;
      
      double64 S_inf = SatFunctions.InflectionPointSaturation(E1) ;
      
      cerr <<"Infelection point: "<< S_inf <<"\n" ;
      
      // checking Infelection point:
      // ---------------------------------------------
      _test( !isnan(S_inf) );
      _test( !( S_inf < 0. || S_inf > 1. ) );
      
      double64 MaxFlowfraction = SatFunctions.MaxFractionalFlowDerivative(E1) ;
      _test( !isnan(MaxFlowfraction) );

      cerr <<"Max Flow fraction: "<< MaxFlowfraction <<"\n" ;
      
      double64 V_sh = SatFunctions.ShockSpeed(E1) ;
      
      cerr <<"Shock Speed: "<< V_sh <<"\n" ;
      
      // checking Speed of shock
      // ---------------------------------------------
      _test( !isnan(V_sh) );
      
      double64 S_sh = SatFunctions.ShockHeight(E1) ;
      
      cerr <<"Shock Saturations (Height): "<< S_sh <<"\n" ;
      
      // checking height of shock; saturation of shock
      // ---------------------------------------------
      _test( !isnan(S_sh) );
      _test( !( S_sh < 0. || S_sh > 1. ) );
    
      
      /**
      double64 GT = SatFunctions.GravityTerm(E1) ;  // should be check in 2D or 3D model
      cerr <<"Gravity Term : "<< GT <<"\n" ;
      
      // checking GravityTerm
      // ---------------------------------------------
      _test( !isnan(GT) );
    
      
      double64 GM = SatFunctions.GravityMultiplier_G(E1) ;
      
      cerr <<"Gravity Multiplier : "<< GM <<"\n" ;
      
      // checking Gravity Multiplier
      // ---------------------------------------------
      _test( !isnan(GM) );
      
      double64 dGM = SatFunctions.GravityMultiplier_dGds(E1, false) ;
      
      cerr <<"Gravity Multiplier derivative : "<< dGM <<"\n" ;
      
      // checking Gravity Multiplier derivative
      // ---------------------------------------------
      _test( !isnan(dGM) );
       
      
      
      
      
      double64 DM_water = SatFunctions.DiffusionMultiplier(E1, 0U) ;
      
      cerr <<"Diffusion Multiplier  of water phase: "<< DM_water <<"\n" ;
      
      // checking Capillary Diffusion Multiplier
      // ---------------------------------------------
      _test( !isnan(DM_water) );
       
  
      
      double64 CDM = SatFunctions.CapillaryDiffusionMultiplier(E1) ;
      
      cerr <<"Capillary Diffusion Multiplier: "<< CDM <<"\n" ;
      
      // checking Gravity Multiplier derivative
      // ---------------------------------------------
      _test( !isnan(CDM) );
       
       */
      
      
      double64 S_tan = SatFunctions.TangentPointSaturation(E1) ;
      
      cerr <<"Tangent Point Saturation : "<< S_tan <<"\n" ;
      
      // checking Tangent Point Saturation:
      // ---------------------------------------------
      _test( !isnan(S_tan) );
      _test( !( S_tan < 0. || S_tan > 1. ) );
      
      
      double64 V_f = SatFunctions.ShockFrontVelocity(E1) ;
      
      cerr <<"ShockFrontVelocity: "<< V_f <<"\n" ;
      
      // checking Speed of shock
      // ---------------------------------------------
      _test( !isnan(V_f) );
      
      
      /// Test has been done in all functions for Flow Functions
      
      
      
      
      outputfile<<Sw<<" , "<<pc<<" , "<<dpcds<<" , "<<krw<<" , "<<krn<<" , "<<dkrwds<<" , "<<dkrnds<<" , "<<Flowfraction<<" , "<<S_sh<<" , "<<V_sh<<endl ; // otherwise we can not visulise with python
      
      cerr <<"\n" ;
      
      
      cerr << "============================================================================\n" ;
      cerr << "    Entering into Prescribed saturation loop, to test the functionality  \n" ;
      cerr << "============================================================================\n" ;
      
      
      // Functionality for the prescribed saturation
      // ---------------------------------------------
      for (double64 S(0.0); S <= 1.0; S+=0.1){
        
        cerr <<"\n" ;
        cerr <<"Prescribed saturation:  "<< S <<"\n" ;
      
        double64 S_at = SatFunctions.EffectiveSaturation_at(E1, S) ;
        
        // water Saturation
        // ---------------------------------------------
        _test(!( S_at < 0. || S_at > 1.  ));
        
        double64 krw = SatFunctions.krw_at(E1, S) ;
        
        cerr <<"Relative permability of Water: "<< krw <<"\n" ;
        
        // checking relative permeability of water
        // ---------------------------------------------
        _test( !isnan(krw) );
        _test(!( krw < 0. || krw > 1. )) ;
        
        double64 dkrwds = SatFunctions.dkrwds_at(E1, S) ;
        
        cerr <<"Derivative of Relative permability of Water: "<< dkrwds <<"\n" ;
        
        // checking 1st derivative of relative permeability of water
        // ---------------------------------------------
        _test( !isnan(dkrwds) );
        
        double64 krn = SatFunctions.krn_at(E1, S) ;
        
        cerr <<"Relative permability of Co2: "<< krn <<"\n" ;
        
        // checking relative permeability of CO2
        // ---------------------------------------------
        _test( !isnan(krn) );
        _test(!( krn < 0. || krn > 1. )) ;
        
        double64 dkrnds = SatFunctions.dkrnds_at(E1, S) ;
        
        cerr <<"Derivative of Relative permability of Co2: "<< dkrnds <<"\n" ;
        
        
        // checking 1st derivative of relative permeability of CO2
        // ---------------------------------------------
        _test( !isnan(dkrnds) );
        
        
        double64 lambda_w = SatFunctions.Mobility_at(E1, 0U, S) ;
        
        cerr <<"Mobility of Water: "<< lambda_w <<"\n" ;
        
        // water Mobility
        // ---------------------------------------------
        _test( !isnan(lambda_w) );
        _test( !( lambda_w < 0. || lambda_w > (1./water_mu)) );   // Mahyar: upper limits of water mobility is 1/mu
        
        double64 dlambda_w = SatFunctions.MobilityDerivative_at(E1, 0U, S) ;
        
        cerr <<"1st derivative of Mobility of Water: "<< dlambda_w <<"\n" ;
        
        // checking 1st derivative of water Mobility
        // ---------------------------------------------
        _test( !isnan(dlambda_w) );
        
        
        double64 tlambda = SatFunctions.TotalMobility_at(E1, S) ;
        
        cerr <<"Total Mobility of Water and Co2: "<< tlambda <<"\n" ;
        
        // checking 1st derivative of total Mobility
        // ---------------------------------------------
        _test( !isnan(tlambda) );
        _test( !( tlambda < 0. || tlambda > (1./water_mu + 1./CO2_mu) ) ); // Mahyar: upper limits of total mobility is sum of 1/mu
        
        
        
        double64 dPlambda =  SatFunctions.MobilityProductDerivative_at(E1, S) ;
        
        cerr <<"1st derivative of Mobility product of Water: "<< dPlambda <<"\n" ;
        
        double64 Flowfraction = SatFunctions.f_at(E1, 0U, S) ;
        
        cerr <<"frcational flow function of water: "<< Flowfraction <<"\n" ;
        
        // checking Flow fraction
        // ---------------------------------------------
        _test( !isnan(Flowfraction) );
        _test( !( Flowfraction < 0. || Flowfraction > 1. ) );
        
        double64 FirstDerivative_FlowFunctions = SatFunctions.dfds_at(E1, Sw) ;
        
        cerr <<"dfds of Water: "<< FirstDerivative_FlowFunctions <<"\n" ;
        
      
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
    //  and create an example of turing bpoint in saturations curve....
    //
    
    cerr << "==========================================================\n" ;
    cerr << "         Start run over saturations  \n" ;
    cerr << "==========================================================\n" ;
    
    runOverSaturationRange(model) ;
    
    
} // run()
  
} // csmp
