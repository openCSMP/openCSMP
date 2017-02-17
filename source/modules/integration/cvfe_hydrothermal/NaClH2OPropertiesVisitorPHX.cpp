#include <cmath>
#include <algorithm>

#include "NaClH2OPropertiesVisitorPHX.h"
#include "ConvertConcentrationUnitsNaCl.h"
#include "CompareFloats.h"


using namespace std;

namespace csmp
{

  template<size_t dim>
  NaClH2OPropertiesVisitorPHX<dim>::NaClH2OPropertiesVisitorPHX( Model<dim>& model )
  //
#include "NaClH2OPropertiesVisitorPHX_initializer_list.hpp"
  //
  {
    // Initialise INDEX variables
    t_key = pref.StorageKey("temperature");
    p_key = pref.StorageKey("fluid pressure");
    hCl_key = pref.StorageKey("enthalpy content liquid");
    hCv_key = pref.StorageKey("enthalpy content vapor");
    tp_key  = pref.StorageKey("previous temperature");
    phi_key = pref.StorageKey("nodal porosity");
    cpr_key = pref.StorageKey("nodal heat capacity rock");
    ncp_key = pref.StorageKey("nodal heat capacity");
    rr_key  = pref.StorageKey("nodal density rock");
    hClp_key  = pref.StorageKey("previous enthalpy content liquid");
    hCvp_key  = pref.StorageKey("previous enthalpy content vapor");
    cpv_key = pref.StorageKey("pore volume");
    mf_key  = pref.StorageKey("fluid mass");
    ml_key  = pref.StorageKey("fluid mass liquid");
    mv_key  = pref.StorageKey("fluid mass vapor");
    mlp_key = pref.StorageKey("previous fluid mass liquid");
    mvp_key = pref.StorageKey("previous fluid mass vapor");
    mt_key  = pref.StorageKey("fluid density");
    mtp_key = pref.StorageKey("previous fluid density");
    rho_bulk_key  = pref.StorageKey("bulk fluid density");
    Htp_key = pref.StorageKey("previous total enthalpy");
    hf_key  = pref.StorageKey("fluid enthalpy");
    sl_key  = pref.StorageKey("saturation liquid");
    sv_key  = pref.StorageKey("saturation vapor");
    rl_key  = pref.StorageKey("density liquid");
    rv_key  = pref.StorageKey("density vapor");
    mul_key = pref.StorageKey("viscosity liquid");
    muv_key = pref.StorageKey("viscosity vapor");
    hl_key  = pref.StorageKey("enthalpy liquid");
    hv_key  = pref.StorageKey("enthalpy vapor");
    cpf_key = pref.StorageKey("fluid heat capacity");
    beta_key  = pref.StorageKey("compressibility");
    //*** new TD May 2011
    beta_p_key  = pref.StorageKey("previous compressibility");
    beta_ref_key  = pref.StorageKey("reference compressibility");
    apc_key = pref.StorageKey("after phasechange counter");
    dpc_key = pref.StorageKey("dangerous phase change");
    //*** end new
    hVl_key = pref.StorageKey("volumetric enthalpy liquid");
    hVv_key = pref.StorageKey("volumetric enthalpy vapor");
    rl_transport_key  = pref.StorageKey("density liquid transport");
    rv_transport_key  = pref.StorageKey("density vapor transport");
    nQ_key  = pref.StorageKey("nodal fluid volume source");
    vol_fac_key = pref.StorageKey("volume factor");
    state_key = pref.StorageKey("fluid state");
    state_p_key = pref.StorageKey("previous fluid state");
    src_h_key = pref.StorageKey("fluid source h");
    src_rate_key  = pref.StorageKey("fluid source rate");
    CT_key  = pref.StorageKey("nodal total compressibility");
    betar_key = pref.StorageKey("nodal compressibility rock");
    mml_key = pref.StorageKey("liquid mass mobility");
    mmv_key = pref.StorageKey("vapor mass mobility");
    mmld_key  = pref.StorageKey("liquid mass mobility density");
    mmvd_key  = pref.StorageKey("vapor mass mobility density");
    eml_key = pref.StorageKey("liquid enthalpy mobility");
    emv_key = pref.StorageKey("vapor enthalpy mobility");
    emld_key  = pref.StorageKey("liquid enthalpy mobility density");
    emvd_key  = pref.StorageKey("vapor enthalpy mobility density");
    rvl_key = pref.StorageKey("relperm viscosity liquid");
    rvv_key = pref.StorageKey("relperm viscosity vapor");
    bfm_key = pref.StorageKey("boundary flow mass");
    bfe_key = pref.StorageKey("boundary flow enthalpy");
    bfs_key = pref.StorageKey("boundary flow salt");
    // time factor based on fluid props
    time_factor_key = pref.StorageKey("time factor");
    ref_h_top_key = pref.StorageKey("reference specific enthalpy top");
    // salt stuff
    wt_key  = pref.StorageKey("salinity");
    sh_key  = pref.StorageKey("saturation halite");
    rh_key  = pref.StorageKey("density halite");
    mh_key  = pref.StorageKey("solid mass halite");
    hh_key  = pref.StorageKey("enthalpy halite");
    hVh_key = pref.StorageKey("volumetric enthalpy halite");
    hCh_key = pref.StorageKey("enthalpy content halite");
    xVl_key = pref.StorageKey("volumetric salinity liquid");
    xCl_key = pref.StorageKey("salt content liquid");
    xClp_key  = pref.StorageKey("previous salt content liquid");
    xVv_key = pref.StorageKey("volumetric salinity vapor");
    xCv_key = pref.StorageKey("salt content vapor");
    xCvp_key  = pref.StorageKey("previous salt content vapor");
    xVh_key = pref.StorageKey("volumetric salinity halite");
    xCh_key = pref.StorageKey("salt content halite");
    xCf_key = pref.StorageKey("salt content fluid");
    xf_key  = pref.StorageKey("fluid salt fraction");
    xl_key  = pref.StorageKey("salt fraction liquid");
    xv_key  = pref.StorageKey("salt fraction vapor");
    xh_key  = pref.StorageKey("salt fraction halite");
    xml_key = pref.StorageKey("liquid salt mobility");
    xmv_key = pref.StorageKey("vapor salt mobility");
    src_wt_key  = pref.StorageKey("fluid source wt");
    msp_key = pref.StorageKey("previous mass salt");

    this->ApplicationLevel(REGION);
    this->ApplicationTarget(NODE);

    Liquid.InitToZero();
    Vapor.InitToZero();
    Bulk.InitToZero();
    Salt.InitToZero();
  }



  template<size_t dim>
  NaClH2OPropertiesVisitorPHX<dim>::~NaClH2OPropertiesVisitorPHX()
  {
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::Visit(Region<dim>* n)
  {

  }

  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::Visit( Node<dim>* n )
  {

    double64 id(n->Idx());

    // ****************************************
    // 1. Read all nodal varriables of interest
    // ****************************************
    ReadAllVariables( n );

    old_state = int(state()+0.01); // should always give the correct result
    state_p() = double64(old_state);

    // ****************************************
    // 2. Initialize some bools
    // ****************************************
    bogus_variables   = false;
    pure_halite       = false;
    CheckForOutOfRange( n );
    CheckBoundaryFlags( n );

    // ****************************************
    // 3. Prepare a few variables for thermal equilibration code
    // ****************************************
    CalculateAbsoluteVariables( );
    UpdateSowatVariables();
    t_diffusion_       = t();

    // ****************************************
    // 4. EQUILIBRATE, core task of the visitor
    // ****************************************
    Equilibrate( n );
    

    // ****************************************
    // 4a. Functions for "compressibility bug"
    // ****************************************
    if(adjust_compressibility_after_phasechange)
      {
        CheckPhaseChange();
        CheckVolumeMismatchCompensation();
      }
    
    
    // ****************************************
    // 5. Prepare variables to be passed back to CSMP
    // ****************************************
    UpdateCSMPVariables( *n );


    // ****************************************
    // 5a. time stepping influence
    // ****************************************
    TimeStepAdjustment();


    // ****************************************
    // 6. Volume Mismatch
    // ****************************************
    VolumeFactorComputations( n );


    PrepareVariablesForStorage();

    // if (essentiallyEqual(sh(),1.0,numeric_limits<double64>::epsilon()))
    //   {
    // 	csmp_error.notice( FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
    // 			"\nPure Halite!");
    // 	int sto;
    // 	cout << "pure halite";
    // 	cin >> sto;
    //   }
    
    if(Bulk.state!=none)
      {
        StorePropertiesAndFlags( *n );
      }
    else 
      {
        cerr << "\nNode: ("<<n->x()<<", "<<n->y()<<", Pressure: "<<p()<<" Pa, total enthalpy: "<<H_current_<<" J, total mass: " << mt() <<" kg, rock temperature: "<<t()<<" oC " << endl;
        ScreenOutputSowatVariables();
        csmp_error.notice( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                        "\nFluid's state is undefined, not storing result from equilibration!");
      }
	
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::Equilibrate( Node<dim>* n )
  {
    if (open_boundaries && t.Flag() == DIRICH && n->AtBoundary() != NOT && n->AtBoundary() != INTERNAL)
      {
        BoundaryIteration();
      }
    else
      {
        equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
      }
    
    if(equilibrator.Fatal())
      {
        cout << "NaClH2OPropertiesVisitorPHX<dim>::Equilibrate received message Fatal() from equilibrator ...\n";
        cout << "for conditions : \n";
        ScreenOutputSowatVariables();
        csmp_error.notice( FATAL_ERROR, 
                        "NaClH2OPropertiesVisitorPHX<dim>::Visit( Node<dim>* n ) -",
                        "received message Fatal() from equilibrator ... teminating!!!");
      }
    current_state = Bulk.state; // that should be type safe
    state()       = double64(current_state);    
  }



  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::SetTimeIncrement( double64 time_increment )
  {
    dt_ = time_increment;
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>:: TopBoundaryHandling( bool TB )
  {
    top_boundary = TB;
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::ReadAllVariables( Node<dim>* n )
  {
    n->Read( t_key,  t );
    n->Read( p_key,  p );
    n->Read( mt_key, mt  );
    n->Read( hCl_key, hCl );
    n->Read( hCv_key, hCv );
    n->Read( tp_key, tp );
    n->Read( phi_key,   phi );
    n->Read( rr_key,  rr );
    n->Read( cpr_key, cpr );
    n->Read( hClp_key,hClp );
    n->Read( hCvp_key,hCvp );   
    n->Read( cpv_key, cpv );   
    n->Read( src_h_key, src_h );
    n->Read( src_rate_key, src_rate );
    n->Read( ml_key, ml  );
    n->Read( mv_key, mv  );
    n->Read( mlp_key,mlp );
    n->Read( mvp_key,mvp );
    n->Read( mtp_key,mtp );
    n->Read( Htp_key, Htp );
    n->Read( sl_key, slp );
    n->Read( sv_key, svp );
    n->Read( betar_key, beta_rock );
    //*** new TD May 2011
    n->Read( beta_p_key, beta_p );
    n->Read( beta_ref_key, beta_ref );
    n->Read( state_key, state ); //*** probably not new
    n->Read( apc_key, after_phasechange_counter);
    n->Read( dpc_key, dangerous_phase_change);
    //*** end new

    // salt
    n->Read( xCl_key, xCl );
    n->Read( xClp_key,xClp );
    n->Read( xCv_key, xCv );
    n->Read( xCvp_key,xCvp );
    n->Read( xCf_key, xCf );
    //	n->Read( xCfp_key,xCfp );
    n->Read( mh_key, mh  );
    n->Read( sh_key, sh );  //??? GG
    //	n->Read( sh_key, shp );
    n->Read( msp_key,msp );
    n->Read( src_wt_key, src_wt );
    n->Read( wt_key, wt );

    // for open boundary calculations
    n->Read( rho_bulk_key, rho_bulk );
    n->Read( mul_key,mulp );
    n->Read( muv_key,muvp );
    n->Read( xl_key,xlp );
    n->Read( xv_key,xvp );
    n->Read( hl_key,hlp );
    n->Read( hv_key,hvp );
    n->Read( hh_key,hhp );
    n->Read( rl_key,rlp );
    n->Read( rv_key,rvp );

    n->Read( time_factor_key, time_factor );
    n->Read( ref_h_top_key, reference_enthalpy_top );

  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::CalculateAbsoluteVariables( )
  {
    pore_volume                 = phi();
    rock_volume                 = 1.0 - phi();
    // energy added to fluid by diffusion
    //dT_diff_       = t() - tp();
    hrock_prev     = rock.Enthalpy(tp());
    hrock_curr     = rock.Enthalpy(t());
    //dh_rock_diff_  = rock_volume * cpr() * rr() * dT_diff_;
    dh_rock_diff_  = rock_volume * rr() * (hrock_curr - hrock_prev);
    // *** This needs very careful re-evaluation !!! ***
    //   dh_fluid_diff_ = pore_volume * cpf() * mt() * dT_diff_;
    dh_fluid_diff_ = 0.;
    dhCl_ = pore_volume * (hCl() - hClp());
    dhCv_ = pore_volume * (hCv() - hCvp());

    // mass, salt, energy added by sources
    if(!first)
      {
        //	    src_rate() /= fv_man.ControlVolume(n.Idx()-1U).CVPoreVolume(); // mt-like property (I hope)
        src_rate() /= cpv();
        // umbenennen pcv

        dmt_src     = src_rate() * dt_ ;//?  * pore_volume; 
        dh_src      = dmt_src * src_h();//?  * pore_volume;
        ds_src      = dmt_src * src_wt();//? * pore_volume;

      }
    else
      {
        dmt_src     = 0.0; 
        dh_src      = 0.0;
        ds_src      = 0.0;
      }	    


    // liquid + vapor mass per pore volume (including dissolved salt)
    dml_ = ml() - mlp();
    dmv_ = mv() - mvp();
    dxCl_    = xCl() - xClp(); // [kg] added by liq per pore volume
    dxCv_    = xCv() - xCvp(); // [kg] added by vap per pore volume
    dx_diff_ = 0.; // no salt diffusion
		
    ms()  = msp();
    ms() += dxCl_;
    ms() += dxCv_;
    ms() += dx_diff_;
    /* //new */	ms() += ds_src*0.01;

    mt() = mtp();
    mt()+= dml_;
    mt()+= dmv_;
    mt()+= dmt_src;

    /*    if(mt()<0.01)
          {
          cout <<"\nNaClH2OPropertiesVisitorPHX: WARNING: Reseting mt() from "<<mt()<<" to 0.01 "<< endl;
          mt() = 0.01;
          }*/

    m_fluid_ = mt() * pore_volume;      // [kg]
    m_rock_  = rr() * rock_volume;      // [kg]

    MS_current_ = ms() * pore_volume;

    // salinity
    wt_ = MS_current_;                 // kg 's of salt
    wt_ /= m_fluid_;                   // mass fraction
    wt_ *= 100.;                       // wt%

    if( (wt_ - 1.0e2) > 0. )
      wt_ = 1.0e2;

    // specific enthalpy of fluid
    h_fluid_  = (hCl() + hCv() + hCh())/ mt();
    h_fluid_  += dh_src / mt();
	
    bfm() = bfe() = bfs() = 0.0;

  }



  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::UpdateSowatVariables()
  {
    t_          = t();
//    p_bar_      = p()/1.e5;
	p_bar_ = p();
	x_ = Weight2XNaCl(wt_);

    if ( (x_ - 1.0e0) > 0. ) x_ = 1.0e0;

    cp_rock_    = rock.HeatCapacity(t_);
    rho_rock_   = rr();
    phi_        = phi();
    tp_         = tp();

// PW 14.12.2015 shouldn't this be in Pa?
    p_current_  = p_bar_;
//	p_current_ = p();

	H_previous_ = Htp();
    H_current_  = Htp();
    H_current_ += dh_fluid_diff_;
    H_current_ += dh_rock_diff_;
    H_current_ += dhCl_;
    H_current_ += dhCv_;
    H_current_ += dh_src*pore_volume; // May 26, 2009, from Philipp's dicovery: *pore_volume is new
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables( Node<dim>& n )
  {
    // Fluid properties
    t()   = t_ = Bulk.t;
    cpr() = cp_rock_ = rock.HeatCapacity(t_);
    ncp() = cpr()*rr()*(1.-phi());

    sl()  = Liquid.s;
    sv()  = Vapor.s;
    sh()  = Salt.s;

    // COARSE FIXES
    if( (abs(sh()-1.0e0)) < 1.0e-10 ) // otherwise NANs in properties from SOWAT
      {
        cout<<"\nPure halite forming at node "<<n.Idx()<<", reseting saturations from: state = "<<Bulk.state<<", sl = "<<sl()<<", sv = "<<sv()<<", sh = "<<sh()<< endl;

        Vapor.InitToZero();
        Liquid.InitToZero();
        Salt.state = H;
        Salt.s = 1.0;
        Salt.mf = 1.0;
        Bulk = Salt;
        sl() = Liquid.s;
        sv() = Vapor.s;
        sh() = Salt.s;

        pure_halite = true;
        ScreenOutputSowatVariables();
        cout<<" to: sl = "<<sl()<<", sv = "<<sv()<<", sh = "<<sh()<< endl;
      }

    if( t() > 799.99 && x_ > 0.999 && Bulk.state == V ) // pure halite in VH field
      {
        cout<<"\nPure halite forming at node "<<n.Idx()<<", in VAPOR field, reseting state to HALITE "<< endl;

        Vapor.InitToZero();
        Liquid.InitToZero();
        Salt.state = H;
        Salt.s = 1.0;
        Salt.mf = 1.0;
        Bulk = Salt;
        sl() = Liquid.s;
        sv() = Vapor.s;
        sh() = Salt.s;

        pure_halite = true;

        cout<<" to: sl = "<<sl()<<", sv = "<<sv()<<", sh = "<<sh()<<" and state = "<<Bulk.state<< endl;
      }

    rl()  = Liquid.rho;
    rv()  = Vapor.rho;
    rh()  = Salt.rho;

    mul() = Liquid.mu;
    muv() = Vapor.mu;

    hf()  = Bulk.h;
    hl()  = Liquid.h;
    hv()  = Vapor.h;
    hh()  = Salt.h;

    cpf() = Bulk.cp;

    wt()  = Bulk.wt;
    xf()  = Bulk.smf;
    xl()  = Liquid.smf;
    xv()  = Vapor.smf;
    xh()  = Salt.smf;

    beta()  = Bulk.beta;

    /*    cout << Bulk.state << endl;
          if (old_state == V && Bulk.state == F)
          {
          cout << "beta(): " << beta() << endl;
          cout << "Vapor.beta: " << Vapor.beta << endl;
          cout << "Liquid.beta: " << Liquid.beta << endl;
          cout << endl;
          }*/
      
    // Transport properties
    ml()  = rl()  * sl();
    mv()  = rv()  * sv();
    hVl() = hl()  * rl();
    hCl() = hVl() * sl();
    hVv() = hv()  * rv();
    hCv() = hVv() * sv();
    hVh() = hh()  * rh();
    hCh() = hVh() * sh();
    xVl() = xl()  * rl();
    xCl() = xVl() * sl();
    xVv() = xv()  * rv();
    xCv() = xVv() * sv();
    xVh() = xh()  * rh();
    xCh() = xVh() * sh();
    xCf() = Bulk.rho * xf();

    //TD: put this into a function and enclose with #ifdef DEBUG
    /*	if( mul()<1.e-10 && muv()<1.0e-10 ){
        bogus_variables = true;
        int cont;
        cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
        cout << "mul, muv = " << mul << "\t" << muv << endl;
        ScreenOutputSowatVariables();
        cin >> cont;
        }
	    
        if(sl()==0.0 && sv()==0.0){
        bogus_variables = true;
        int cont;
        cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
        cout << "sl, sv, sh = " << sl << "\t" << sv << endl;
        ScreenOutputSowatVariables();
        cin >> cont;
        }
	
        if(sl()<0.0 || sv()<0.0){
        bogus_variables = true;
        int cont;
        cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
        cout << "sl, sv, sh = " << sl << "\t" << sv << endl;
        ScreenOutputSowatVariables();
        cin >> cont;
        }
	    
        if(ml()==0.0 && mv()==0.0){
        bogus_variables = true;
        int cont;
        cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
        cout << "ml, mv = " << ml << "\t" << mv << endl;
        ScreenOutputSowatVariables();
        cin >> cont;
        }

        if(beta() < 0.)
        {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
        "\nWARNING: NEGATIVE BETA.");
        cout <<"\n ("<<n.x()<<", "<<n.y()<<"): beta() = "<<beta()<< endl;
        if(Bulk.state == L || Bulk.state == F)  beta() = Liquid.beta;
        if(Bulk.state == V)  beta() = Vapor.beta;
        if(Bulk.state == VL) beta() = TwoPhasePureWaterCompressibility(Liquid.cp, Vapor.cp);
        cout <<", reseting to "<<beta()<< endl;
        ScreenOutputSowatVariables();
        }

        if(muv() < 0.)
        {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
        "\nWARNING: NEGATIVE VISCOSITY VAPOR.");
        cout <<"\n ("<<n.x()<<", "<<n.y()<<"): muv() = "<<muv()<<" reseting to 1.0e-5"<< endl;
        muv() = 1.0e-5;
        }
    */
    //*** April 29, 2009 : needs attention forr pure halite (state == 8 ) case

    if(Bulk.state != 8)
      {
        if( mul()<0.0 )
          mul() = fabs(mul());
        if( muv()<0.0 )
          muv() = fabs(muv());
        if( mul()<1.e-10 && muv()<1.0e-10 )
          {

            bogus_variables = true;
            //	    int cont;
            cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
            cout << "mul, muv = " << mul << "\t" << muv << endl;
            ScreenOutputSowatVariables();
            //	    cin >> cont;
          }

        if(sl()==0.0 && sv()==0.0 && sh()==0.0)
          {
            bogus_variables = true;
            int cont;
            cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
            cout << "sl, sv, sh = " << sl << "\t" << sv << "\t" << sh << endl;
            ScreenOutputSowatVariables();
            //	    cin >> cont;
          }
	    
        if(sl()<0.0 || sv()<0.0 || sh()<0.0)
          {
            bogus_variables = true;
            int cont;
            cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
            cout << "sl, sv, sh = " << sl << "\t" << sv << "\t" << sh << endl;
            ScreenOutputSowatVariables();
            //	    cin >> cont;
          }
	    
        if(ml()==0.0 && mv()==0.0 && mh()==0.0)
          {
            bogus_variables = true;
            int cont;
            cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:\n";
            cout << "ml, mv, mh = " << ml << "\t" << mv << "\t" << mh << endl;
            ScreenOutputSowatVariables();
            //	    cin >> cont;
          }
        // 	if(bogus_variables){
        // 	    cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::UpdateCSMPVariables: Warning: Bogus variables coming from SOWAT!:"<< endl;
      }
    if(beta() < 0.)
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX_TOPBC<dim>::Visit",
                        "\nWARNING: NEGATIVE BETA.");
        cout <<"\n ("<<n.x()<<", "<<n.y()<<"): beta() = "<<beta()<< endl;
        if(Bulk.state == L || Bulk.state == F || Bulk.state == LH)  beta() = Liquid.beta;
        if(Bulk.state == V || Bulk.state == VH)  beta() = Vapor.beta;
        if(Bulk.state == VL|| Bulk.state == VLH) beta() = TwoPhasePureWaterCompressibility(Liquid.cp, Vapor.cp);
        if(Bulk.state == H) beta() = Salt.beta;
        if(Bulk.state == 9) beta() = Salt.beta;

        cout <<", reseting to "<<beta()<< endl;
        ScreenOutputSowatVariables();
      }
    if(muv() < 0.)
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX_TOPBC<dim>::Visit",
                        "\nWARNING: NEGATIVE VISCOSITY VAPOR.");
        cout <<"\n ("<<n.x()<<", "<<n.y()<<"): muv() = "<<muv()<<" reseting to 1.0e-5"<< endl;
        muv() = 1.0e-5;
      }
	
  }



  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::BoundaryHandling( Node<dim>* n )
  {
    h_fluid_ = Bulk.h;
    Bulk   = equilibrator.ReportBulkProperties(t_,p_bar_,x_,h_fluid_);
    Liquid = equilibrator.ReportLiquidProperties(t_,p_bar_,x_,h_fluid_);
    Vapor  = equilibrator.ReportVaporProperties(t_,p_bar_,x_,h_fluid_);
    Salt   = equilibrator.ReportSaltProperties(t_,p_bar_,x_,h_fluid_);
    UpdateCSMPVariables( *n );

    mt()  = mtp()  = Bulk.rho;
    //    H_current_     = Bulk.h * m_fluid_ + rr()*cpr()*rock_volume*t();
    H_current_     = Bulk.h * m_fluid_ + rr()*rock_volume*rock.Enthalpy(t_);
    Htp()          = H_current_;
    rl_transport() = rl();
    rv_transport() = rv();
    rho_bulk()     = Bulk.rho;
    nQ()           = 0.;
    vol_fac()      = 1.0;
  }


  template<size_t dim>
  double64 NaClH2OPropertiesVisitorPHX<dim>::TwoPhasePureWaterCompressibility(double64 cpl, double64 cpv)
  {
    // compressibility of liquid/vapor mixture after Grant and Sory, WRR 15(3) p. 684-686
    // compressibility modified as to add the vapor contribution cpv * rv * sv * phi as well
    double64 product, b; 
    product = ( rl() - rv() ) / ( ( hv() - hl() ) * rl() * rv() );
    b  = ( 1.0 - phi() ) * cpr() * rr();
    b += phi() * rl() * cpl * sl();
    b += phi() * rv() * cpv * sv();
    b *= product * product;
    b *= ( t() + kelvin ) * 1.0 / phi();

    return b;
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::StorePropertiesAndFlags( Node<dim>& n )
  {


    n.Store( mt_key, mt );
    n.Store( mtp_key,mtp );
    n.Store( Htp_key, Htp );
    n.Store( msp_key, msp );

    //   if (n.Idx() == 2 && sh() > 0.)
    //       cout << endl;

    if( bogus_variables )
      {
        cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::StorePropertiesAndFlags: Warning: Bogus variables coming from SOWAT!"<< endl;
      }
    else
      {
        // temperature
        n.Store( t_key,  t );
        n.Store( sl_key, sl );
        n.Store( sv_key, sv );
        n.Store( sh_key, sh );
        //*** new TD May 2011
        n.Store( state_key, state ); // really new?
        n.Store( apc_key, after_phasechange_counter );
        n.Store( dpc_key, dangerous_phase_change );
        //*** end new

        // densities
        n.Store( rho_bulk_key, rho_bulk );
        n.Store( rl_key, rl );
        n.Store( rv_key, rv );
        n.Store( rh_key, rh );
        n.Store( mf_key, mf );
        n.Store( ml_key, ml );
        n.Store( mv_key, mv );
        n.Store( mh_key, mh );

        //viscosities
        n.Store( mul_key,mul );
        n.Store( muv_key,muv );

        //enthalpy variables
        n.Store( hf_key, hf );
        n.Store( hl_key, hl );
        n.Store( hVl_key,hVl );
        n.Store( hCl_key,hCl );
        n.Store( hv_key, hv );
        n.Store( hVv_key,hVv );
        n.Store( hCv_key,hCv );
        n.Store( hh_key, hh );
        n.Store( hVh_key,hVh );
        n.Store( hCh_key,hCh );
        n.Store( cpf_key,cpf );
        // JPW Nov 2010
        n.Store( cpr_key,cpr );
        n.Store( ncp_key,ncp );

        // expansivities
        n.Store( beta_key,beta );
        n.Store( beta_p_key,beta_p );
        n.Store( beta_ref_key,beta_ref );
        n.Store( CT_key,  CT );

        // extras
        n.Store( rl_transport_key, rl_transport );
        n.Store( rv_transport_key, rv_transport );
        n.Store( vol_fac_key, vol_fac );
        n.Store( nQ_key,    nQ );

        n.Store( phi_key,   phi );
	    
        n.Store( mml_key,   mml );
        n.Store( mmv_key,   mmv );
        n.Store( mmld_key,   mmld );
        n.Store( mmvd_key,   mmvd );
        n.Store( eml_key,   eml );
        n.Store( emv_key,   emv );
        n.Store( emld_key,   emld );
        n.Store( emvd_key,   emvd );
        n.Store( xml_key,   xml );
        n.Store( xmv_key,   xmv );
        n.Store( rvl_key,   rvl );
        n.Store( rvv_key,   rvv );
        n.Store( bfm_key,   bfm );
        n.Store( bfe_key,   bfe );
        n.Store( bfs_key,   bfs );

        // salinity variables
        //    if(wt.Flag()!=DIRICH)
        //     {
        // n.Store( wt_key, wt );
        n.Store( xf_key, xf );
        n.Store( xl_key, xl );
        n.Store( xVl_key,xVl );
        n.Store( xCl_key,xCl );
        n.Store( xv_key, xv );
        n.Store( xVv_key,xVv );
        n.Store( xCv_key,xCv );
        n.Store( xh_key, xh );
        n.Store( xVh_key,xVh );
        n.Store( xCh_key,xCh );
        n.Store( xCf_key,xCf );

        n.Store( time_factor_key, time_factor );

      }

  }




  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::StoreInitialPropertiesAndFlags( Node<dim>& n )
  {
    int temp;

    // temperature
    n.Store( t_key,  t );
    n.Store( sl_key, sl );
    n.Store( sv_key, sv );
    n.Store( sh_key, sh );
    //*** new TD May 2011
    n.Store( state_key, state ); // really new?
    n.Store( state_p_key, state_p );
    n.Store( apc_key, after_phasechange_counter );
    n.Store( dpc_key, dangerous_phase_change );
    // new june 2011
    n.Store( cpr_key,cpr );
    n.Store( ncp_key,ncp );

    //*** end new

    // densities
    n.Store( mt_key, mt );
    n.Store( mtp_key,mtp );
    n.Store( rho_bulk_key, rho_bulk );
    n.Store( rl_key, rl );
    n.Store( rv_key, rv );
    n.Store( rh_key, rh );
    n.Store( mf_key, mf );
    n.Store( ml_key, ml );
    n.Store( mv_key, mv );
    n.Store( mh_key, mh );

    //viscostities
    n.Store( mul_key,mul );
    n.Store( muv_key,muv );

    //enthalpy variables
    n.Store( hf_key, hf );
    n.Store( hl_key, hl );
    n.Store( hVl_key,hVl );
    n.Store( hCl_key,hCl );
    n.Store( hv_key, hv );
    n.Store( hVv_key,hVv );
    n.Store( hCv_key,hCv );
    n.Store( hh_key, hh );
    n.Store( hVh_key,hVh );
    n.Store( hCh_key,hCh );
    n.Store( cpf_key,cpf );

    // total salt mass
    n.Store( msp_key, msp );

    // salinity variables
    n.Store( wt_key, wt );
    n.Store( xf_key, xf );
    n.Store( xl_key, xl );
    n.Store( xVl_key,xVl );
    n.Store( xCl_key,xCl );
    n.Store( xv_key, xv );
    n.Store( xVv_key,xVv );
    n.Store( xCv_key,xCv );
    n.Store( xh_key, xh );
    n.Store( xVh_key,xVh );
    n.Store( xCh_key,xCh );
    n.Store( xCf_key,xCf );
	
    // total enthalpy
    n.Store( Htp_key, Htp );

    // expansivities
    n.Store( beta_key,beta );
    n.Store( beta_p_key,beta_p );
    n.Store( beta_ref_key,beta_ref );
    n.Store( CT_key,CT );

    // extras
    n.Store( rl_transport_key, rl_transport );
    n.Store( rv_transport_key, rv_transport );
    n.Store( vol_fac_key, vol_fac );
    n.Store( nQ_key,    nQ );

    n.Store( phi_key,   phi );


    // new JPW
    n.Store( mml_key,   mml );
    n.Store( mmv_key,   mmv );
    n.Store( mmld_key,   mmld );
    n.Store( mmvd_key,   mmvd );
    n.Store( eml_key,   eml );
    n.Store( emv_key,   emv );
    n.Store( emld_key,   emld );
    n.Store( emvd_key,   emvd );
    n.Store( xml_key,   xml );
    n.Store( xmv_key,   xmv );
    n.Store( rvl_key,   rvl );
    n.Store( rvv_key,   rvv );
    n.Store( bfm_key,   bfm );
    n.Store( bfe_key,   bfe );
    n.Store( bfs_key,   bfs );

    // all transported variables fluid mass, density, and volumetric enthalpy 
    // have the same flag as temperature, salt transport properties have the same flag as salinity

    if (!open_boundaries)
      {
        n.Status( ml_key,  t.Flag() );
        n.Status( mv_key,  t.Flag() );
        n.Status( mh_key,  t.Flag() );
        n.Status( rl_key , t.Flag() );
        n.Status( rv_key , t.Flag() );
        n.Status( rh_key , t.Flag() );
        n.Status( rl_transport_key , t.Flag() );
        n.Status( rv_transport_key , t.Flag() );

        n.Status( hVl_key, t.Flag() );
        n.Status( hVv_key, t.Flag() );
        n.Status( hVh_key, t.Flag() );
        n.Status( hCl_key, t.Flag() );
        n.Status( hCv_key, t.Flag() );
        n.Status( hCh_key, t.Flag() );

        n.Status( xVl_key, wt.Flag() );
        n.Status( xVv_key, wt.Flag() );
        n.Status( xVh_key, wt.Flag() );
        n.Status( xCl_key, wt.Flag() );
        n.Status( xCv_key, wt.Flag() );
        n.Status( xCh_key, wt.Flag() );
        n.Status( xCf_key, wt.Flag() );

      } // end !open_boundaries

    if ( open_boundaries &&
         !const_top_pressure &&
         t.Flag() == DIRICH &&
         n.AtBoundary() != NOT &&
         n.AtBoundary() != INTERNAL )
      n.Store( ref_h_top_key, reference_enthalpy_top );

  }

  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX()
  {
    // Only use for calculating properties when calculating static pressure

    H2ONaClFluidProperties   fluid( t_, p_bar_, x_, h_fluid_, cp_rock_, rho_rock_, phi_, false);
//	  H2ONaClFluidProperties   fluid(t_, p_current_, x_, h_fluid_, cp_rock_, rho_rock_, phi_, false);
	//    fluid(tdummy, p_current_eq, x_current_eq, hdummy, cp_rock, rho_rock, phi, verbose),
	//equilibrator(m_rock_, cp_rock_, rho_rock_, phi_, m_fluid_, wt_, tp_, p_current_, H_current_, H_previous_, fixed_temperature, t_fixed, t_diffusion_, verbose_eq),

    cout <<"\nNaClH2OPropertiesVisitorPHX<dim>::InitialPropertiesFromPTX(): Initialising fluid properties"<< endl;

    typename std::deque<Node<dim> >::iterator it;

    int count(0);
    for ( it = pmesh.NodesBegin(); it != pmesh.NodesEnd(); it++ )
      {
        //	cout << count << endl;
        count++;
        // 1. reading input variables (fluid)
        ReadAllVariables( &(*it) );
        tp() = t();
        mtp() = 0.01;
        ml() = mlp() = mv() = mvp() = src_rate() = 0.0;
        cpr() = rock.HeatCapacity(t()); 
        ncp() = cpr()*rr()*(1.-phi());

        // Fixed temperature or not
        if(t.Flag()==DIRICH)
          {
            fixed_temperature = true;
            t_fixed = t();
          }
        else
          {
            fixed_temperature = false;
          }

        //	cout << "a\n";
        pore_volume = phi();
        rock_volume = 1.0 - phi();
        //	h_fluid_ = 2.086e6;

        // warning if variables below minimum of lookup table
        if ( p() < 101325.0 || t() < 5.0 )
          {
            cerr << "\nNode: " << it->Idx() << ", Pressure: " << p() << " Pa, T: " << t() << " oC" << endl;
            csmp_error.notice( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::CalculateInitialPropertiesFromPT",
                            "\nPressure or temperature below minimum values of lookup table, erroneous results are possible...!");
          }
        // exit if negative values are encountered
        if ( p() < 0.0 || t() < 0.0 )
          {
            cerr<<"\nNode: "<<it->Idx()<<", Pressure: "<< p() <<" Pa, T: "<< t() <<" oC"<<endl;
            csmp_error.notice( FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::CalculateInitialPropertiesFromPT",
                            "\nNegative input variable, terminating...!");
          }
        if ( p() > 5000.0e5 || t() > 1000.0 )
          {
            cerr<<"\nNode: " << it->Idx() << ", Pressure: " << p() << " Pa, T: " << t() <<" oC" <<endl;
            csmp_error.notice( FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::CalculateInitialPropertiesFromPT",
                            "\nTooLarge input variable, terminating...!");
          }

        // Absolute variables
        //	cout << "b\n";
        CalculateAbsoluteVariables( );
        wt_      = wt();
        h_fluid_ = 2086000.; // critical enthalpy of water
		UpdateSowatVariables();
        //	cout << "c\n";
        //ScreenOutputSowatVariables();

		// This is code from main_equilibrator_test.cpp
        Bulk       = fluid.BulkProperties();
        //	cout << "d\n";
        Liquid     = fluid.ReportLiquidProperties();
        //	cout << "e\n";
        Vapor      = fluid.ReportVaporProperties();
        //	cout << "f\n";
        Salt       = fluid.ReportSaltProperties();
        //	cout << "g\n";
        h_fluid_   = Bulk.h;
	
        // // Bulk     = equilibrator.Bulk();
        // // Liquid   = equilibrator.Liquid();
        // // Vapor    = equilibrator.Vapor();
        // // Salt     = equilibrator.Salt();

        // Bulk   = equilibrator.ReportBulkProperties(t_,p_bar_,x_,h_fluid_);
        // Liquid = equilibrator.ReportLiquidProperties(t_,p_bar_,x_,h_fluid_);
        // Vapor  = equilibrator.ReportVaporProperties(t_,p_bar_,x_,h_fluid_);
        // Salt   = equilibrator.ReportSaltProperties(t_,p_bar_,x_,h_fluid_);
        //	cout << "h\n";
        UpdateCSMPVariables( *it );
        beta_p() = beta();
        //*** new TD May 2011
        after_phasechange_counter() = 0.0;
        dangerous_phase_change() = 0.0;
        //***end new
        //	cout << "i\n";
        mt()   = mtp() = Bulk.rho;
        m_fluid_       = mt() * pore_volume;

        H_current_     = Bulk.h * m_fluid_ + rr()*rock.Enthalpy(t())*rock_volume;
        Htp()          = H_current_;
        ms()  = msp() = xCl() + xCv() + xCh();
        mh()  = Salt.mf * mt();

        rl_transport() = rl();
        rv_transport() = rv();
        rho_bulk()     = Bulk.rho;
        nQ()           = 0.;
        vol_fac()      = 1.0;

        // nodal total compressibility
        CT()  = beta()*pore_volume;
        CT() += beta_rock()*rock_volume;
        CT() *= mt();

        // new JPW
        // calculating nodal variables
        // liquid/vapor mass/enthalpy mobilities
    
        mml() = mmv() = mmld() = mmvd() = 0.0;
        eml() = emv() = rvl() = rvv() = 0.0;
        emld() = emvd() = 0.0;
        bfm() = bfe() = bfs() = xml() = xmv() = 0.0;
        //	krl = EffectiveLiquidSaturationHalitePresent( sl(), sv() );
        //	krv = 1.0 - krl;
        krl = RelativePermeabilityLiquid( sl(), sv() );
        krv = 1.0 - sh() - krl;
        //	cout << "j\n";
        if (sl() > 0.0) mml() = krl * rl_transport() / mul();
        if (sv() > 0.0) mmv() = krv * rv_transport() / muv();
        //	if (sl() > 0.0) mmld() = mml() * rl();
        //	if (sv() > 0.0) mmvd() = mmv() * rv();
        if (sl() > 0.0) mmld() = mml() * rl_transport();
        if (sv() > 0.0) mmvd() = mmv() * rv_transport();
        if (sl() > 0.0) eml() = krl * hVl() / mul();
        if (sv() > 0.0) emv() = krv * hVv() / muv();
        //	if (sl() > 0.0) emld() = eml() * rl();
        //	if (sv() > 0.0) emvd() = emv() * rv();
        if (sl() > 0.0) emld() = eml() * rl_transport();
        if (sv() > 0.0) emvd() = emv() * rv_transport();
        if (sl() > 0.0) xml() = krl * xVl() / mul();
        if (sv() > 0.0) xmv() = krv * xVv() / muv();
        if (sl() > 0.0) rvl() = krl / mul();
        if (sv() > 0.0) rvv() = krv / muv();

        reference_enthalpy_top() = h_fluid_;
    
        //*** new TD May 2011
        state()   = Bulk.state;
        state_p() = state();
        //*** end new
        //   cout << "k\n";
        StoreInitialPropertiesAndFlags( *it );
        //   cout << "l\n";
      }
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::ScreenOutputSowatVariables()
  {
    cout.precision(15);

    cout <<"\nSOWAT parameters:"<< endl;
    cout <<"\nBulk.state = "<<Bulk.state<<", previous_state = "<<previous_state<< endl;

    cout <<"double64 rho_rock   = "<<rr()<<";\n";
    cout <<"double64 phi        = "<<phi()<<";\n";
    cout <<"double64 h_fluid    = "<<h_fluid_<<"; // this is one passed to FRE"<<endl;
    cout <<"double64 m_fluid    = "<<m_fluid_<<";\n";
    cout <<"double64 hf         = "<<hf()<<"; // this is one coming back from FRE"<<endl;
    cout <<"double64 m_rock     = "<<m_rock_<<";\n";
    cout <<"double64 t          = "<<t_<<";\n";
    cout <<"double64 p          = "<<p_bar_<<";\n";
    cout <<"double64 x          = "<<x_<<";\n";
    cout <<"double64 cp_rock    = "<<cp_rock_<<";\n";
    cout <<"double64 wt         = "<<wt_<<";\n";
    cout <<"double64 t_previous = "<<tp_<<";\n";
    cout <<"double64 p_current  = "<<p_current_<<";\n";
    cout <<"double64 H_current  = "<<H_current_<<";\n";
    cout <<"double64 H_previous = "<<H_previous_<< ";\n"; 
    cout <<"\n    x - 1.0    = "<<(x_-1.0)<<endl;
    cout <<"\n    wt - 100.  = "<<(wt_-100.0)<<endl;

    cout.precision(8);

    cout <<"Fluid variables:"<< endl;
    cout <<"rl() = "<<rl()<<", rv() = "<<rv()<<", rh() = "<<rh()<< endl;
    cout <<"sl() = "<<sl()<<", sv() = "<<sv()<<", sh() = "<<sh()<< endl;
    cout <<"hl() = "<<hl()<<", hv() = "<<hv()<<", hh() = "<<hh()<< endl;
    cout <<"mul() = "<<mul()<<", muv() = "<<muv()<< endl;
    cout <<"rho_bulk() = "<<rho_bulk()<<", Bulk.rho = "<<Bulk.rho<< endl;
    cout <<"volume factor_LHS = "<<volume_factor_LHS<<", volume factor_RHS = "<<volume_factor_RHS<< endl;
    cout <<"mt() = "<<mt()<<", mtp() = "<<mtp()<<", mf() = "<<mf()<< endl;
    cout <<"dml_ = "<<dml_<<", dmv_ = "<<dmv_<<", dhCl_ = "<<dhCl_<<", dhCv_ = "<<dhCv_<<endl;
    cout <<"ml() = "<<ml()<<", mv() = "<<mv()<<", mh() = "<<mh()<<", mlv_ = "<<mlv_<< endl;
    cout <<"hCl() = "<<hCl()<<", hCv() = "<<hCv()<<", hCh() = "<<hCh()<< endl;
    cout <<"hVl() = "<<hVl()<<", hVv() = "<<hVv()<<", hVh() = "<<hVh()<< endl;
    cout <<"t() = "<<t()<<", p() = "<<p()<<", tp() = "<<tp()<< endl;
    cout <<"cpf() = "<<cpf()<<", cpr() = "<<cpr()<< endl;
    cout <<", beta() = "<<beta()<< endl;
    cout <<"nQ() = "<<nQ()<<", Htp() = "<<Htp()<< endl;

  }

  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::Output_PTXState()
  {
    cout.precision(15);

    cout <<"\nBulk.state = "<<Bulk.state<<", previous_state = "<<previous_state<< endl;
    cout <<"    t          = "<<t_<<endl;
    cout <<"    p          = "<<p_bar_<<endl;
    cout <<"    x          = "<<x_<<endl;
    cout <<"    wt         = "<<wt_<<endl;

    cout.precision(8);

  }



  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::CheckForOutOfRange( Node<dim>* n )
  {
    // Out of range checks:
    if(hCl()<0.) 
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                        "\nNegative hCl:");
        cout <<hCl()<< endl;
      }
    if(hCv()<0.) 
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                        "\nNegative hCv:");
        cout <<hCv()<< endl;
      }
    if(xCl()<0.) 
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX_TOPBC<dim>::Visit",
                        "\nNegative xCl:");
        cout <<xCl()<< endl;
      }
    if(xCv()<0.) 
      {
        csmp_error.notice( WARNING, "NaClH2OPropertiesVisitorPHX_TOPBC<dim>::Visit",
                        "\nNegative xCv:");
        cout <<xCv()<< endl;
      }
    if ( p() < 101325.0 || t() < 5.0 )
      {
        cerr << "\nNode: ("<<n->x()<<", "<<n->y()<<", Pressure: "<<p()<<" Pa, total enthalpy: "<<H_current_<<" J, total mass: " <<
          mt() <<" kg, rock temperature: "<<t()<<" oC " << endl;
        csmp_error.notice( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                        "\nPressure or temperature below minimum values of lookup table, erroneous results are possible...!");
      }
    if ( mt() < 0. )
      {
        cerr << "\nNode: ("<<n->x()<<", "<<n->y()<<", mt: "<<mt()<<", total enthalpy: "<<H_current_<<" J, pressure: " <<
          p() <<" Pa, rock temperature: "<<t()<<" oC " << endl;
        cout <<"\ndml_ = "<<dml_<<", dmv_ = "<<dmv_<<", mtp() = "<<mtp()<<", mt() = "<<mt()<< endl;
        csmp_error.notice( ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit",
                        "\nmt() below zero...!");
      }
  }


  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::CheckBoundaryFlags( Node<dim>* n )
  {
    if(t.Flag() == DIRICH && n->AtBoundary() != NOT && n->AtBoundary() != INTERNAL)
      {
        if (!open_boundaries)
          {
            fixed_temperature = true;
            t_fixed           = t();
          }
        else
          fixed_temperature = false;	    
      }
    else
      fixed_temperature = false;
  }


  // JPW Nov 2010
  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::TemperatureDependentHeatCapacityRock( double64 cpr_min_ext, double64 t_min_ext,
                                                                                      double64 cpr_max_ext, double64 t_max_ext )
  {
    T_dependent_cpr = true;
    //      equilibrator.TemperatureDependentHeatCapacityRock(cpr_min_ext,t_min_ext,cpr_max_ext,t_max_ext);
  }

  // JPW  Nov 2010
  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::WithOpenBoundaries( double64 reference_specific_enthalpy,
                                                                    double64 reference_salinity )
  {
    open_boundaries = true;
    const_top_pressure = true;
    ref_spec_h = reference_specific_enthalpy;
    ref_sal    = reference_salinity*0.01;
  }

  // JPW  Nov 2010
  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::WithOpenBoundaries( double64 reference_salinity )
  {
    open_boundaries = true;
    const_top_pressure = false;
    ref_sal    = reference_salinity*0.01;
  }

  // JPW  Nov 2010
  /*  template< size_t dim>
      double64 NaClH2OPropertiesVisitorPHX<dim>::BoundaryFlow()
      {

      double64 bfm_l(0.), bfm_v(0.), bfm_h(0.);

      bfe() = 0.;
      bfs() = 0.;

      if (definitelyGreaterThan(bfm(),0.,numeric_limits<double64>::epsilon()) &&
      definitelyLessThan(sh(),1.,numeric_limits<double64>::epsilon()))
      {
      krl = EffectiveLiquidSaturationHalitePresent( slp(), svp() );
      krv = 1.0 - krl;
      double64 mob_l(0.), mob_v(0.), mob(0.);
      if (slp()>0.) mob_l = krl*rlp()/mulp();
      if (svp()>0.) mob_v = krv*rvp()/muvp();
      mob = mob_l+mob_v;
      if (mob>0.) mob_l /= mob;
      if (mob>0.) mob_v /= mob;

      if (!definitelyGreaterThan(bfm(),(mlp()+mvp()),numeric_limits<double64>::epsilon()))
      {
      bfm_l = bfm()*mob_l;
      bfm_v = bfm()*mob_v;
      if (definitelyGreaterThan(bfm_v,mvp(),numeric_limits<double64>::epsilon()))
      {
      bfm_l += bfm_v-mvp();
      bfm_v = mvp();
      }
      else if (definitelyGreaterThan(bfm_l,mlp(),numeric_limits<double64>::epsilon()))
      {
      bfm_v += bfm_l-mlp();
      bfm_l = mlp();
      }
      }
      else
      {
      bfm_l = mlp();
      bfm_v = mvp();
      //            if (bfm()<mtp())
      bfm_h = mtp()-bfm_l-bfm_v;
      }

      bfe() = phi()*(bfm_l*hlp()+bfm_v*hvp()+bfm_h*hhp());
      bfs() = bfm_l*xlp()+bfm_v*xvp()+bfm_h;

      if (bfm()>mtp())
      {
      bfm_v = (bfm()-mtp())*mob_v;
      bfm_l = (bfm()-mtp())*mob_l;
      if (definitelyGreaterThan(bfm_v,dmv_,numeric_limits<double64>::epsilon()))
      {
      bfm_l += bfm_v-dmv_;
      bfm_v = dmv_;
      }
      else if (definitelyGreaterThan(bfm_l,dml_,numeric_limits<double64>::epsilon()))
      {
      bfm_v += bfm_l-dml_;
      bfm_l = dml_;
      }
      if (dml_>0.0)
      {
      bfe() += phi()*bfm_l*dhCl_/dml_;
      bfs() += bfm_l*dxCl_/dml_;
      }
      if (dmv_>0.0)
      {
      bfe() += phi()*bfm_v*dhCv_/dmv_;
      bfs() += bfm_v*dxCv_/dmv_;
      }
      }

      }
      if (definitelyLessThan(bfm(),0.,numeric_limits<double64>::epsilon()))
      {
      bfe() = bfm()*ref_spec_h*phi();
      bfs() = bfm()*ref_sal;
      }
      mt()  -= bfm();
      m_fluid_ = mt() * pore_volume;
      H_current_ -= bfe();
      if (bfs() < ms())
      ms() -= bfs();
      else
      {
      bfs() = ms();
      ms() = 0.;
      }
      wt() = ms()/mt()*100.;
      if (definitelyGreaterThan(wt(),100.,numeric_limits<double64>::epsilon()))
      wt() = 100.;
      wt_ = wt();

      if (mt() < 0. || H_current_ < 0. || ms() < 0.)
      {
      cout << "subtrating too much" << endl;
      cout << "mt(): " << mt() << endl;
      cout << "H_current_: " << H_current_ << endl;
      cout << "ms(): " << ms() << endl;
      cout << "bfm(): " << bfm() << endl;
      cout << "bfe(): " << bfe() << endl;
      cout << "bfs(): " << bfs() << endl;
      }

      if (time_tracking) tracker.Start("equilibration");
      equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
      if (time_tracking) tracker.Stop("equilibration");

      //    cout << "Sl: " << Liquid.s << endl;

      return mt()/Bulk.rho;

      }
  */
  // JPW  Oct 2011 correctipon
  template<size_t dim>
  double64 NaClH2OPropertiesVisitorPHX<dim>::BoundaryFlow()
  {

    double64 bfm_l(0.), bfm_v(0.), bfm_h(0.);
    double64 old_liquid(0.), old_vapor(0.), old_salt(0.);
    double64 new_liquid(0.), new_vapor(0.);

    if (definitelyGreaterThan(mtp(),mlp()+mvp(),numeric_limits<double64>::epsilon()))
      old_salt   = mtp() - mlp() - mvp();
    old_liquid = mlp() + std::min(0.,dml_);
    old_vapor  = mvp() + std::min(0.,dmv_);
    new_liquid = std::max(0.,dml_);
    new_vapor  = std::max(0.,dmv_);
    double64 total(old_salt+old_liquid+old_vapor+new_liquid+new_vapor);


    if (!essentiallyEqual(mt(),total,numeric_limits<double64>::epsilon()*total))
      {
        cout << "old_salt: " << old_salt << endl;
        cout << "old_liquid: " << old_liquid << endl;
        cout << "old_vapor: " << old_vapor << endl;
        cout << "new_liquid: " << new_liquid << endl;
        cout << "new_vapor: " << new_vapor << endl;
        cout << "total: " << total << endl;
        cout << "mt(): " << mt() << endl;
        cout << "diff: " << mt()-total << endl;
        cout << "." << endl;
      }

    bfe() = 0.;
    bfs() = 0.;

    if (definitelyGreaterThan(bfm(),0.,numeric_limits<double64>::epsilon()) &&
        definitelyLessThan(sh(),1.,numeric_limits<double64>::epsilon()))
      {
        krl = EffectiveLiquidSaturationHalitePresent( slp(), svp() );
        krv = 1.0 - krl;
        double64 mob_l(0.), mob_v(0.), mob(0.);
        if (slp()>0.) mob_l = krl*rlp()/mulp();
        if (svp()>0.) mob_v = krv*rvp()/muvp();
        mob = mob_l+mob_v;
        if (mob>0.) mob_l /= mob;
        if (mob>0.) mob_v /= mob;

        if (!definitelyGreaterThan(bfm(),(old_liquid+old_vapor),numeric_limits<double64>::epsilon()))
          {
            bfm_l = bfm()*mob_l;
            bfm_v = bfm()*mob_v;
            if (definitelyGreaterThan(bfm_v,old_vapor,numeric_limits<double64>::epsilon()))
              {
                bfm_l += bfm_v-old_vapor;
                bfm_v = old_vapor;
              }
            else if (definitelyGreaterThan(bfm_l,old_liquid,numeric_limits<double64>::epsilon()))
              {
                bfm_v += bfm_l-old_liquid;
                bfm_l = old_liquid;
              }
          }
        else
          {
            bfm_l = old_liquid;
            bfm_v = old_vapor;
            bfm_h = std::min(old_salt,bfm()-old_liquid-old_vapor);
          }

        bfe() = phi()*(bfm_l*hlp()+bfm_v*hvp()+bfm_h*hhp());
        bfs() = bfm_l*xlp()+bfm_v*xvp()+bfm_h;

        double64 old_mass(bfm_l+bfm_v+bfm_h);

        if (definitelyGreaterThan(bfm(),old_mass,numeric_limits<double64>::epsilon()))
          {
            bfm_v = (bfm()-old_mass)*mob_v;
            bfm_l = (bfm()-old_mass)*mob_l;
            if (definitelyGreaterThan(bfm_v,new_vapor,numeric_limits<double64>::epsilon()))
              {
                bfm_l += bfm_v-new_vapor;
                bfm_v = new_vapor;
              }
            else if (definitelyGreaterThan(bfm_l,new_liquid,numeric_limits<double64>::epsilon()))
              {
                bfm_v += bfm_l-new_liquid;
                bfm_l = new_liquid;
              }
            if (new_liquid>0.0)
              {
                //                bfe() += phi()*bfm_l*dhCl_/new_liquid;
                bfe() += bfm_l*dhCl_/new_liquid;
                bfs() += bfm_l*dxCl_/new_liquid;
              }
            if (dmv_>0.0)
              {
                //                bfe() += phi()*bfm_v*dhCv_/new_vapor;
                bfe() += bfm_v*dhCv_/new_vapor;
                bfs() += bfm_v*dxCv_/new_vapor;
              }
          }

      }
    if (definitelyLessThan(bfm(),0.,numeric_limits<double64>::epsilon()))
      {
        bfe() = bfm()*ref_spec_h*phi();
        bfs() = bfm()*ref_sal;
      }
    mt()  -= bfm();
    m_fluid_ = mt() * pore_volume;
    H_current_ -= bfe();
    if (bfs() < ms())
      ms() -= bfs();
    else
      {
        bfs() = ms();
        ms() = 0.;
      }
    wt() = ms()/mt()*100.;
    if (definitelyGreaterThan(wt(),100.,numeric_limits<double64>::epsilon()))
      wt() = 100.;
    wt_ = wt();

    if (mt() < 0. || H_current_ < 0. || ms() < 0.)
      {
        cout << "subtrating too much" << endl;
        cout << "mt(): " << mt() << endl;
        cout << "H_current_: " << H_current_ << endl;
        cout << "ms(): " << ms() << endl;
        cout << "bfm(): " << bfm() << endl;
        cout << "bfe(): " << bfe() << endl;
        cout << "bfs(): " << bfs() << endl;
      }

    equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);

    return mt()/Bulk.rho;

  }

  // JPW  Nov 2010
  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::BoundaryIteration( )
  {
    //    cout << "enter. " << endl;

    double64 mt_backup, ms_backup, H_backup;
    double64 bfm_min, bfm_max;
    double64 vol, vol_min, vol_max;
    double64 crit;
    bool outflow(false);

    mt_backup = mt();
    ms_backup = ms();
    H_backup = H_current_;

    if (!const_top_pressure) ref_spec_h = reference_enthalpy_top();

    if ( essentiallyEqual( mt(), rho_bulk(), numeric_limits<double64>::epsilon()) )
      {
        bfm() = 0;
        bfe() = 0;
        bfs() = 0.;
        vol = 1.;
        equilibrator.ThreePhaseProperties(Bulk,Liquid,Vapor,Salt);
      }
    else
      {

        // determine outflow or inflow
        // get first guess based on old desity

        //***********************
        // old code:
        /*    if (definitelyGreaterThan(mt()-rho_bulk(), 0.0, numeric_limits<double64>::epsilon()))
              outflow = true;

              bfm() = 0;

              vol = BoundaryFlow();
              if (outflow && definitelyLessThan( vol, 1., numeric_limits<double64>::epsilon() ))
              {
              //      int temp;
              cout << "Problem at initial open top" << endl;
              cout << "Outflow but initial vol is smaller than one" << endl;
              cout << "mt: " << mt() <<endl;
              cout << "rho_bulk: " << rho_bulk <<endl;
              cout << "Bulk.rho: " << Bulk.rho <<endl;
              cout << "t: " << t() <<endl;
              cout << "p: " << p() <<endl;
              cout << "wt: " << wt() <<endl;
              //      cin >> temp;
              }
              else if (!outflow && definitelyGreaterThan( vol, 1., numeric_limits<double64>::epsilon() ))
              {
              //      int temp;
              cout << "Problem at initial open top" << endl;
              cout << "Inflow but initial vol is greater than one" << endl;
              cout << "mt: " << mt() <<endl;
              cout << "rho_bulk: " << rho_bulk <<endl;
              cout << "Bulk.rho: " << Bulk.rho <<endl;
              cout << "t: " << t() <<endl;
              cout << "p: " << p() <<endl;
              cout << "wt: " << wt() <<endl;
              //      cin >> temp;
              }*/
        // end old code
        //***********************

        bfm() = 0;
        vol = BoundaryFlow();

        if (definitelyGreaterThan( vol, 1., numeric_limits<double64>::epsilon()))
          outflow = true;

        if (outflow)
          {
            vol_min = vol;
            bfm_min = bfm();
          }
        else
          {
            vol_max = vol;
            bfm_max = bfm();
          }

        bfm() = (mt()-rho_bulk())/2.;
        //    if (bfm() > rho_bulk()*(slp()+svp()))
        //        bfm() = rho_bulk()*(slp()+svp());
        //    if (bfm() > (rlp()*slp()+rvp()*svp()))
        //        bfm() = (rlp()*slp()+rvp()*svp());
        //    if (bfm() > (ml()+mv()))
        //        bfm() = (ml()+mv());
        //    bfm() = mt()/100.;
        if (outflow && definitelyLessThan(bfm(),0.,numeric_limits<double64>::epsilon()))
          bfm() *= -1.;
        if (!outflow && definitelyGreaterThan(bfm(),0.,numeric_limits<double64>::epsilon()))
          bfm() *= -1.;
        vol = BoundaryFlow();

        // determine bfm_min and bfm_max
        if (outflow)
          {
            if (definitelyLessThan(vol,1.,numeric_limits<double64>::epsilon()))
              // maximum has been found
              {
                vol_max = vol;
                bfm_max = bfm();
              }
            else
              // more mass has to be taken out
              {
                bool find_vol_max(true);
                int count_vol_max(0);
                while (find_vol_max)
                  {
                    count_vol_max++;
                    bfm() *= 2.;
                    //          if (bfm() > (ml()+mv()))
                    //              bfm() = (ml()+mv());
                    mt() = mt_backup;
                    ms() = ms_backup;
                    H_current_ = H_backup;

                    if (bfm()>mt()) bfm() = mt() - 1.e-4;
                    //          cout << "outflow bfm(). " << endl;
                    vol = BoundaryFlow();

                    if (definitelyLessThan(vol,1.,numeric_limits<double64>::epsilon()) || count_vol_max > 100)
                      // found the maximum
                      {
                        vol_max = vol;
                        bfm_max = bfm();
                        find_vol_max = false;
                      }
                  } // end find_vol_max
              } // end vol < 1

            // minimum should be 0 if nothing is removed
            /*     bfm() = 0.;
                   mt() = mt_backup;
                   ms() = ms_backup;
                   H_current_ = H_backup;

                   //     cout << "bfm() = 0. " << endl;
                   vol = BoundaryFlow();

                   vol_min = vol;
                   bfm_min = bfm();*/
          }
        else // inflow
          {
            if (definitelyGreaterThan(vol,1.,numeric_limits<double64>::epsilon()))
              {
                vol_min = vol;
                bfm_min = bfm();
              }
            else
              {
                bool find_vol_min(true);
                while (find_vol_min)
                  {
                    bfm() *= 2.;
                    mt() = mt_backup;
                    ms() = ms_backup;
                    H_current_ = H_backup;

                    //      cout << "inflow bfm(). " << endl;
                    vol = BoundaryFlow();

                    if (definitelyGreaterThan(vol,1.,numeric_limits<double64>::epsilon()))
                      {
                        vol_min = vol;
                        bfm_min = bfm();
                        find_vol_min = false;
                      }
                  }
              }

            /*        bfm() = 0.;
                      mt() = mt_backup;
                      ms() = ms_backup;
                      H_current_ = H_backup;

                      //        cout << "bfm() = 0. " << endl;
                      vol = BoundaryFlow();

                      vol_max = vol;
                      bfm_max = bfm();*/
          }


        if (definitelyGreaterThan( bfm_min, bfm_max, numeric_limits<double64>::epsilon() ) ||
            definitelyLessThan( vol_min, 1., numeric_limits<double64>::epsilon() ) ||
            definitelyGreaterThan( vol_max, 1., numeric_limits<double64>::epsilon() ))
          {
            cout << "Problem in boundary flow calculations." << endl;
            cout << "bfm_min: " << bfm_min << endl;
            cout << "bfm_max: " << bfm_max << endl;
            cout << "bfm_max-bfm_min: " << bfm_max-bfm_min << endl;
            cout << "vol_min: " << vol_min << endl;
            cout << "vol_max: " << vol_max << endl;
            cout << "t: " << t() <<endl;
            cout << "p: " << p() <<endl;
            cout << "wt: " << wt() <<endl;
            //        int temp;
            //        cin >> temp;
          }

        //    cout << "iteration starts" << endl;

        crit = abs(1.-vol);

        int count(0);

        // Bisection iteration
        //   while (crit > 1.e-4 &&
        //          essentiallyEqual( bfm_min, bfm_max, numeric_limits<double64>::epsilon()))
        while (crit > 1.e-4)
          {
            count++;

            bfm()  = bfm_min*(1.-vol_max);
            bfm() += bfm_max*(vol_min-1.);
            bfm() /= (vol_min - vol_max);

            //   cout << "bfm (guess): " << bfm() << endl;

            mt() = mt_backup;
            ms() = ms_backup;
            H_current_ = H_backup;

            vol = BoundaryFlow();

            crit = abs(1.-vol);

            if (vol>1.)
              {
                bfm_min = bfm();
                vol_min = vol;
              }
            if (vol<1.)
              {
                bfm_max= bfm();
                vol_max = vol;
              }
            if (count > 100)
              {
                cout << "Couldn't find good value for boundary flow after 100 iterations" << endl;
                cout << "bfm_min: " << bfm_min << endl;
                cout << "bfm_max: " << bfm_max << endl;
                cout << "bfm_max-bfm_min: " << bfm_max-bfm_min << endl;
                cout << "vol_min: " << vol_min << endl;
                cout << "vol_max: " << vol_max << endl;
                cout << "t: " << t() <<endl;
                cout << "p: " << p() <<endl;
                cout << "wt: " << wt() <<endl;
                //       int temp;
                //	    cin >> temp;
                crit = 0.;
              }
          }

        /*    if (vol > 1.01 || vol < 0.99)
              {
              cout << "Problem" << endl;
              cout << "t: " << t() <<endl;
              cout << "p: " << p() <<endl;
              cout << "wt: " << wt() <<endl;
              int temp;
              cin >> temp;
              }*/
      }
    //    cout << "exit. " << endl;

  }

  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::TimeStepAdjustment()
  {
    time_factor_p = 1.0;
    // next two unused for the moment
    time_factor_h = 1.0;
    time_factor_s = 1.0;      

    const double64 threshold_p(0.1);

    expected_dp = 1.0/mt()/Bulk.beta * ( mt()-Bulk.rho ); // if bulk.rho < mt -> pressure increase (squeeze back); if bulk.rho > mt() -> pressure decrease (expand)
    time_factor() =  1.0;

    if( fabs( expected_dp/p() ) > threshold_p )
      {
        time_factor() = threshold_p / fabs(expected_dp/p());

        /*    	if( time_factor > 1.0 )
                {
                cout << "time_factor() > 0 !!!\nCan't be!!!\n";
                cout << "time_factor() = " << time_factor() << endl;
                cout << "threshold_p   = " << threshold_p   << endl;
                cout << "expected_dp   = " << expected_dp   << endl;
                cout << "p()           = " << p   << endl;
                //    	    char mychar;
                //    	    cin >> mychar;
                }
                if( time_factor < 1e-6 )
                {
                cout << "WARNING: time_factor() < 0.01 !!!\n";
                cout << "time_factor() = " << time_factor() << endl;
                cout << "threshold_p   = " << threshold_p   << endl;
                cout << "expected_dp   = " << expected_dp   << endl;
                cout << "p()           = " << p   << endl;
                char mychar;
                //  	    cin >> mychar;
                }*/
      }
    else
      time_factor() = 1.0;
  }



  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::VolumeFactorComputations( Node<dim>* n )
  {
    rl_transport() = rl();
    rv_transport() = rv();

    // calculate mf
    mh()  = Salt.mf * mt();
    mf()  = mt() - mh();

    if(mf() < 1.0e-10)
      {
        cout <<"\nNode: ("<<n->x()<<", "<<n->y()<<", mf goes to zero: mf = "<<mf()
             << ", rhol = " << rl() << ", rhov = " << rv() << ", phi = " << phi() << endl;
        Output_PTXState();
      }

    // liquid + vapor mass
    mlv_               = ml() + mv();
    if (mlv_ > 0.) volume_factor_RHS  = mf()/mlv_;
    else volume_factor_RHS  = 0.;

    //    if(pure_halite)
    //      {
    //	cout<<"\nPure halite case: reseting volume_factor_RHS to 1.0"<< endl;
    //	volume_factor_RHS = 1.0;
    //	Output_PTXState();
    //      }

    volume_factor_LHS  = mt()/Bulk.rho;
    // in case of pure water, volume_factor_RHS = volume_factor_LHS?

    if(volume_factor_RHS < 0.0)
      {
        cerr<<"\nNode: "<< n->Idx()<<", volume_Factor_RHS less than zero: " << volume_factor_RHS << " for SoWat conditions:\n";
        ScreenOutputSowatVariables();
        csmp_error.notice( FATAL_ERROR, "NaClH2OPropertiesVisitorPHX<dim>::Visit()",
                        "\nNegative Volume Factor, terminating...!");
      }

    //    volume_factor = volume_factor_RHS;
    volume_factor = volume_factor_LHS;

    // new source term
    // Philipp Weis, 4 Feb 2010

    // version without dt -> has to be treated accordingly with the PointSource-operator
    //	nQ   = (mt()-Bulk.rho)/dt_;
    nQ()   = (mt()-Bulk.rho);

    //+++++++++++++++++++
    // Big Baustelle JPW

    /*    if (p.Flag() == PLAIN)
          {

          if(p() > 38.0e6 && p() < 40.e6 &&
          t() > 590. && t() < 605.)
          {
          double64 factor(fabs(p()-39.e6)/1.0e6);
          if (nQ() > 0. && nQ() > factor*factor*mt() )//&&
          {
          nQ() = factor*factor*mt();
          }
          }
          if (p()+expected_dp-1.0e5 < 0.)//&&
          {
          nQ() = -(p()-1.0e5)*mt()*Bulk.beta;
          }
          }*/

    // new quick fix to avoid extreme pressure changes
    /*    double64 old_nQ = nQ();
    //    if (p() > 25.e6 && p.Flag() == PLAIN && sl() < 1. )
    if (p() > 25.e6 && p.Flag() == PLAIN)
    {
    if (nQ() > threshold_p*p()*mt()*Bulk.beta)
    {
    nQ() = threshold_p*p()*mt()*Bulk.beta;
    }
    if (nQ() < -threshold_p*p()*mt()*Bulk.beta && sl() < 1.)
    {
    nQ() = -threshold_p*p()*mt()*Bulk.beta;
    }
    }*/

    //+++++++++++++++++++

    // The PointSource-operator doesn't multiply the variable with the Volume
    // so it has to be done here.
    // The old version multiplied with phi() which gives phi*Volume = PoreVolume
    // NOTE: if a different PDE-operator is used the old verion might be needed
    // (that's the case for NumIntegral_NT_op_N_dV for example) 

    //	nQ() *=  phi(); //old version
    //	nQ() *=  fv_man.ControlVolume(n->Idx()-1U).CVPoreVolume(); //new version
    //    nQ() *=  cpv(); //new version without FV Manager
    nQ() *=  cpv()*(1.-sh()); //new version without FV Manager
    //    if (pure_halite) nQ = 0.0;
	
    if(Bulk.rho>0. && mt()>0.)
      {
        // Compress LHS transport variables
        ml() *= volume_factor_LHS;
        mv() *= volume_factor_LHS;

        hCl() *= volume_factor_LHS;
        hCv() *= volume_factor_LHS;

        xCl() *= volume_factor_LHS;
        xCv() *= volume_factor_LHS;

        // Compress RHS transport variables
        rl_transport() *= volume_factor_RHS;
        rv_transport() *= volume_factor_RHS;

        hVl() *= volume_factor_RHS;
        hVv() *= volume_factor_RHS;

        xVl() *= volume_factor_RHS;
        xVv() *= volume_factor_RHS;

      }

    vol_fac() = volume_factor;
  }




  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::PrepareVariablesForStorage()
  {
    state()   = static_cast<double64>(Bulk.state);

    Htp()     = H_current_;
    msp()     = ms();
    mtp()     = mt();
    rho_bulk()= Bulk.rho;

    // nodal total compressibility
    CT()  = beta()*pore_volume*(1.-sh());
    CT() += beta_rock()*(rock_volume+pore_volume*sh());
    //    CT() *= mt();
    CT() *= rho_bulk();

    // new JPW
    // calculating nodal variables
    // liquid/vapor mass/enthalpy mobilities
    
    mml() = mmv() = mmld() = mmvd() = 0.0;
    eml() = emv() = rvl() = rvv() = 0.0;
    emld() = emvd() = xml() = xmv() = 0.0;
    //    krl = EffectiveLiquidSaturationHalitePresent( sl(), sv() );
    //    krv = 1.0 - krl;
    krl = RelativePermeabilityLiquid( sl(), sv() );
    krv = 1.0 - sh() - krl;

    if (sl() > 0.0) mml() = krl * rl_transport() / mul();
    if (sv() > 0.0) mmv() = krv * rv_transport() / muv();
    //    if (sl() > 0.0) mmld() = mml() * rl();
    //    if (sv() > 0.0) mmvd() = mmv() * rv();
    if (sl() > 0.0) mmld() = mml() * rl_transport();
    if (sv() > 0.0) mmvd() = mmv() * rv_transport();
    if (sl() > 0.0) eml() = krl * hVl() / mul();
    if (sv() > 0.0) emv() = krv * hVv() / muv();
    //    if (sl() > 0.0) emld() = eml() * rl();
    //    if (sv() > 0.0) emvd() = emv() * rv();
    if (sl() > 0.0) emld() = eml() * rl_transport();
    if (sv() > 0.0) emvd() = emv() * rv_transport();
    if (sl() > 0.0) xml() = krl * xVl() / mul();
    if (sv() > 0.0) xmv() = krv * xVv() / muv();
    if (sl() > 0.0) rvl() = krl / mul();
    if (sv() > 0.0) rvv() = krv / muv();
  }






  //*** new TD May 2011

  template<size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::CheckPhaseChange()
  {
    if( (old_state == VL  && (current_state == L || current_state == F || current_state == LH))
        ||
        (old_state == VH  && (current_state == L || current_state == F || current_state == LH || current_state == VL || current_state == H))
        ||
        (old_state == VLH && (current_state == L || current_state == F || current_state == LH || current_state == VL || current_state == H))
        )
      {
        dangerous_phase_change()    = 1.0; // = true
        after_phasechange_counter() = 1.0;
        beta_ref()                  = beta_p();
        //   cout << "dangerous phase change!" << endl;
        return;
      }
    else if( old_state == current_state && after_phasechange_counter() > 0.0 )
      {
        dangerous_phase_change()     = 1.0; // = true
        after_phasechange_counter() += 1.0;
        if(after_phasechange_counter() > max_after_phasechange_counter)
          {
            dangerous_phase_change()    = 0.0; // false
            after_phasechange_counter() = 0.0;
          }
        return;
      }
    else return;
  }

  template< size_t dim>
  void NaClH2OPropertiesVisitorPHX<dim>::CheckVolumeMismatchCompensation()
  {
    if( ( (mt()-Bulk.rho)/Bulk.beta/mt() / p() < -0.05 || p()+(mt()-Bulk.rho)/Bulk.beta/mt() < 1.0e5 )
        && (dangerous_phase_change() > 0.5) )
      {
        //      cout << "beta(): " << beta();
        beta() += (max_after_phasechange_counter - after_phasechange_counter())/max_after_phasechange_counter * beta_ref();
        beta() *= 0.5;
        //   cout << " changed to: " << beta() << endl;

      }
    return;
  }
  //*** end new

  template class NaClH2OPropertiesVisitorPHX<1>;
  template class NaClH2OPropertiesVisitorPHX<2>;
  template class NaClH2OPropertiesVisitorPHX<3>;

}

