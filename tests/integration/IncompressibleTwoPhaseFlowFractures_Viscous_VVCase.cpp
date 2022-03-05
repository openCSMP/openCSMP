#include "IncompressibleTwoPhaseFlowFractures_Viscous_VVCase.h"
#include "ModelComparator.h"
#include "VSet.h"


using namespace std;

namespace csmp {
template<uint32_t dim>
IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::IncompressibleTwoPhaseFlowFractures_Viscous_VVCase(const char* prefix, const char* twophase_model, const char* explicit_implicit, const char* first_second_order):
tpncfvt_(NULL),
relperm_model_(NULL),
model_(NULL),
steady_state_pressure_solver_(NULL)
{
    ///////////////////////////////////////////////////////////////////
    /// Set the name of the test and the name of variables file

    string test_name("IncompressibleTwoPhaseFlowFractures_Viscous_VVCase");
    variables_file_prefix_ = test_name;
    prefix_ = prefix;
    input_file_name_ = prefix;

    switch(dim){
    case 1:
        test_name += "_1D";
        break;
    case 2:
        test_name += "_2D";
        break;
    case 3:
        test_name += "_3D";
    }

    test_name = test_name +"_"+explicit_implicit+"_"+first_second_order;

    this->setName(test_name);

    ////////////////////////////////////////////////////////////////////

	if(string("implicit") == explicit_implicit){

        explicit_scheme_ = false;
        nonlinear_scheme_ = true;
        //nonlinear_scheme_ = false; // old scheme

	}else{

        explicit_scheme_ = true;
        nonlinear_scheme_ = true;    // doesn't matter false or true for explicit scheme

	}

	if(string("2nd order") == first_second_order){

        second_order_in_space_ = true;
        second_order_in_time_  = false;
        with_lsm_grad_limiter_ = false;

	}else{

        second_order_in_space_ = false;
        second_order_in_time_  = false;
        with_lsm_grad_limiter_ = false;

	}

    //====================================================================

    relperm_model_name_ = twophase_model;

    if(string("Corey Model") == twophase_model)
        relperm_model_type_ = COREY_MODEL;
    else if(string("Linear Two Phase Model") == twophase_model)
        relperm_model_type_ = LINEAR_TWOPHASE_MODEL;
    else if  (string("Brooks Corey Model") == twophase_model)
        relperm_model_type_ = BROOKS_COREY_MODEL;
    else if (string("Van Genuchten Model") == twophase_model )
        relperm_model_type_ = VAN_GENUCHTEN_MODEL;
    else{

        cout<< "\n Your Two Phase Model :"<< twophase_model<<" is not in the list!"<<endl;
        cout<< "\n The list of available Two Phase Models:"
            << "\n 01) Linear Two Phase Model"
            << "\n 02) Brooks Corey Model"
            << "\n 03) Corey Model"
            << "\n 04) Van Genuchten Model"
            << "\n";

        throw csmp::Exception( ERROR, "IncompressibleTwoPhaseFlowFractures_Viscous_VVCase",
                               "Two phase model not recognized" );
    }

    //====================================================================


    with_capillary_forces_ = false;
    with_gravity_forces_ = false;


    //Debug options
    debug_ = false;   
    assign_dirichlet_saturation_on_the_right_boundaries_ = false;
    lu_solver_ = false;
}

template<uint32_t dim>
IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::~IncompressibleTwoPhaseFlowFractures_Viscous_VVCase()
{
    delete model_;
    delete relperm_model_;
    delete tpncfvt_;
    delete steady_state_pressure_solver_;
}


//******************************************************************************************************************************
//
// Load Model Setup

template<uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::LoadModel()
{

    variables_file_cutted_prefix_ = variables_file_prefix_.substr(0, variables_file_prefix_.find("VVCase"))+"VVCase";

    switch(dim){
        case 2:
            model_ = dynamic_cast<Model<dim>*> ( new ANSYS_Model2D(input_file_name_.data(),(variables_file_cutted_prefix_+".txt").c_str(),false,true) );
            break;
        case 3:
            model_ = dynamic_cast<Model<dim>*> ( new ANSYS_Model3D(input_file_name_.data(),variables_file_cutted_prefix_.c_str(),(variables_file_cutted_prefix_+".txt").c_str()) );
    }

}


//******************************************************************************************************************************
//
// Variables Setup

template<uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::VariablesSetup(){

    // fluid and rock properties
    porosity_ ="porosity";
    visc_n_ = "viscosity oil";
    visc_w_ = "viscosity water";
    rho_n_ = "density oil";
    rho_w_ = "density water";
    total_mobility_ = "total mobility";
    permeability_  = "permeability";
    brooks_corey_parameter_ = "brooks corey parameter";
    entry_pressure_ = "entry pressure";
    frac_apperture_ = "fracture aperture";

    // flow properties
    fluid_pressure_ = "fluid pressure";
    cap_diffusivity_ = "capillary diffusivity";
    sw_ = "saturation water";
    swr_ = "residual saturation water";
    snr_ = "residual saturation oil";
    sn_ = "saturation oil";
    nodal_fluid_volume_source_ = "nodal fluid volume source";
    fluid_volume_source_  = "fluid volume source";

    total_velocity_ = "total velocity";
    velocity_n_ = "velocity oil";
    velocity_w_ = "velocity water";
    total_pore_velocity_ = "total pore velocity";
    total_volume_flux_ = "total volume flux";
    thickness_ = "thickness";

    no_flow_bc_ = fluid_pressure_.c_str();

    porosity_idx_ = model_->Database().StorageKey( porosity_.c_str());
    rho_w_idx_ =  model_->Database().StorageKey(rho_w_.c_str());
    rho_n_idx_ = model_->Database().StorageKey(rho_n_.c_str());

    total_mobility_idx_ = model_->Database().StorageKey( total_mobility_.c_str() );
    permeability_idx_ =  model_->Database().StorageKey( permeability_.c_str() );

    fluid_pressure_idx_=  model_->Database().StorageKey(fluid_pressure_.c_str());
    cap_diffusivity_idx_ = model_->Database().StorageKey( cap_diffusivity_.c_str() );
    sw_idx_ = model_->Database().StorageKey(sw_.c_str());
    sn_idx_ = model_->Database().StorageKey(sn_.c_str());
    fluid_volume_source_idx_ = model_->Database().StorageKey( fluid_volume_source_.c_str() );
    nodal_fluid_volume_source_idx_ = model_->Database().StorageKey( nodal_fluid_volume_source_.c_str() );

    total_velocity_idx_ = model_->Database().StorageKey( total_velocity_.c_str() );
    velocity_n_idx_ = model_->Database().StorageKey( velocity_n_.c_str() );
    velocity_w_idx_ = model_->Database().StorageKey( velocity_w_.c_str() );
    thickness_idx_ = model_->Database().StorageKey(thickness_.c_str());

}


//******************************************************************************************************************************
//
// Model Setup

template<uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::ModelSetup(){


    inlet_boundary_name_ = "LEFT";
    outlet_boundary_name_ = "RIGHT";

    const double frac_aperture=0.001;
    const double matrix_porosity=0.15;
    const double fracture_porosity=1.0*frac_aperture;
    const double visc_nonwet=1.0e-3;
    const double visc_wet=1.0e-3;
    const double dens_wet(1000.);
    const double dens_nonwet(800.);
    const double Km(1.0e-15);
    const double Kf(8.33e-8);// acording to power law with Kf=a^2/12 with a=0.001
    const double swres=0.0;
    const double snres=0.0;

    const double entry_pres(0.0);

    double lambda=2.;     //Brooks Corey

    const double wet_sat(0.0);
    const double nonwet_sat(1.0);
    const double fluid_pres(0.0);
    const double fluid_pres_LB(1.0e+7);
    const double fluid_pres_RB(1.0e+5);

    csmp::Region<dim>& mat_ref( model_->Region( "MATRIX" ) );
    csmp::Region<dim>& frac_ref( model_->Region( "FRACS" ) );

	// Model Configuring
    mat_ref.InputPropertyValue( porosity_.c_str(), makeScalar( PLAIN, matrix_porosity ) );
    frac_ref.InputPropertyValue( porosity_.c_str(), makeScalar( PLAIN, fracture_porosity ) );
    model_->InputPropertyValue( visc_w_.c_str(), makeScalar( PLAIN, visc_wet ) );
    model_->InputPropertyValue( visc_n_.c_str(), makeScalar( PLAIN, visc_nonwet ) );
    model_->InputPropertyValue( rho_w_.c_str(), makeScalar( PLAIN, dens_wet ) );
    model_->InputPropertyValue( rho_n_.c_str(), makeScalar( PLAIN, dens_nonwet ) );

    mat_ref.InputPropertyValue( permeability_.c_str(), makeScalar( PLAIN, Km ) );
    frac_ref.InputPropertyValue( permeability_.c_str(), makeScalar( PLAIN, Kf ) );
    model_->InputPropertyValue( entry_pressure_.c_str(), makeScalar( PLAIN, entry_pres ) );
    model_->InputPropertyValue( brooks_corey_parameter_.c_str(), makeScalar( PLAIN, lambda ) );   //Brooks Corey
    model_->InputPropertyValue( thickness_.c_str(), makeScalar( PLAIN, 1.0 ) );          
    frac_ref.InputPropertyValue( thickness_.c_str(), makeScalar( PLAIN, frac_aperture ) );          

    model_->InputPropertyValue( swr_.c_str(), makeScalar( PLAIN, swres ) );
    model_->InputPropertyValue( snr_.c_str(), makeScalar( PLAIN, snres ) );
    model_->InputPropertyValue( fluid_pressure_.c_str(), makeScalar( PLAIN, fluid_pres ) );

    // -----------------------------
    // Initial conditions
    // -----------------------------

    model_->InputPropertyValue( sw_.c_str(), makeScalar( PLAIN, wet_sat ) );
    model_->InputPropertyValue( sn_.c_str(), makeScalar( PLAIN, nonwet_sat ) );
    model_->InputPropertyValue( nodal_fluid_volume_source_.c_str(), makeScalar( PLAIN,0.0) );     // no sources&sinks
    model_->InputPropertyValue( fluid_volume_source_.c_str(), makeScalar( PLAIN,0.0) );			// no sources&sinks

    // -----------------------------
    // Boundary conditions
    // -----------------------------

	model_->Boundary( inlet_boundary_name_.c_str() ).InputPropertyValue( fluid_pressure_.c_str(), makeScalar( DIRICH, fluid_pres_LB ) );
	model_->Boundary( outlet_boundary_name_.c_str() ).InputPropertyValue( fluid_pressure_.c_str(), makeScalar( DIRICH, fluid_pres_RB ) );

	model_->Boundary( inlet_boundary_name_.c_str() ).InputPropertyValue( sw_.c_str(), makeScalar( DIRICH, 1.0 ) );
	model_->Boundary( inlet_boundary_name_.c_str() ).InputPropertyValue( sn_.c_str(), makeScalar( DIRICH, 0.0 ) );

	if(assign_dirichlet_saturation_on_the_right_boundaries_){
		model_->Boundary( outlet_boundary_name_.c_str() ).InputPropertyValue( sw_.c_str(), makeScalar( DIRICH, 0.0 ) );
		model_->Boundary( outlet_boundary_name_.c_str() ).InputPropertyValue( sn_.c_str(), makeScalar( DIRICH, 1.0 ) );
	}
	
	
    switch( relperm_model_type_ )
	{

	case LINEAR_TWOPHASE_MODEL:
		
        relperm_model_ = new LinearTwoPhaseModel<dim> ( model_->Database(), permeability_.c_str(),
                                                                          visc_n_.c_str(),visc_w_.c_str(),
                                                                          rho_n_.c_str(), rho_w_.c_str(),
                                                                          entry_pressure_.c_str(),
                                                                          sw_.c_str(),snr_.c_str(), swr_.c_str());
        lambda =0;
		break;

	default:

        relperm_model_ = new BrooksCorey<dim>( model_->Database(), permeability_.c_str(),
                                                                 visc_n_.c_str(), visc_w_.c_str(),
                                                                 rho_n_.c_str(), rho_w_.c_str(),
                                                                 brooks_corey_parameter_.c_str(), entry_pressure_.c_str(),
                                                                 sw_.c_str(),snr_.c_str(), swr_.c_str());

	}

	//Update saturations and total mobility
    UpdateSaturations( *relperm_model_ );
    ComputeTotalMobility( *relperm_model_ );

}




//******************************************************************************************************************************
//
// Velocities calculations


template<uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::UpdateSaturations( TwoPhaseModel<dim>& saturationFunctions )
 {
    ScalarVariable  sc, mob_t, thickness;

    csmp::Region<dim>& mref( model_->Region( "Model" ) );

    // updating the saturation of water
    const auto nodesEnd( mref.NodesEnd() );
    for ( auto it = mref.NodesBegin(); it != nodesEnd; ++it )
      {
         sc() = 1. - (*it)->Read( sn_idx_ );
         (*it)->Store(  sw_idx_, sc );
      }

} // updateSaturationsAndComputeTotalMobility

template<uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::ComputeTotalMobility( TwoPhaseModel<dim>& saturationFunctions )
 {
    ScalarVariable  sc, mob_t, thickness;

    csmp::Region<dim>& mref( model_->Region( "Model" ) );

    // computing the total mobility
    const auto elementsEnd( mref.ElementsEnd() );
    for ( auto it = mref.ElementsBegin(); it != elementsEnd; ++it )
      {
         // setting up the relative permeability model
         saturationFunctions.Initialize( *(*it) );
         saturationFunctions.InitializeForBaryCenter( *(*it) );
         saturationFunctions.EffectiveSaturation();

         // total mobility
         mob_t() = saturationFunctions.TotalMobility();
         (*it)->Store(  total_mobility_idx_, mob_t );
      }

  printRangeOfVariable( *model_, model_->Database().Name( total_mobility_idx_ ) );

} // updateSaturationsAndComputeTotalMobility


/* Incompressible Two Phase Flow Verification Test

  =================================
  2D:
  Mesh:       UnitSquareFracs_xzplane || UnitSquareFracs_yzplane || UnitSquareFracs_orthogonal || UnitSquareFracs_irregular
  Model:      1x1 m^2 2D
  Test:       vset_compare Test
  BC:         Dirichlet ("fluid pressure") on the LEFT & RIGHT boundaries
  Criterion:  Comparison with vset obtained with previous csmp version: saturation water

  =================================

  =================================
  3D:
  Mesh:       UnitCubeFracs_xzplane || UnitCubeFracs_yzplane || UnitCubeFracs_orthogonal || UnitCubeFracs_irregular
  Model:      1x1x1 m^3 3D
  Test:       vset_compare Test
  BC:         Dirichlet ("fluid pressure") on the LEFT & RIGHT boundaries
  Criterion:  Comparison with vset obtained with previous csmp version: saturation water

  =================================
  */

template <uint32_t dim>
void IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<dim>::run()
{

    //##################################################

    // Phi dS/dt + div (u* fw(S))= qw

    // S is the water saturation
    // fw is the water fractional flow function
    // u is the total velocity
    // qw is water volumetric source = 0
    //##################################################


    //*/****************************************************************************************************
    // Simulation parameters

    double dt;
    int32_t save_frequency;

    // max time increment: 100 days
    const double MAX_TIME(86400. * 100.);

    switch(dim){

    case 1:
			
		cout<<"IncompressibleTwoPhaseFlowFractures_Viscous_VVCase for 2D and 3D only"<<endl;

    case 2:
			
			// Time Settings
			// CFL (30.-200. s)
            dt = 1.0e3;
            max_time_ = 14.*dt;

            if(explicit_scheme_){

                time_increment_ = 0.25*dt;
                cfl_multiplier_ = 0.4; //notice: for TwoPhaseExplicitNodeCenteredFiniteVolumetransport cfl_multiplier is set to 0.4 (if no capillary effect) and to 0.1 (if there is capillary effect)

                // Output Settings
                save_frequency= 4;

            }else{

                time_increment_ = 0.25*dt;
                cfl_multiplier_ = 0.1;

                // Output Settings
                save_frequency= 4;

            }
            
            break;
    case 3:

			// Time Settings
            // CFL (20.-50. s)
            dt = 1.0e3;
            max_time_ = 14.*dt;

            if(explicit_scheme_){

                time_increment_ = 0.1*dt;
                cfl_multiplier_ = 0.4; //notice: for TwoPhaseExplicitNodeCenteredFiniteVolumetransport cfl_multiplier is set to 0.4 (if no capillary effect) and to 0.1 (if there is capillary effect)

                // Output Settings
                save_frequency= 10;

            }else{

                time_increment_ = 0.25*dt;
                cfl_multiplier_ = 1.0;

                // Output Settings
                save_frequency= 4;

            }
    
    }



    //*/****************************************************************************************************
    // Geometry

    cout<<"\n"<<this->getName()<<" :"<<" Building model.."<<endl;

	LoadModel();

	VariablesSetup();

	ModelSetup();

    cout<<"\n"<<this->getName()<<" :"<<" Finished Building Model..."<<endl;

    //*/*****************************************************************************************************


	// Setting the propertires for output

	list<string> vtuInitialOutputProps;
    vtuInitialOutputProps.push_back(sw_.c_str());
    vtuInitialOutputProps.push_back(sn_.c_str());
    vtuInitialOutputProps.push_back(fluid_pressure_.c_str());
    vtuInitialOutputProps.push_back(total_velocity_.c_str());
    
	list<string> vtuOutputProps;
    vtuOutputProps.push_back(sw_.c_str());
    vtuOutputProps.push_back(sn_.c_str());
    vtuOutputProps.push_back(fluid_pressure_.c_str());
    vtuOutputProps.push_back(total_velocity_.c_str());

    VTU_Interface<dim> vtu( *model_ );
    //vtu.OmitZeroInFileName(true);

    // ------------------------------------------------------------------------------------------
    // 1.0 Setting up an FE algorithm to solve the diffusion equation 0 = div(lambda_t grad p)
    //
    //     lambda_t = total mobility
    //     p = fluid pressure
    // ------------------------------------------------------------------------------------------

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Steady state pressure solver

    if(lu_solver_){

         steady_state_pressure_solver_ = new PDE_Integrator<dim,Region>  ( );

    }else{

        #ifdef CSMP_WITH_SAMG_SOLVER

        steady_state_pressure_solver_ = new PDE_Integrator<dim,Region>  ( new SAMG_Solver(&steady_state_pressure_solver_settings_));

        // Solver Settings

        // Select SAMG instance
        steady_state_pressure_solver_settings_.SetSolverInstance(0);
        // Re-use solver setup from previous timestep, without internal checks forcing a new setup when required
        // iswit(4) is required for first call (without purging memory) to store previous SAMG setup information
        steady_state_pressure_solver_settings_.Set_iswit(4);

        steady_state_pressure_solver_settings_.Set_ncgtyp(5);

        // nxtyp: 1,2=ILU, 0,5 = Gauss Seidel
        if( dim!=1 )
            steady_state_pressure_solver_settings_.Set_nxtyp(0);
        else
            steady_state_pressure_solver_settings_.Set_nxtyp(1);

        steady_state_pressure_solver_settings_.Set_ndefault(40);

        steady_state_pressure_solver_settings_.Set_iout1(0);
        steady_state_pressure_solver_settings_.Set_iout2(0);

        // SAMG solution criterion
        steady_state_pressure_solver_settings_.Set_eps(0.);
        steady_state_pressure_solver_settings_.Set_rel_eps(1.E-10);

        // Agressive first level coarsening nredlev(1) for decreased setup time and reduced no. of cycles
        //samg_settings.Set_nred(1);

        // Pre-adjust SAMG coarse matrix size relative to original size, based on solver output
        steady_state_pressure_solver_settings_.Set_a_cmplx(2);

        // Pre-adjust SAMG mesh complexity, based on solver output
        steady_state_pressure_solver_settings_.Set_g_cmplx(1.5);
        steady_state_pressure_solver_settings_.Set_w_avrge(2);

        // First approximation u=0
        steady_state_pressure_solver_settings_.Set_itypu(1);

      #else

        steady_state_pressure_solver_ = new PDE_Integrator<dim,Region>  ( new CSMP_DEFAULT_LINEAR_SOLVER() );

      #endif
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //LHS:
    NumIntegral_dNT_op_dN_dV<dim,Element<dim> >  conductance( model_->Database(), total_mobility_.c_str(), fluid_pressure_.c_str(),fluid_pressure_.c_str() );
    //RHS:
    NumIntegral_NT_op_N_dV<dim,Element<dim> >    src( model_->Database(),  fluid_volume_source_.c_str(), fluid_pressure_.c_str() );

	// operation to compute velocity
    VelocityAndVolumeFlux<dim,Element<dim> >  velo( *model_,
                                      total_mobility_.c_str(),
                                      porosity_.c_str(),
                                      fluid_pressure_.c_str(), false,
                                      total_velocity_.c_str(),
                                      total_pore_velocity_.c_str(),
                                      total_volume_flux_.c_str());

    steady_state_pressure_solver_->Add( &conductance );
    steady_state_pressure_solver_->Add( &src );

	// add velocity calculation as post process
    steady_state_pressure_solver_->AddPostProcess( &velo );
    
    model_->Apply(*steady_state_pressure_solver_);

	// operation to compute velocity
    //ComputeTotalAndPhaseVelocities(&velo);

    // show results
    printRangeOfVariable( *model_, fluid_pressure_.c_str() );
    printRangeOfVariable( *model_, total_velocity_.c_str() );

    // average total velocity
    double vel(model_->Region("Model").Average(total_velocity_.c_str()));
    cout<<"Velocity="<<vel<<endl;

    // Output to VTU files
    vtu.OutputDataToVTU(  std::string(this->getName()+"_Initial").c_str(), vtuInitialOutputProps,"Model", 0.0 );

    // -------------------------------------------------------------
    // 2.0 Construct the finite volume grid and transport algorithms
    // -------------------------------------------------------------

    if(explicit_scheme_){

        //EXPLICIT SCHEME

        tpncfvt_=dynamic_cast<TwoPhaseExplicitNodeCenteredFVTransport<dim,ExplicitStencilProcessor>*>
                (new TwoPhaseExplicitNodeCenteredFVTransport<dim,ExplicitStencilProcessor>(
                                                                                             "Model",*model_,
                                                                                             porosity_.c_str(),
                                                                                             cap_diffusivity_.c_str(),
                                                                                             sw_.c_str(),
                                                                                             sn_.c_str(),
                                                                                             total_velocity_.c_str(),
                                                                                             nodal_fluid_volume_source_.c_str(),
                                                                                             second_order_in_space_,
                                                                                             with_capillary_forces_,
                                                                                             with_gravity_forces_,
                                                                                             no_flow_bc_.c_str()));

        dynamic_cast<TwoPhaseExplicitNodeCenteredFVTransport<dim,ExplicitStencilProcessor>*>(tpncfvt_)->DisableCapillarySpreading();

    }else{

        // IMPLICIT SCHEME

        tpncfvt_=dynamic_cast< TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>
                (new TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor> (         "Model", *model_,
                                                                                             second_order_in_space_,
                                                                                             second_order_in_time_,
                                                                                             porosity_.c_str(),
                                                                                             visc_n_.c_str(),
                                                                                             visc_w_.c_str(),
                                                                                             rho_n_.c_str(),
                                                                                             rho_w_.c_str(),
                                                                                             sw_.c_str(),
                                                                                             sn_.c_str(),
                                                                                             total_velocity_.c_str(),
                                                                                             nodal_fluid_volume_source_.c_str(),
                                                                                             with_capillary_forces_,
                                                                                             with_gravity_forces_,
                                                                                             no_flow_bc_.c_str(),
                                                                                             nonlinear_scheme_));

        #ifdef CSMP_WITH_SAMG_SOLVER
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().SetSolverInstance(1);
        // minimal screen output
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().Set_iout1( 0 );
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().Set_iout2( 0 );
        // solution accuracy
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().Set_eps(0.);
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().Set_rel_eps(1.E-10);
        // iswit(4) required to keep setup after first call (i.e., not to purge memory in DLL)
        dynamic_cast<TwoPhaseImplicitNodeCenteredFVTransport<dim,StencilProcessor>*>(tpncfvt_)->GetSolverSettings().Set_iswit(4);
        #endif
    }

    if(with_lsm_grad_limiter_)
        tpncfvt_->WithLsmGradientLimiter(*model_);

    // CFL increment:

    tpncfvt_->CFL_Multiplier(cfl_multiplier_);

    courant_increment_ = tpncfvt_->AnisotropicCourantIncrement( *relperm_model_, MAX_TIME );


    // -----------------------
    // 3.0 Transient loop
    // -----------------------


    //*/****************************************************************************************************

    double model_time;
    int32_t save_counter,time;
    model_time= 0.;
    time = 0;
    save_counter = 1;

    clock_t start = clock();

    cout<<"\nVelocity="<<vel<<endl;
    cout<<"\nMax Time="<<max_time_<<endl;
    cout<<"\nTime Increment="<<time_increment_<<endl;
    cout<<"\nCourant Increment="<<courant_increment_<<endl;
    cout<<"\nCFL Multiplier="<<cfl_multiplier_<<endl;

    while ( model_time < max_time_ )
    {

        // increment time
        model_time += time_increment_;

        // compute advection of phases
        tpncfvt_->TransportPhase(*relperm_model_, time_increment_ );

		//update pressure and velocity field
		UpdateSaturations( *relperm_model_ );
		ComputeTotalMobility( *relperm_model_ );
	    model_->Apply(*steady_state_pressure_solver_);

		// output variables
        if ( save_counter == save_frequency ) {

      time = model_time / dt;

            
			// output to VTU files
            vtu.OutputDataToVTU(  this->getName().c_str(), vtuOutputProps,"Model", time );
            save_counter = 0;

        }
        cout<<"\nModel Time: "<<model_time<<"\tTime Increment: "<<time_increment_<<endl;
        cout<<"\nMax Time: "<<max_time_<<endl;
        save_counter++;

        // runtime info
        cout <<"\nmain: RUNTIME : "<< model_time/dt << endl << endl;

    }

    if(dim==2){

        model_->OutputToBinaryFile(this->getName().c_str());//csmp binary results.


        ModelComparator<dim> comparitor;

        double shouldBeZero( comparitor.CompareVSetsRenumberedNodes((this->getName()+".vset").c_str(), (this->getName()+"_Comparison.vset").c_str(),
                                                        "saturation water", "saturation water", 
														"IncompressibleTwoPhaseFlowFractures_Viscous_VVCase.txt", "IncompressibleTwoPhaseFlowFractures_Viscous_VVCase.txt", true ) );

        VSet<dim> vset_comparison;
        double model_time_comparison(model_time);
        vset_comparison.InputFrom( (this->getName()+"_Comparison.vset").c_str(), model_time_comparison );
        const Model<dim> model_comparison( vset_comparison, "IncompressibleTwoPhaseFlowFractures_Viscous_VVCase.txt", true );
        VTU_Interface<dim> vtu_comparison( model_comparison);
        vtu_comparison.OutputDataToVTU( "IncompressibleTwoPhaseFlow_Viscous_Comparison", vtuOutputProps,"Model", time );


        cout<<"L2 norm of difference for comparison: "<<shouldBeZero<<endl;


        _equal( shouldBeZero, 0., 1.0E-2 );

	}
	
	clock_t end = clock();

    //*/****************************************************************************************************

    unsigned long millisec ((end - start) * 1000 / CLOCKS_PER_SEC);

    cout<<"\nIncompressibleTwoPhaseFlowFractures_Viscous_VVCase:"<<endl;
    cout<<"\nElapsed Time = "<<millisec<<" ms ("<<(double)(millisec)/1000.<<" sec; "<<(double)(millisec)/60000.<<" min; "<<(double)(millisec)/3600000.<<" hours)"<<endl;

    // terminate
    delete steady_state_pressure_solver_;
    cout << "\nThat's it..."<< endl;


} // end

template class IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<2U>;
template class IncompressibleTwoPhaseFlowFractures_Viscous_VVCase<3U>;

} // namespace csmp

