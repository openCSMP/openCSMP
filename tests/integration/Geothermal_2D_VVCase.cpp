#include "Geothermal_2D_VVCase.h"

using namespace std;


namespace csmp
{

  Geothermal_2D_VVCase::Geothermal_2D_VVCase(const char* prefix)
  {
    this->setName("Geothermal_2D_VVCase");
    prefix_=prefix;
  }

/** Geothermal 2D Test Case
=================================
Mesh:       Linear Triangles
Model:      2000x1000 m 2D
Test:       Temperature
BC:         constant pressure, constant temperature
Criterion:
=================================
*/

  void Geothermal_2D_VVCase::outputToVTU( Model<3U>& model,std::string model_name, const list<string>& props, size_t timestep )
  {
      static VTU_Interface<3U> vtu(model);
      vtu.OutputDataToVTU( ( string(model_name) + "_Properties" ).c_str(), props, model.Region("Model"), timestep);
  }

  void Geothermal_2D_VVCase::ComputeMassConductivity (Model<3U>& model)
  {
      PropertyHandle<3U> rho_ph ( model, "density liquid", SCALAR, NODE );
      PropertyHandle<3U> kappa_ph ( model, "conductivity", SCALAR, ELEMENT );
      PropertyHandle<3U> lambda_ph ( model, "mass conductivity", SCALAR, ELEMENT );
      
      ConductivityVisitor<3U> conductivity_visitor( model, "conductivity", "permeability", "fluid viscosity" );
      model.Accept(conductivity_visitor);

      lambda_ph = 0.;
      lambda_ph += rho_ph;
      lambda_ph *= kappa_ph;
  }

  void Geothermal_2D_VVCase::ComputeMassGravityTerm (Model<3U>& model)
  {
      ComputeGravityTermVisitor<3U> gravity_visitor ( model, "gravity term", "permeability", 
                                                       "fluid viscosity", "density liquid" );
      model.Accept (gravity_visitor);
      
      Index gravityVectorKey  = model.Database().StorageKey ("gravity term");
      Index massGravityVectorKey  = model.Database().StorageKey ("mass gravity term");
      Index fluidDensityKey  = model.Database().StorageKey ("density liquid");
      
      ScalarVariable rho;
      VectorVariable<3U> gravityVector;
      
      const vector<Element< 3U>*>::const_iterator modelElementsEnd( model.Region("Model").ElementsEnd() );
      
      for( vector<Element<3U>*>::const_iterator it( model.Region("Model").ElementsBegin() ); 
           it != modelElementsEnd; ++it )
      {
        for (size_t ip=0;ip<(*it)->IntegrationPoints (); ++ip)
        {
          (*it)->Read(ip, gravityVectorKey, gravityVector);
          (*it)->PropertyValueAtIntegrationPoint( fluidDensityKey, ip, rho );
          (*it)->Store( ip, massGravityVectorKey, gravityVector*rho() );  
        }
      }
  }
  
  void Geothermal_2D_VVCase::run()
  {
    bool& globalVerbose( GlobalVerbose::Instance().globalVerbose );
    const  size_t DIM(3U);

    globalVerbose=true;

    // ---------------------------------
    // Finite Element Mesh Construction
    std::string geometry_name ("2000x1000_mesh");
	std::string regions_name (this->getName()); 
	std::string config_name (this->getName()); 
	std::string vars_name (this->getName()+".txt");
    
    
    ANSYS_Model3D model(geometry_name.c_str(), regions_name.c_str(), vars_name.c_str() , true, true, true);
    const PropertyDatabase<DIM>&  pd_ref(model.Database()); //reference to the models property database.
    
    printModelDimensions<DIM>(model, true );

    list<string> output_props;
    list<string> region_names;

    output_props.push_back("fluid pressure");
    output_props.push_back("velocity");
    output_props.push_back("temperature");
    output_props.push_back("fluid density");
    output_props.push_back("density liquid");
    output_props.push_back("volumetric enthalpy liquid");
    output_props.push_back("enthalpy content liquid");
    output_props.push_back("nodal fluid volume source");


    region_names.push_back("Model");

    // -----------------------------------------------------------------------
    // 1. Assigning material properties, initial conditions, and boundary
    //    conditions from *-configuration.txt files.
    // -----------------------------------------------------------------------
    InputDataManager<DIM>  model_configuration;
    ComputationalSettings  run_settings;
    
    // boolean flags set reading to: 2) default prop.values, 3) group prop.values, 4) essential conditions for box-shaped model
    
    model_configuration.ConfigureFromFile( model, config_name.c_str(), false,    // groupname from parameter range
                                                               true,    // default property values
                                                               true,    // regional property values
                                                               false,    // boundary conditions for box-shaped model
                                                               true,    // essential conditions for groups
                                                               true,    // csmp::Boundary properties
                                                               run_settings );
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Visitors ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Thermal Visitor
    ThermalVisitor<DIM>   thermal_equilibrator( model);
	SourceVisitor<DIM>    source_calculator(model);

    //! get nodal rock properties
    model.ExtrapolateElementToNodeProperty("porosity", "nodal porosity");
    model.ExtrapolateElementToNodeProperty("density rock", "nodal density rock");
    model.ExtrapolateElementToNodeProperty("heat capacity rock", "nodal heat capacity rock");
    model.ExtrapolateElementToNodeProperty("compressibility rock", "nodal compressibility rock");


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! FE  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //! to quickly replace SAMG with LU solver when necessary, uncomment
    #ifdef CSMP_WITH_SAMG_SOLVER
    //LUdcmp_Solver solver;
    SAMG_Solver solver;
    #else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    #endif

    //! steady state pressure
    PDE_Integrator<3U, Region>  steady_state_pressure( &solver );

    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  p_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_dNT_op_dV<DIM,Element<DIM> >     gravity( pd_ref, "mass gravity term", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<DIM,Element<DIM> >    zero_fluid_src( pd_ref, "fluid pressure" );

    VelocityAndVolumeFlux<DIM,Element<DIM> >  velocity( model, "conductivity", "porosity", "fluid pressure", "density liquid", false );
    //VelocityAndVolumeFlux<DIM,Element<DIM> >  velocity( model, "conductivity", "porosity", "fluid pressure", false );

    steady_state_pressure.Add( &p_conductance );                           
    steady_state_pressure.Add( &gravity );    
    steady_state_pressure.Add(&zero_fluid_src);
    steady_state_pressure.AddPostProcess( &velocity );
                                                                           
    //! transient pressure
    PDE_Integrator<3U, Region>  transient_pressure( &solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  pt_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    pt_conductance.MultiplyWithTimeIncrement(true);
    
    NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> > pt_capacitance_lhs( pd_ref, "total compressibility", "fluid pressure", "fluid pressure" );
    pt_capacitance_lhs.LumpedFormulation(true);
    
    NumIntegral_NT_op_N_dV<DIM,Element<DIM> >    pt_capacitance_rhs( pd_ref, "total compressibility",  "fluid pressure" );
    pt_capacitance_rhs.LumpedFormulation(true);
  
    NumIntegral_dNT_op_dV<DIM,Element<DIM> >     pt_gravity( pd_ref, "mass gravity term", "fluid pressure" );
    pt_gravity.MultiplyWithTimeIncrement(true);
    pt_gravity.AddAccumulateLater();
    
    // nodal fluid volume source for mass correction term
    PointSource_rhsop <DIM,Element<DIM> > fluid_src( pd_ref, "nodal fluid volume source", "fluid pressure" );
    fluid_src.MultiplyWithTimeIncrement(false);
    fluid_src.LumpedFormulation(true);
    fluid_src.AddAccumulateLater();
    
    VelocityAndVolumeFlux<DIM,Element<DIM> >  t_velocity( model, "conductivity", "porosity", "fluid pressure", "density liquid", false );
    //VelocityAndVolumeFlux<DIM,Element<DIM> >  t_velocity( model, "conductivity", "porosity", "fluid pressure", false );
    
    transient_pressure.Add( &pt_conductance );
    transient_pressure.Add( &pt_capacitance_lhs );
    transient_pressure.Add( &pt_capacitance_rhs );
    transient_pressure.Add( &pt_gravity );
    transient_pressure.Add( &fluid_src );
    transient_pressure.AddPostProcess( &t_velocity );
    
    //! temperature diffusion
    PDE_Integrator<DIM, Region>  temperature_diffusion( &solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >   t_conductance( pd_ref, "thermal conductivity", "temperature", "temperature" );
    
    NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> >  t_capacitance_lhs( pd_ref, "total heat capacity", "temperature", "temperature" );
    t_capacitance_lhs.LumpedFormulation(true);
    t_capacitance_lhs.MultiplyWithTimeIncrement(true);
    
    NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  t_capacitance_rhs( pd_ref, "total heat capacity", "temperature" );
    t_capacitance_rhs.MultiplyWithTimeIncrement(true);
    
    temperature_diffusion.Add( &t_conductance );
    temperature_diffusion.Add( &t_capacitance_lhs );
    temperature_diffusion.Add( &t_capacitance_rhs );
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! FV //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ExplicitMassBasedTransport<DIM, MassBasedStencilProcessor> mass_advection ("Model", model,
                                            "porosity",          //porosity    
                                            "fluid density",     //advected property lhs 
                                            "density liquid",    //advected property rhs 
                                            "velocity",          //transport velocity
                                            "nodal heat source", //nodal source
                                            false);              //second order in space

    ExplicitMassBasedTransport<DIM, MassBasedStencilProcessor> enthalpy_advection ("Model", model,
                                                    "porosity",                    //porosity    
                                                    "enthalpy content liquid",     //advected property lhs  
                                                    "volumetric enthalpy liquid",  //advected property  rhs
                                                    "velocity",                    //transport velocity
                                                    "nodal heat source",           //nodal source
                                                    false);                        //second order in space

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //! time
    double global_time        = 0.;
    double max_time           = 60*3.15e7; // 60 years
    double time_increment     = 3.15e7; // 1 year
	double time_increment_advection, time_advection;
    size_t time_step            = 0;
      
    //! Initial equalibration
    thermal_equilibrator.SetInitialProperties(&model );
    model.InterpolateNodeToElementProperty("nodal total heat capacity", "total heat capacity");
    model.InterpolateNodeToElementProperty("nodal total compressibility", "total compressibility");
	model.InterpolateNodeToElementProperty("density liquid", "density liquid element");
    
    ComputeMassConductivity(model);
    ComputeMassGravityTerm (model);
    model.Apply (steady_state_pressure);
    
    thermal_equilibrator.SetInitialProperties(&model );
    model.InterpolateNodeToElementProperty("nodal total heat capacity", "total heat capacity");
    model.InterpolateNodeToElementProperty("nodal total compressibility", "total compressibility");
	model.InterpolateNodeToElementProperty("density liquid", "density liquid element");
  
    //!important
    //! set mt and hCl (advected properties) to Dirich at the boundaries where p, t are dirichlet
    model.Boundary("LEFT").ChangePropertyStatus("fluid density", DIRICH);
    model.Boundary("RIGHT").ChangePropertyStatus("fluid density", DIRICH);
    model.Boundary("LEFT").ChangePropertyStatus("enthalpy content liquid", DIRICH);
    model.Boundary("RIGHT").ChangePropertyStatus("enthalpy content liquid", DIRICH);
  
    //! output initial equilibrated state
    outputToVTU( model, config_name, output_props, time_step );
    time_step++;
	cin.get();
  
    while ( global_time <= max_time )
    {
      temperature_diffusion.TimeIncrement(1./time_increment);
      model.Apply (temperature_diffusion);

      time_advection = 0.;
	  while(time_advection < time_increment)
	  {
		time_increment_advection = mass_advection.AnisotropicCourantIncrement();
		if (time_advection + time_increment_advection > time_increment)
			time_increment_advection = time_increment - time_advection;

	    mass_advection.AdvectVariable(time_increment_advection);
		enthalpy_advection.AdvectVariable(time_increment_advection);

		time_advection += time_increment_advection;
	  }


      model.Accept(thermal_equilibrator);
	  model.Accept(source_calculator);

      model.InterpolateNodeToElementProperty("nodal total heat capacity", "total heat capacity");
      model.InterpolateNodeToElementProperty("nodal total compressibility", "total compressibility");
	  model.InterpolateNodeToElementProperty("density liquid", "density liquid element");
      
      //! conductivity and gravity term depend on "density liquid"
      ComputeMassConductivity(model);
      ComputeMassGravityTerm (model);

      transient_pressure.TimeIncrement( time_increment );
      model.Apply (transient_pressure);

	  outputToVTU( model, config_name, output_props, time_step );

      global_time += time_increment;
      time_step ++;
    }

    
    return;
  }

}//end csmp
