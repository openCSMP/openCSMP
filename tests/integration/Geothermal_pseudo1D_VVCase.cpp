#include "Geothermal_pseudo1D_VVCase.h"

using namespace std;


namespace csmp
{

  Geothermal_pseudo1D_VVCase::Geothermal_pseudo1D_VVCase(const char* prefix)
  {
    this->setName("Geothermal_pseudo1D_VVCase");
    prefix_=prefix;
  }

/** Geothermal pseudo 1D Test Case
=================================
Mesh:       Linear Triangles
Model:      2000x1000 m 2D
Test:       Temperature
BC:         constant pressure, constant temperature
Criterion:  comparison with TOUGH
=================================
*/

  void Geothermal_pseudo1D_VVCase::outputToVTU( Model<2U>& model,string model_name, const list<string>& props, size_t timestep )
  {
      static VTU_Interface<2U> vtu(model);
      vtu.OutputDataToVTU( ( model_name + "_Properties" ).c_str(), props, model.Region("Model"), timestep);
  }

  void Geothermal_pseudo1D_VVCase::ComputeMassConductivity (Model<2U>& model)
  {
      PropertyHandle<2U> rho_ph ( model, "density liquid", SCALAR, NODE );
      PropertyHandle<2U> kappa_ph ( model, "conductivity", SCALAR, ELEMENT );
      PropertyHandle<2U> lambda_ph ( model, "mass conductivity", SCALAR, ELEMENT );
      
      ConductivityVisitor<2U> conductivity_visitor( model, "conductivity", "permeability", "fluid viscosity" );
      model.Accept(conductivity_visitor);

      lambda_ph = 0.;
      lambda_ph += rho_ph;
      lambda_ph *= kappa_ph;
  }


  bool Geothermal_pseudo1D_VVCase::Compare (Model<2U>& model, string file)
  {

    FILE* fin;
    double x, y, val, res, tol;
    size_t ind(0);
    size_t flag(0);
    size_t nPoints;
    ScalarVariable result( PLAIN, 0.0 );
    map<size_t, std::vector<double> > points;
    vector<double>  values;
    
    fin = fopen(file.c_str(), "r");
    if (fin == NULL)
    {
      printf("\n File could not be opened");
      return false;
    }

    // TOUGH simulation results
    // reads (x_coord, y_coord, temperature) data from a text file
    while (flag != EOF)
    {
      flag = fscanf( fin, "%lf\t%lf\t%lf\n", &x, &y, &val);
      points[ind].push_back(x);
      points[ind].push_back(y);
      values.push_back(val);
      ind ++;
    }

    fclose(fin);

    nPoints = points.size();

    // finds CSMP computed temperature values at given points
    PropertyAtPointVisitor<2U> pAt(model, points, "temperature");
    model.Accept(pAt);

    res = 0.;
    tol = 1.;

    // computes C-norm of the difference (sum of the absolute values)
    for (size_t ind2=0; ind2<nPoints; ind2++ )
    {
      pAt.PropertyValueAt(ind2, result);
      res += fabs(result() - values[ind2]);
    }

    // divided by the number of points
    res /= nPoints;
    
    //printf("\n res = %f", res);
    
    // if residual is smaller than tolerance, test passed
    if (res < tol)
      return true;
    else
      return false;
  }
  
  void Geothermal_pseudo1D_VVCase::run()
  {
    bool& globalVerbose( GlobalVerbose::Instance().globalVerbose );
    const  size_t DIM(2U);

    globalVerbose=true;

    // ---------------------------------
    // Finite Element Mesh Construction
    std::string geometry_name ("2000x1000_mesh"); 
	std::string regions_name (this->getName()); 
	std::string config_name (this->getName()); 
	std::string vars_name (this->getName()+".txt");
    
    ANSYS_Model2D model(geometry_name.c_str(), regions_name.c_str(), vars_name.c_str() , false, true );
    model.InstantiateFiniteVolumes();
    const PropertyDatabase<DIM>&  pd_ref(model.Database()); //reference to the models property database.
    
    printModelDimensions<DIM>(model, true );

    // output properties if necessary
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
    //output_props.push_back("conductivity");
    //output_props.push_back("mass conductivity");
    //output_props.push_back("permeability");
    //output_props.push_back("total compressibility");
    //output_props.push_back("porosity");
    //output_props.push_back("total heat capacity");
    //output_props.push_back("thermal conductivity");

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
                                                               true,    // boundary conditions for box-shaped model
                                                               true,    // essential conditions for groups
                                                               false,    // csmp::Boundary properties
                                                               run_settings );
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! Visitors ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Thermal Visitor
    ThermalVisitor<DIM>   thermal_equilibrator( model);
	SourceVisitor<DIM>    source_calculator(model);

    //! get nodal rock properties
    model.ExtrapolateCellToNodeProperty("porosity", "nodal porosity");
    model.ExtrapolateCellToNodeProperty("density rock", "nodal density rock");
    model.ExtrapolateCellToNodeProperty("heat capacity rock", "nodal heat capacity rock");
    model.ExtrapolateCellToNodeProperty("compressibility rock", "nodal compressibility rock");


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //! FE  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //! to quickly replace SAMG with LU solver when necessary, uncomment
    #ifdef CSMP_WITH_SAMG_SOLVER
    //LUdcmp_Solver solver;
    //SAMG_Settings settings;
    //SAMG_Solver solver(&settings);
    SAMG_Solver solver;

    #else
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    #endif
    
    //! steady state pressure
    PDE_Integrator<DIM, Region>  steady_state_pressure( solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  p_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<DIM,Element<DIM> >    zero_fluid_src( pd_ref, "fluid pressure" );

    VelocityAndVolumeFlux<DIM,Element<DIM> >  velocity( model, "conductivity", "porosity", "fluid pressure", false );

    steady_state_pressure.Add( &p_conductance );                           
    steady_state_pressure.Add(&zero_fluid_src);
	steady_state_pressure.AddPostProcess( &velocity );
                                                                           
    //! transient pressure
    PDE_Integrator<DIM, Region>  transient_pressure( solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >  pt_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    pt_conductance.MultiplyWithTimeIncrement(true);
    
    NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> > pt_capacitance_lhs( pd_ref, "total compressibility", "fluid pressure", "fluid pressure" );
    pt_capacitance_lhs.LumpedFormulation(true);
    
    NumIntegral_NT_op_N_dV<DIM,Element<DIM> >    pt_capacitance_rhs( pd_ref, "total compressibility",  "fluid pressure" );
    pt_capacitance_rhs.LumpedFormulation(true);
  
    // nodal fluid volume source for mass correction term
    PointSource_rhsop <DIM,Element<DIM> > fluid_src( pd_ref, "nodal fluid volume source", "fluid pressure" );
    fluid_src.MultiplyWithTimeIncrement(false);
    fluid_src.LumpedFormulation(true);
    fluid_src.AddAccumulateLater();
    
    VelocityAndVolumeFlux<DIM,Element<DIM> >  t_velocity( model, "conductivity", "porosity", "fluid pressure", false );
    
    transient_pressure.Add( &pt_conductance );
    transient_pressure.Add( &pt_capacitance_lhs );
    transient_pressure.Add( &pt_capacitance_rhs );
    transient_pressure.Add( &fluid_src );
    transient_pressure.AddPostProcess( &t_velocity );
    
    //! temperature diffusion
    PDE_Integrator<DIM, Region>  temperature_diffusion( solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >   t_conductance( pd_ref, "thermal conductivity", "temperature", "temperature" );
    
    NumIntegral_NT_lhsop_N_dV<DIM,Element<DIM> >  t_capacitance_lhs( pd_ref, "total heat capacity", "temperature", "temperature" );
    t_capacitance_lhs.LumpedFormulation(true);
    t_capacitance_lhs.MultiplyWithTimeIncrement(true);
    
    NumIntegral_NT_op_N_dV<DIM,Element<DIM> >  t_capacitance_rhs( pd_ref, "total heat capacity", "temperature" ); 
    t_capacitance_rhs.MultiplyWithTimeIncrement(true);
    t_capacitance_rhs.LumpedFormulation(true);

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
    double max_time           = 3.15e9*30; //3000 years  
    double time_increment     = 3.15e7; //1 year
	double time_increment_advection, time_advection;
	size_t time_step            = 0;
      

    
    //! Initial equalibration
    thermal_equilibrator.SetInitialProperties(&model );
    model.InterpolateNodeToCellProperty("nodal total heat capacity", "total heat capacity");
    model.InterpolateNodeToCellProperty("nodal total compressibility", "total compressibility");
	model.InterpolateNodeToCellProperty("density liquid", "density liquid element");
	
	   
    ComputeMassConductivity(model);
    model.Apply (steady_state_pressure);


    thermal_equilibrator.SetInitialProperties(&model );
    model.InterpolateNodeToCellProperty("nodal total heat capacity", "total heat capacity");
    model.InterpolateNodeToCellProperty("nodal total compressibility", "total compressibility");
	model.InterpolateNodeToCellProperty("density liquid", "density liquid element");
  
    //!important
    //! set mt and hCl (advected properties) to Dirich at the boundaries where p, t are dirichlet
    model.Boundary("LEFT").ChangePropertyStatus("fluid density", DIRICH);
    model.Boundary("RIGHT").ChangePropertyStatus("fluid density", DIRICH);
    model.Boundary("LEFT").ChangePropertyStatus("enthalpy content liquid", DIRICH);
    model.Boundary("RIGHT").ChangePropertyStatus("enthalpy content liquid", DIRICH);
  
      
    //! output initial equilibrated state
    time_step++;
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

      model.InterpolateNodeToCellProperty("nodal total heat capacity", "total heat capacity");
      model.InterpolateNodeToCellProperty("nodal total compressibility", "total compressibility");
	  model.InterpolateNodeToCellProperty("density liquid", "density liquid element");
      
      //! conductivity depends on "density liquid"
      ComputeMassConductivity(model);

      transient_pressure.TimeIncrement( time_increment );
      model.Apply (transient_pressure);

      // uncomment for an output
	  //if (time_step % 100 == 0)

      if (time_step == (size_t)1600)
      {
        if (Compare(model, this->getName()+"_Comparison_TOUGH_data_1600yr.txt"))
            _succeed();
        outputToVTU( model, config_name, output_props, time_step );
      }
      
      if (time_step == (size_t)1900)
      {
        if (Compare(model, this->getName()+"_Comparison_TOUGH_data_1900yr.txt"))
          _succeed();
        outputToVTU( model, config_name, output_props, time_step );
      }
      
      if (time_step == (size_t)3000)
      {
        if (Compare(model, this->getName()+"_Comparison_TOUGH_data_3000yr.txt"))
            _succeed();
        outputToVTU( model, config_name, output_props, time_step );
      }
      
      global_time += time_increment;
      time_step ++;
	  	  
    }

    std::cout << "Geothermal_Pseudo1D_VVCase Completed" << std::endl;
    
    return;
  }

}//end csmp
