#include "Geothermal_1D_VVCase.h"
#include "LineElementMesher.h"
#include "ModelTopology.h"
#include "Model.h"
#include "TextInterface.h"
#include "VTU_Interface.h"
#include "CSMP_highLevelUtilities.h"
#include "LinearSolver.h"
#include "GlobalVerbose.h"

// finite volumes
#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"

// finite elements
#include "PDE_Integrator.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "NumIntegral_SetRHS_to_Zero.h"
#include "PointSource_rhsop.h"
#include "VelocityAndVolumeFlux.h"

// visitors
//#include "ComputeGravityTermVisitor.h"
#include "ThermalVisitor.h"
#include "SourceVisitor.h"
#include "ConductivityVisitor.h"
#include "PropertyAtPointVisitor.h"

// utilities and monitoring
#include "InputDataManager.h"
#include "ComputationalSettings.h"

using namespace std;

namespace csmp {

Geothermal_1D_VVCase::Geothermal_1D_VVCase()
  {
    this->setName("Geothermal_1D_VVCase");
  }



void Geothermal_1D_VVCase::run()
  {
      bool& globalVerbose( GlobalVerbose::Instance().globalVerbose );
      enum { DIM=1 };

      globalVerbose=true;

     // 1. Creating a 1D mesh with 100 elements
     // --------------------------------------------
     VSet<DIM>    mesh_container;
     const double model_length{ 2000. };
     const int    n_elmts{ 100 };
     LineElementMesher<DIM>().BuildUniformMesh( mesh_container, model_length, n_elmts );

     // creating a model topological region
     set<string>    fem_types; fem_types.insert("ISOPARAMETRIC_LINEAR_BAR");
     ModelTopology  mesh_topology(true);  // isoparametric

     // making a rock region of elements
     vector<size_t>  elms;
     for ( size_t i{0U}; i<mesh_container.Elements(); i++ ) elms.push_back(i);
     mesh_topology.AddDomain( "ROCK", fem_types, elms );
     elms.erase( elms.begin(), elms.end() );
     // mesh_container.Out();

     // 2. Building the 1D Region named 'model'
     // --------------------------------------------
     string geometry_name("2000_mesh");
     string regions_name( this->getName());
     string config_name( this->getName());
     string vars_name( this->getName()+"-variables.txt");

     Model<DIM>  model( mesh_topology, mesh_container, vars_name.c_str(), true );
     const PropertyDatabase<DIM>&  pd_ref(model.Database()); //reference to the model's property database.
     printModelDimensions<DIM>(model, true );
     //model.Mesh().Out();
    
     // output properties if necessary
     list<string> output_props, region_names;

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
    
    model_configuration.ConfigureFromFile( model, config_name.c_str(),
                                           false, // do not read region name from parameter range
                                           true,  // read default property values
                                           true,  // regional property values
                                           false, // boundary conditions for box-shaped model
                                           true,  // essential conditions for groups
                                           true,  // csmp::Boundary properties
                                           run_settings );
                                           
    // extra "box" boundary conditions for 1D model (CNR1 = LEFT, CNR2 = RIGHT )
    Node<DIM>* const cnr1 = cornerFlaggedNode( model, CNR1 ); // LEFT
    Node<DIM>* const cnr2 = cornerFlaggedNode( model, CNR2 ); // LEFT
    assert( cnr1 );
    assert( cnr2 );
    
    const csmp::Index pf_key = model.Database().StorageKey("fluid pressure");
    const csmp::Index T_key  = model.Database().StorageKey("temperature");
    
    cnr1->Store( pf_key, makeScalar(DIRICH,202650.) );
    cnr1->Store( T_key, makeScalar(DIRICH,80.) );

    cnr2->Store( pf_key, makeScalar(DIRICH,101325.) );
    cnr2->Store( T_key, makeScalar(DIRICH,50.) );
    
    //! Thermal Visitor
    ThermalVisitor<DIM>   thermal_equilibrator( model);
	  SourceVisitor<DIM>    source_calculator(model);

    //! get nodal rock properties
    model.ExtrapolateCellToNodeProperty("porosity", "nodal porosity");
    model.ExtrapolateCellToNodeProperty("density rock", "nodal density rock");
    model.ExtrapolateCellToNodeProperty("heat capacity rock", "nodal heat capacity rock");
    model.ExtrapolateCellToNodeProperty("compressibility rock", "nodal compressibility rock");

    // -----------------------------------------------------------------------
    // 2. Pressure solver
    // -----------------------------------------------------------------------
    //! steady state pressure
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    PDE_Integrator<DIM,Region> steady_state_pressure( solver );
    
    NumIntegral_dNT_op_dN_dV<DIM,Element<DIM> >   p_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    NumIntegral_SetRHS_to_Zero<DIM,Element<DIM> > zero_fluid_src( pd_ref, "fluid pressure" );
    VelocityAndVolumeFlux<DIM,Element<DIM> >      velocity( model, "conductivity", "porosity", "fluid pressure", false );

    steady_state_pressure.Add( &p_conductance );
    steady_state_pressure.Add(&zero_fluid_src);
	  steady_state_pressure.AddPostProcess( &velocity );
                                                                           
    //! transient pressure
    PDE_Integrator<DIM,Region>  transient_pressure( solver );
    
    NumIntegral_dNT_op_dN_dV<DIM>  pt_conductance( pd_ref, "mass conductivity", "fluid pressure", "fluid pressure" );
    pt_conductance.MultiplyWithTimeIncrement(true);
    NumIntegral_NT_lhsop_N_dV<DIM> pt_capacitance_lhs( pd_ref, "total compressibility", "fluid pressure", "fluid pressure" );
    pt_capacitance_lhs.LumpedFormulation(true);
    NumIntegral_NT_op_N_dV<DIM>    pt_capacitance_rhs( pd_ref, "total compressibility",  "fluid pressure" );
    pt_capacitance_rhs.LumpedFormulation(true);
  
    // nodal fluid volume source for mass correction term
    PointSource_rhsop<DIM> fluid_src( pd_ref, "nodal fluid volume source", "fluid pressure" );
    fluid_src.MultiplyWithTimeIncrement(false);
    fluid_src.LumpedFormulation(true);
    fluid_src.AddAccumulateLater();
    
    VelocityAndVolumeFlux<DIM>  t_velocity( model, "conductivity", "porosity", "fluid pressure", false );
    
    transient_pressure.Add( &pt_conductance );
    transient_pressure.Add( &pt_capacitance_lhs );
    transient_pressure.Add( &pt_capacitance_rhs );
    transient_pressure.Add( &fluid_src );
    transient_pressure.AddPostProcess( &t_velocity );
    
    //! temperature diffusion
    PDE_Integrator<DIM, Region>     temperature_diffusion( solver );
    NumIntegral_dNT_op_dN_dV<DIM>   t_conductance( pd_ref, "thermal conductivity", "temperature", "temperature" );
    NumIntegral_NT_lhsop_N_dV<DIM>  t_capacitance_lhs( pd_ref, "total heat capacity", "temperature", "temperature" );
    t_capacitance_lhs.LumpedFormulation(true);
    t_capacitance_lhs.MultiplyWithTimeIncrement(true);
    
    NumIntegral_NT_op_N_dV<DIM>  t_capacitance_rhs( pd_ref, "total heat capacity", "temperature" );
    t_capacitance_rhs.MultiplyWithTimeIncrement(true);
    
    temperature_diffusion.Add( &t_conductance );
    temperature_diffusion.Add( &t_capacitance_lhs );
    temperature_diffusion.Add( &t_capacitance_rhs );
    

    // -----------------------------------------------------------------------
    // 2. Thermal advection / diffusion
    // -----------------------------------------------------------------------
    model.InstantiateFiniteVolumes(); // since there are no variables here with a FV indicative placement
    
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
    //! time
    double global_time        = 0.;
    double max_time           = 3.1536e9*30; //3000 years
    double time_increment     = 3.1536e7; //1 year
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
    const csmp::Index rhof_key = model.Database().StorageKey("fluid density");
    const csmp::Index Hl_key   = model.Database().StorageKey("enthalpy content liquid");
    cnr1->Status( rhof_key, DIRICH );
    cnr2->Status( rhof_key, DIRICH );
    cnr1->Status( Hl_key, DIRICH );
    cnr2->Status( Hl_key, DIRICH );  
      
    //! output initial equilibrated state
    // uncomment for an initial output
    //OutputToVTU( model, config_name, output_props, time_step );
    TextInterface  text_output;
    time_step++;
  
    while ( global_time <= max_time )
      {
        temperature_diffusion.TimeIncrement(1./time_increment);
        model.Apply(temperature_diffusion);
        
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
        // if (time_step % 100 == 0)
        //	outputToVTU( model, config_name, output_props, time_step );

        if (time_step == 1600UL ) {
            text_output.OutputDataAsTextColumns( model, "pressure1600", "fluid pressure" );
            text_output.OutputDataAsTextColumns( model, "temperature1600", "temperature" );
            if (Compare(model, this->getName()+"_Comparison_TOUGH_data_1600yr.txt"))
              _succeed();
          }
        
        if (time_step == 1900UL ) {
            text_output.OutputDataAsTextColumns( model, "pressure1900", "fluid pressure" );
            text_output.OutputDataAsTextColumns( model, "temperature1900", "temperature" );
            if (Compare(model, this->getName()+"_Comparison_TOUGH_data_1900yr.txt"))
              _succeed();
          }
        
        if (time_step == 3000UL ) {
            text_output.OutputDataAsTextColumns( model, "pressure3000", "fluid pressure" );
            text_output.OutputDataAsTextColumns( model, "temperature3000", "temperature" );
            if (Compare(model, this->getName()+"_Comparison_TOUGH_data_3000yr.txt"))
              _succeed();
          }
            
        global_time += time_increment;
        time_step ++;
      }
    
  } // end run



void Geothermal_1D_VVCase::OutputToVTU( Model<1U>& model,string model_name, const list<string>& props, size_t timestep )
  {
     VTU_Interface<1U> vtu(model);
     vtu.OutputDataToVTU( (model_name + "_Properties" ).c_str(), props, model.Region("Model"), timestep);
  }


void Geothermal_1D_VVCase::ComputeMassConductivity (Model<1U>& model)
  {
      PropertyHandle<1U> rho_ph ( model, "density liquid", SCALAR, NODE );
      PropertyHandle<1U> kappa_ph ( model, "conductivity", SCALAR, ELEMENT );
      PropertyHandle<1U> lambda_ph ( model, "mass conductivity", SCALAR, ELEMENT );
      
      ConductivityVisitor<1U> conductivity_visitor( model, "conductivity", "permeability", "fluid viscosity" );
      model.Accept(conductivity_visitor);

      lambda_ph = 0.;
      lambda_ph += rho_ph;
      lambda_ph *= kappa_ph;
  }



bool Geothermal_1D_VVCase::Compare( Model<1U>& model, string file )
  {
    FILE* fin;
    double x, y, val, res, tol;
    size_t ind(0);
    size_t flag(0);
    size_t nPoints;
    ScalarVariable result( PLAIN, 0.0 );
    map<size_t, vector<double> > points;
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
    PropertyAtPointVisitor<1U>  pAt(model, points, "temperature");
    model.Accept(pAt);

    res = 0.;
    tol = 1.;

    // computes C-norm of the difference (sum of the absolute values)
    for (size_t ind2=0; ind2<nPoints; ind2++ )
    {
      pAt.PropertyValueAt(ind2, result);
      res += fabs(result() - values[ind]);
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
  
} //end csmp
