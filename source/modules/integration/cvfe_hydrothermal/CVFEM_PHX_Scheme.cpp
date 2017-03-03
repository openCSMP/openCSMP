#include "CVFEM_PHX_Scheme.h"

using namespace std;

namespace csmp {


/** custom constructor */
template<size_t dim>
CVFEM_PHX_Scheme<dim>::CVFEM_PHX_Scheme( Model<dim>& model_ref, bool with_gravity )
   : model(model_ref),
     p_ref(model.Database()),
     verbose(false),
     move_on(false),
     check_consistency(false),
     open_top(false),
     time_tracking(true),
     brick_wall_limiter(true),
     diff_mass(model,names.diff_mass_variables[0].c_str(), SCALAR, NODE ),
     mass(model,names.diff_mass_variables[1].c_str(), SCALAR, NODE ),
     fluid_density(model,names.diff_mass_variables[2].c_str(), SCALAR, NODE ),
     diff_enthalpy(model,names.diff_enthalpy_variables[0].c_str(), SCALAR, NODE ),
     enthalpy(model,names.diff_enthalpy_variables[1].c_str(), SCALAR, NODE ),
     total_enthalpy(model,names.diff_enthalpy_variables[2].c_str(), SCALAR, NODE ),
     bfm(model,names.boundary_corrections[0].c_str(), SCALAR, NODE ),
     bfe(model,names.boundary_corrections[1].c_str(), SCALAR, NODE ),
     bfs(model,names.boundary_corrections[2].c_str(), SCALAR, NODE ),
     dhc(model,names.conduction_visitor_variables[1].c_str(), SCALAR, NODE ),
     pres_grad(model,names.pres_grad_variables[0].c_str(),
                  names.pres_grad_variables[1].c_str(),
                  names.pres_grad_variables[2].c_str(),
                  names.pres_grad_variables[3].c_str(),
                  names.pres_grad_variables[4].c_str()),
     pore_visitor(model,names.pore_visitor_variables[0].c_str(),
                     names.pore_visitor_variables[1].c_str(),
                     names.pore_visitor_variables[2].c_str()),
     fv_transport_vapor(model,
                        names.lhs_vapor_vector,
                        names.rhs_vapor_vector,
                        names.fv_transport_vapor_variables[0].c_str(),
                        names.fv_transport_vapor_variables[1].c_str(),
                        names.fv_transport_vapor_variables[2].c_str()),
     fv_transport_liquid(model,
                        names.lhs_liquid_vector,
                        names.rhs_liquid_vector,
                        names.fv_transport_liquid_variables[0].c_str(),
                        names.fv_transport_liquid_variables[1].c_str(),
                        names.fv_transport_liquid_variables[2].c_str()),
     upwind_control( model,fv_transport_liquid,fv_transport_vapor,
                        names.upwind_control_variables[0].c_str(),
                        names.upwind_control_variables[1].c_str(),
                        names.densities,
                        names.relperm_visc,
                       names.saturations,
                       names.time_control,
                        names.velocities,
                        names.pore_velocities),
     transport( model,upwind_control,
                fv_transport_vapor,fv_transport_liquid),
    #ifdef CSMP_WITH_SAMG_SOLVER
    T_SAMG_solver(&T_settings),
    T_FE_SAMG(&T_SAMG_solver),
    p_SAMG_solver(&p_settings),
    p_FE_SAMG(&p_SAMG_solver),
    #else
     /// add extra functionality for alternative solver if needed
     T_FE_SAMG(&T_LINEAR_solver),
     p_FE_SAMG(&p_LINEAR_solver),
    #endif
     capacitance_lhs(p_ref,names.capacitance_lhs_variables[0].c_str(),
                           names.capacitance_lhs_variables[1].c_str(),
                           names.capacitance_lhs_variables[2].c_str()),
     conductance(p_ref,names.conductance_variables[0].c_str(),
                       names.conductance_variables[1].c_str(),
                       names.conductance_variables[2].c_str()),
     capacitance_rhs(p_ref,names.capacitance_rhs_variables[0].c_str(),
                           names.capacitance_rhs_variables[1].c_str()),
     heat_bottom(p_ref,names.heat_bottom_variables[0].c_str(),
                       names.heat_bottom_variables[1].c_str()),
     capacitance_lhs_p(p_ref,names.capacitance_lhs_p_variables[0].c_str(),
                             names.capacitance_lhs_p_variables[1].c_str(),
                             names.capacitance_lhs_p_variables[2].c_str()),
     conductance_p_upwind_liquid(p_ref,upwind_control,fv_transport_liquid,
                                       names.conductance_p_upwind_liquid_variables[0].c_str(),
                                       names.conductance_p_upwind_liquid_variables[1].c_str(),
                                       names.conductance_p_upwind_liquid_variables[2].c_str(),
                                       names.conductance_p_upwind_liquid_variables[3].c_str(),
                                       names.conductance_p_upwind_liquid_variables[4].c_str()),
     conductance_p_upwind_vapor(p_ref,upwind_control,fv_transport_vapor,
                                       names.conductance_p_upwind_vapor_variables[0].c_str(),
                                       names.conductance_p_upwind_vapor_variables[1].c_str(),
                                       names.conductance_p_upwind_vapor_variables[2].c_str(),
                                       names.conductance_p_upwind_vapor_variables[3].c_str(),
                                       names.conductance_p_upwind_vapor_variables[4].c_str()),
     capacitance_rhs_p(p_ref,names.capacitance_rhs_p_variables[0].c_str(),
                             names.capacitance_rhs_p_variables[1].c_str()),
     grav_liq(p_ref,upwind_control,fv_transport_liquid,
                    names.grav_liq_variables[0].c_str(),
                    names.grav_liq_variables[1].c_str(),
                    names.grav_liq_variables[2].c_str(),
                    names.grav_liq_variables[3].c_str()),
     grav_vap(p_ref,upwind_control,fv_transport_vapor,
                    names.grav_vap_variables[0].c_str(),
                    names.grav_vap_variables[1].c_str(),
                    names.grav_vap_variables[2].c_str(),
                    names.grav_vap_variables[3].c_str()),
     source_p(p_ref,names.source_p_variables[0].c_str(),
                    names.source_p_variables[1].c_str()),
     source_p2(p_ref,names.source_p2_variables[0].c_str(),
                     names.source_p2_variables[1].c_str()),
     equilibrator_properties(model),
//     mass_visitor(model,names.mass_visitor_variables[0].c_str(),
//                        names.mass_visitor_variables[1].c_str()),
//     enthalpy_visitor(model,names.enthalpy_visitor_variables[0].c_str(),
//                            names.enthalpy_visitor_variables[1].c_str()),
//     conduction_visitor(model,names.conduction_visitor_variables[0].c_str(),
//                            names.conduction_visitor_variables[1].c_str()),
     pressure_limiter_fluid(p_ref,"fluid pressure",1.01325e5,490.e6),
     pressure_limiter_transport(p_ref,"fluid pressure",0.,1000.e6),
     largest_timestep(60.*60.*24.*365.),
     timestep(0)
 { 

  equilibrator_properties.SetAdjustCompressibilityAfterPhaseChangeBoolTo(true);

#ifdef CSMP_WITH_SAMG_SOLVER
  T_settings.Set_iout1( -1 );
//  T_settings.Set_iout1( 2 );
  T_settings.Set_iout2( -1 );
//  T_settings.Set_iout2( 0 );
  T_settings.Set_idmp( -1 );
//  T_settings.Set_idmp( 8 );
  T_settings.Set_ncgrad(0);

  p_settings.Set_iout1( -1 );
  p_settings.Set_iout2( -1 );
  p_settings.Set_idmp( -1 );
//  p_settings.Set_iout1( 2 );
//  p_settings.Set_iout2( 0 );
  p_settings.Set_eps(-1.e-5);
  p_settings.Set_nred(1);
  p_settings.Set_a_cmplx(3.);
  p_settings.Set_g_cmplx(2.);
  p_settings.Set_w_avrge(2.);
#else
    /// add extra functionality for alternative solver if needed
#endif

  // create PropertyHandles for full reset variables
  
  for (size_t i = 0U; i < names.full_reset_variables.size(); i++)
     reset_properties.push_back(
        new PropertyHandle<dim>(model,("reset "+names.full_reset_variables[i]).c_str(),SCALAR,NODE));

  dt = cfl_dt = current_dt = control_dt = old_dt = largest_timestep;

  upwind_control.Gravity( with_gravity );
  if (with_gravity) transport.WithGravityComponentLiquidAndVapor();

  capacitance_lhs.LumpedFormulation(true);
  capacitance_rhs.LumpedFormulation(true);

  heat_bottom.AddAccumulateLater();

  conductance.MultiplyWithTimeIncrement(true);
  heat_bottom.MultiplyWithTimeIncrement(true);

  capacitance_lhs_p.LumpedFormulation(true);
  capacitance_rhs_p.LumpedFormulation(true);

  source_p.AddAccumulateLater();
  source_p2.AddAccumulateLater();

  conductance_p_upwind_liquid.MultiplyWithTimeIncrement(true);
  conductance_p_upwind_vapor.MultiplyWithTimeIncrement(true);
  grav_liq.MultiplyWithTimeIncrement(true);
  grav_liq.AddAccumulateLater();
  grav_vap.MultiplyWithTimeIncrement(true);
  grav_vap.AddAccumulateLater();
  source_p2.MultiplyWithTimeIncrement(true);

  // pressure

  if (with_gravity) p_FE_SAMG.Add( &grav_liq );
  if (with_gravity) p_FE_SAMG.Add( &grav_vap );
  p_FE_SAMG.Add( &conductance_p_upwind_vapor );
  p_FE_SAMG.Add( &conductance_p_upwind_liquid );
  p_FE_SAMG.Add( &source_p );

  p_FE_SAMG.Add( &capacitance_lhs_p );

  p_FE_SAMG.Add( &capacitance_rhs_p );
  p_FE_SAMG.Add( &source_p2 );

  //temperature

  T_FE_SAMG.Add( &capacitance_lhs );
  T_FE_SAMG.Add( &conductance );

  T_FE_SAMG.Add( &capacitance_rhs );
  T_FE_SAMG.Add( &heat_bottom );

// CVFEM_Visitors are not ported yet!!!
// Do we want to keep functionality?

  // mass visitor
//  mass_visitor.Add( p_ref, &conductance_p_upwind_vapor );
//  mass_visitor.Add( p_ref, &conductance_p_upwind_liquid );
//  if (with_gravity) mass_visitor.Add( &grav_liq );
//  if (with_gravity) mass_visitor.Add( &grav_vap );
//  mass_visitor.Add( &source_p2 );

  // enthalpy visitor
//  enthalpy_visitor.Add( p_ref, &conductance_p_upwind_liquid, names.enthalpy_visitor_liquid_variables[0].c_str());
//  enthalpy_visitor.Add( p_ref, &conductance_p_upwind_vapor, names.enthalpy_visitor_vapor_variables[0].c_str());
// //  enthalpy_visitor.Add( p_ref, &conductance );
//  if (with_gravity) enthalpy_visitor.Add( p_ref,&grav_liq, names.enthalpy_visitor_liquid_variables[0].c_str());
//  if (with_gravity) enthalpy_visitor.Add( p_ref,&grav_vap, names.enthalpy_visitor_vapor_variables[0].c_str());
// //  enthalpy_visitor.Add( &heat_bottom );
//  enthalpy_visitor.Add( p_ref, &source_p2, names.enthalpy_visitor_source_variables[0].c_str() );

  // conduction visitor
//  conduction_visitor.Add( p_ref, &conductance );
//  conduction_visitor.Add( &heat_bottom );

  // calculate and store volume and pore volume
  model.Accept( pore_visitor );

    cout <<"\n\nCVFEM_PHX_Scheme<"<< typeid(double64).name() <<","<< dim;
    cout <<">: Constructed successfully."<< endl;

  } // end constructor 

/** custom destructor */
template<size_t dim>
CVFEM_PHX_Scheme<dim>::~CVFEM_PHX_Scheme()
 {
  for (size_t i = 0U; i < names.full_reset_variables.size(); i++)
     delete reset_properties[i];
 }

/** modifying maximum size of time step */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::SetLargestTimeStep(double64 timestep)
    {

      largest_timestep = dt = cfl_dt = current_dt = control_dt = old_dt = timestep;
      upwind_control.SetLargestTimeStep( largest_timestep );
      transport.SetLargestTimeStep( largest_timestep );

    } // end SetLargestTimeStep

/** adjusting timestep */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::ChangeTimeStepTo(double64 timestep)
    {

      dt = cfl_dt = current_dt = control_dt = old_dt = timestep;

    } // end SetLargestTimeStep

/** initialize fluid properties from current PTX conditions */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::InitialFluidPropertiesFromPTX()
    {

       model.Accept( pore_visitor );
       equilibrator_properties.InitialPropertiesFromPTX();
       model.CopyReplace(names.diff_mass_variables[2].c_str(),names.diff_mass_variables[1].c_str());
       model.CopyReplace(names.diff_enthalpy_variables[2].c_str(),names.diff_enthalpy_variables[1].c_str());

    } // end InitialFluidPropertiesFromPT

/** preparation before transient calculations, calculating pressure gradient and updwind nodes from current status */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::PrepareTransientCalculations()
    {

      model.Accept( pres_grad );
      transport.UpdateProjection();
      model.Accept( upwind_control );
      upwind_control.Reset();

    } // end PrepareTransientCalculations


/** main function to apply CVFEM scheme in transeint calculations, returns time step used for calculations */
template<size_t dim>
double64 CVFEM_PHX_Scheme<dim>::Apply()
    {

      timestep++;
      double64 reset_dt;

      AdvanceTransientVariables();

      current_dt  = std::min( old_dt*1.01,largest_timestep );

      reset_dt = current_dt;

//      p_FE_SAMG.TimeIncrement( 0.0 );
//      model.Pass(p_FE_SAMG);
//      model.Pass(pressure_limiter_fluid);
//      equilibrator_properties.SetTimeIncrement(0.0);
//      model.Accept(equilibrator_properties);
      
      bool vol_loop(true);     
      while (vol_loop)
      {

      // pressure equation and finite volume calculations
      AdvectionDiffusionLoops();

      // heat conduction
      T_FE_SAMG.TimeIncrement( current_dt );
      model.Apply( T_FE_SAMG );

// CVFEM Visitors are not ported yet!!!!
//      if (check_consistency) ApplyCVFEM_Visitors();

      // check for pressure below atmospheric
      model.MinMaxOf( "fluid pressure", min_value, max_value );

      if ( (min_value < 1.02e5 || max_value > 4.9e8) && brick_wall_limiter)
        {
//          cout << "\nApplying pressure limiter (fluid)" << endl;
          model.Apply(pressure_limiter_fluid);
        }

      FluidRockEquilibration();

// CVFEM Visitors not ported yet!!!!
//      if (check_consistency) CheckForConsistency();

       model.MinMaxOf( "volume factor", min_value, max_value );
       cout << "\nvolume factor range: " << min_value << " - " << max_value << endl;

       // optional slot for full reset
/*      if (min_value < 0.1 || max_value > 10.0)
        {
         current_dt *= 0.5;
         cout << "\ncurrent_dt: " << current_dt << endl;

         if (current_dt < 1.)
         {
           cout << "\ntimestep < 1. in Vol-Loop (current_dt: " << current_dt << ")." << endl;
         }
         cin >> temp;

           FullReset();
         }
        else*/
       vol_loop = false;

      }

     old_dt = current_dt;

     return current_dt;

    } // end Apply

/** book keeping of variables for transient calculations */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::AdvanceTransientVariables( )
  {
    for (size_t i = 0; i < names.transient_variables.size(); i++)
       model.CopyReplace( names.transient_variables[i].c_str(),
                       ("previous " + names.transient_variables[i]).c_str());
    for (size_t i = 0; i < names.full_reset_variables.size(); i++)
       model.CopyReplace( names.full_reset_variables[i].c_str(),
                       ("reset " + names.full_reset_variables[i]).c_str());
  } // end AdvanceTransientVariables

/** reset variables for transient pressure calculations */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::ResetVariables()
  {
    for (size_t i = 0; i < names.transient_variables.size(); i++)
       model.CopyReplace( ("previous " + names.transient_variables[i]).c_str(),
                       names.transient_variables[i].c_str());
  } // end ResetVariables

/** reset variables for transient calculations after thermal equilibration */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::FullReset()
  {
    for (size_t i = 0; i < names.transient_variables.size(); i++)
       model.CopyReplace( ("previous " + names.transient_variables[i]).c_str(),
                       names.transient_variables[i].c_str());
    for (size_t i = 0; i < names.full_reset_variables.size(); i++)
       model.CopyReplace( ("reset " + names.full_reset_variables[i]).c_str(),
                       names.full_reset_variables[i].c_str());
    model.Accept( pres_grad );
    transport.UpdateProjection();
  } // end FullReset

/** outer loop including advection and pressure diffion until mass-based time step criterion is met */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::AdvectionDiffusionLoops()
    {

      move_on = false;
      while( ! move_on )
       {

        PressureLoop();

		control_dt = transport.AdvectMassConserved( current_dt );
        
		if ( control_dt > current_dt )
         {
          cout << "\ncontrol_dt > current_dt adouble64er advection step" << endl;
          cin >> temp;
         }
        else if ( current_dt > control_dt )
         {
          ResetVariables();
          current_dt = control_dt*0.9;
         }
        else move_on = true;

      } // end move_on;

    } // end AdvectionDiffusionLoops

/** inner loop including pressure diffion until cfl-based time step criterion is met */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::PressureLoop()
    {
      pressure_loop = true;
      while( pressure_loop )
       {

         //********
	     // control volume finite element calculations
        // changes pressure

         // determine upwind nodes
        upwind_control.Reset();
        model.Accept( upwind_control );

        p_FE_SAMG.TimeIncrement( current_dt );
        model.Apply( p_FE_SAMG );

        // check for pressure below atmospheric
	     model.MinMaxOf( "fluid pressure", min_value, max_value );
        if ( (min_value < 0.1e6 || max_value > 1.e9) && brick_wall_limiter)
           {
               cerr << "\ncommented out." << endl;
//             model.Pass(pressure_limiter_transport);
           }
        else if ( !brick_wall_limiter && min_value < 0.0)
           {
            cout << "\n\n\nPressure below zero!!" << endl;
            cout << "min_value: " << min_value << ", max_value: " << max_value << endl;
            cout << "\ntimestep: " << timestep << endl;
            cin >> temp;
          }

         // calculate K * grad P
         model.Accept( pres_grad );

         // project pressure gradient onto facets
         transport.UpdateProjection();

         // calculate velocities for Finite Volume calculations
         // or CVFEM_Visitors and determine CFL criterion
         upwind_control.WithVelocity( true );
         model.Accept( upwind_control );

         model.MinMaxOf( "courant liquid", cfl_min_l, cfl_max_l );
         model.MinMaxOf( "courant vapor",  cfl_min_v, cfl_max_v );
       
         cfl_dt = std::min( cfl_min_l,cfl_min_v );

         if ( current_dt > cfl_dt)
          {
            ResetVariables();
            current_dt = cfl_dt*0.9;
          }
         else
           pressure_loop = false;

      }
    } // end PressureLoop

/** application of CVFEM_visitors for consistency check - currently not in use */
/*
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::ApplyCVFEM_Visitors()
    {
    
     model.InputUniformScalarValue(names.conduction_visitor_variables[1].c_str(),0.0);
     model.AssignBoundaryFlags(TOP,names.conduction_visitor_variables[1].c_str(), DIRICH );
//     model.AssignBoundaryFlags(LEdouble64,names.conduction_visitor_variables[1].c_str(), DIRICH );
//     model.AssignBoundaryFlags(RIGHT,names.conduction_visitor_variables[1].c_str(), DIRICH );

      //********
      // CVFEM_Visitor
      mass_visitor.SetTimeIncrement( current_dt );
      enthalpy_visitor.SetTimeIncrement( current_dt );
      conduction_visitor.SetTimeIncrement( current_dt );
      model.Accept( mass_visitor );
      model.Accept( enthalpy_visitor );

      model.Accept( conduction_visitor );
      
      enthalpy += dhc;

    
    } // end ApplyCVFEM_Visitors
*/

/** thermal quilibration between fluid and rock */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::FluidRockEquilibration()
    {

     //********
	  // equilibration
     equilibrator_properties.SetTimeIncrement( current_dt );
	  model.Accept( equilibrator_properties );
    
    } // end FluidRockEquilibration


/** performing consistency check between CVFEM_visitors and FV calculations - currently not in use */
/*template<size_t dim>
void CVFEM_PHX_Scheme<dim>::CheckForConsistency()
    {
    
     if (open_top)
       {
        mass     -= bfm;
        enthalpy -= bfe;
       }

     model.InputUniformScalarValue(names.diff_mass_variables[0].c_str(),0.0);
     model.InputUniformScalarValue(names.diff_enthalpy_variables[0].c_str(),0.0);

     diff_mass += mass;
     diff_mass -= fluid_density;
     diff_enthalpy += enthalpy;
     diff_enthalpy -= total_enthalpy;

	  model.MinMaxOf( "diff mass", min_value, max_value );
	  if (min_value < -1.e-6 || max_value > 1.e-6)
	    {
          cout << "\ndiff mass: " << endl;
          cout << "min_value: " << min_value << endl;
          cout << "max_value: " << max_value << endl;
          cin >> temp;	  
	    }
		*/
/*	  model.MinMaxOf( "diff enthalpy", min_value, max_value );
	  if (min_value < -1.e5 || max_value > 1.e5)
	    {
          cout << "\ndiff enthalpy: " << endl;
          cout << "min_value: " << min_value << endl;
          cout << "max_value: " << max_value << endl;
          cin >> temp;	  
       }*/
/*    } // end CheckForConsistency
*/    


/** modify calculations of temperature-dependent heat capacity of the rock */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::TemperatureDependentHeatCapacityRock( double64 cpr_min_ext, double64 t_min_ext,
                                               double64 cpr_max_ext, double64 t_max_ext )
  {

    equilibrator_properties.TemperatureDependentHeatCapacityRock(cpr_min_ext,t_min_ext,cpr_max_ext,t_max_ext);

  } // end TemperatureDependentHeatCapacityRock

/** switch on open top, specifying temperature, pressure and salinity of inflowing fluid */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::OpenBoundaries( double64 T_gradC, double64 p_Pa, double64 wt )
  {

    open_top = true;
    Brine  brine(T_gradC, p_Pa/1.0e5, Weight2XNaCl(wt));
    equilibrator_properties.WithOpenBoundaries( brine.Enthalpy(), wt );

  } // end OpenBoundaries

/** switch on open top, specifying salinity of inflowing fluid */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::OpenBoundaries( double64 wt )
  {

    open_top = true;
    equilibrator_properties.WithOpenBoundaries( wt );

  } // end OpenBoundaries

/** switch consostency check on or off */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::CheckConsistency( bool check )
  {
    check_consistency = check;
  } // end CheckConsistency

/** switch brick wall limiter for fluid pressure on or off */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::SetBrickWallLimiterTo( bool limit )
  {
    brick_wall_limiter = limit;
  } // end SetBrickWallLimiterTo

/** add further variables for FV calculations */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::AddAdvectionVariable( const char* balanced_variable,
                                                     const char* new_lhs_liquid, const char* new_rhs_liquid,
                                                     const char* new_lhs_vapor,  const char* new_rhs_vapor  )
  {
    names.AddAdvectionVariable(balanced_variable, new_lhs_liquid,new_rhs_liquid,
                               new_lhs_vapor,new_rhs_vapor);
    transport.AddAdvectionVariable(new_lhs_liquid,new_rhs_liquid,
                               new_lhs_vapor,new_rhs_vapor);
  } // end AddAdvectionVariable

/** modifying cfl criterion */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::Adjust_CFL_Criterion(double64 scale_factor, bool take_pore_velocity)
    {
      upwind_control.Adjust_CFL_Criterion(scale_factor,take_pore_velocity);
    } // end Adjust_CFL_Criterion


/** access function to transient fluxes */
template<size_t dim>
void CVFEM_PHX_Scheme<dim>::GetFacetFluxFromInsideNodeToOutsideNode( Element<dim>& e, unsigned int facet_idx,
                                                                     double64& flux_liquid, double64& flux_vapor )
   {
    flux_liquid = fv_transport_liquid.GetFacetFlux(e,facet_idx,0U);
    flux_vapor = fv_transport_vapor.GetFacetFlux(e,facet_idx,0U);
   }


template class CVFEM_PHX_Scheme<1U>;
template class CVFEM_PHX_Scheme<2U>;
template class CVFEM_PHX_Scheme<3U>;
 
} // end namespace csmp
