#include "TwoPhaseTransportPHX.h"
#include "Region.h"

using namespace std;

namespace csmp {


/** custom constructor */
template<size_t dim>
TwoPhaseTransportPHX<dim>::TwoPhaseTransportPHX( Model<dim>& model,
                          UpwindControlVisitor<dim>& upwind_visitor,
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid)
   : model_ref(model),
   fv_transport_vapor( fv_vapor ),
   fv_transport_liquid( fv_liquid ),
   UpwindVisitor( upwind_visitor ),
   cfl_dt(0.0),
   upper_shell_T(450.),
   lower_shell_T(350.)
{ 

    tffi_key  =  model.Database().StorageKey( "total fluid flux in" );
    tffo_key  =  model.Database().StorageKey( "total fluid flux out" );
    mlfi_key  =  model.Database().StorageKey( "magmatic liquid flux in" );
    mlfo_key  =  model.Database().StorageKey( "magmatic liquid flux out" );
    mvfi_key  =  model.Database().StorageKey( "magmatic vapor flux in" );
    mvfo_key  =  model.Database().StorageKey( "magmatic vapor flux out" );
    slfi_key  =  model.Database().StorageKey( "shell liquid flux in" );
    slfo_key  =  model.Database().StorageKey( "shell liquid flux out" );
    svfi_key  =  model.Database().StorageKey( "shell vapor flux in" );
    svfo_key  =  model.Database().StorageKey( "shell vapor flux out" );
    T_key     =  model.Database().StorageKey( "temperature" );

    cout <<"\n\nTwoPhaseTransportPHX<"<< typeid(double).name() <<","<< dim;
    cout <<">: Constructed successfully."<< endl;

    
  } // end constructor 

/** default destructor */
template<size_t dim>
TwoPhaseTransportPHX<dim>::~TwoPhaseTransportPHX()
 {
 }

/** modify maximum size of time step */
template<size_t dim>
void TwoPhaseTransportPHX<dim>::SetLargestTimeStep( const double& max_time_step )
 {
    fv_transport_liquid.SetMaximumTimeStep( max_time_step );
    fv_transport_vapor.SetMaximumTimeStep( max_time_step );
 }

/** main function to coordinate two-phase flow */
template<size_t dim>
double  TwoPhaseTransportPHX<dim>::AdvectMassConserved( const double& time_increment )
 {
 
   cfl_dt               = time_increment;

   // project velocities to face normals
   UpdateProjection();
    
   // determine facet fluxes and flux out
   // get suggestions for time step adjustments
   DetermineFacetFluxes();
 
   CheckDryFiniteVolumesAndAdjustTimestep( );
   
   // determine flux in and perform and write fluxes
   PerformFacetFluxes( );

   // track accumulated fluxes
   TrackFluxes();

   // return the time step used for the advection 
   cfl_dt *= time_step_factor;

   return cfl_dt;
 
 } // end AdvectMassConserved


/** calculating facet fluxes for both phases */
template<size_t dim>
void  TwoPhaseTransportPHX<dim>::DetermineFacetFluxes()
  {
    fv_transport_liquid.DetermineFacetFlux( cfl_dt, UpwindVisitor.UpwindMatrices( fv_transport_liquid.GetDensityKey()) );
    fv_transport_vapor.DetermineFacetFlux(  cfl_dt, UpwindVisitor.UpwindMatrices( fv_transport_vapor.GetDensityKey()) );  
  } // end DetermineLiquidAndVaporFluxes


/** checking mass-based time step criterion and adjusting time step if necessary */
template<size_t dim>
void  TwoPhaseTransportPHX<dim>::CheckDryFiniteVolumesAndAdjustTimestep()
  {
  
  time_step_factor = 1.0;

  unsigned int idx;
  double LHS, Outflow, temp_factor;
     
  // time step is only cut if the primary variable (liquid + vapor mass) is running dry
  for ( typename vector<Node<dim>*>::const_iterator 
       fvit=model_ref.Region("Model").NodesBegin(); 
       fvit!=model_ref.Region("Model").NodesEnd(); fvit++ )
       { 
         idx = (*fvit)->Idx();

		 LHS      = fv_transport_vapor.GetMainPropertyLHS(idx);
         Outflow  = fv_transport_vapor.GetFluxOut(idx);
         if (Outflow != 0.) temp_factor = LHS / Outflow;
         else temp_factor = 1.;

         if ( temp_factor < time_step_factor)
           {
             time_step_factor = temp_factor;
           }

         LHS     = fv_transport_liquid.GetMainPropertyLHS(idx);
         Outflow = fv_transport_liquid.GetFluxOut(idx);
         if (Outflow != 0.) temp_factor = LHS / Outflow;
         else temp_factor = 1.;

         if ( temp_factor < time_step_factor)
           {
             time_step_factor = temp_factor;
           }
           
       }	             
  } // end CheckDryFiniteVolumesAndAdjustTimestep

/** finalize FV calculations */
template<size_t dim>
void  TwoPhaseTransportPHX<dim>::PerformFacetFluxes( )
  {

    fv_transport_vapor.AdjustAndPerformFacetFlux(  time_step_factor );
    fv_transport_liquid.AdjustAndPerformFacetFlux( time_step_factor );

  } // end PerformFacetFluxes

/** active gravity component */
template<size_t dim>
void TwoPhaseTransportPHX<dim>::WithGravityComponentLiquidAndVapor()
   {
      fv_transport_liquid.WithGravityComponent();
      fv_transport_vapor.WithGravityComponent();
   }

/** update projection of velocity onto facet normal */
template<size_t dim>
void TwoPhaseTransportPHX<dim>::UpdateProjection( )
   {
      fv_transport_liquid.UpdateProjection(  );
      fv_transport_vapor.UpdateProjection(  );
   }

/** special function to track fluxes of magmatic fluids */
template<size_t dim>
void TwoPhaseTransportPHX<dim>::TrackFluxes( )
   {

  // CAREFUL Hard-coded index!!!!
  unsigned int idx_magmatic_mass(3);
//  unsigned int idx_magmatic_mass(0);

  unsigned int idx;
  for ( typename vector<Node<dim>*>::const_iterator
       fvit=model_ref.Region("Model").NodesBegin();
       fvit!=model_ref.Region("Model").NodesEnd(); fvit++ )
       {

       if ( (*fvit)->Status(  tffi_key ) != DIRICH )
       {
         idx = (*fvit)->Idx();
         (*fvit)->Read(tffi_key, tffi );
         (*fvit)->Read(tffo_key, tffo );
         (*fvit)->Read(mlfi_key, mlfi );
         (*fvit)->Read(mlfo_key, mlfo );
         (*fvit)->Read(mvfi_key, mvfi );
         (*fvit)->Read(mvfo_key, mvfo );
         (*fvit)->Read(slfi_key, slfi );
         (*fvit)->Read(slfo_key, slfo );
         (*fvit)->Read(svfi_key, svfi );
         (*fvit)->Read(svfo_key, svfo );
         (*fvit)->Read(T_key, temperature );

         tffi += fv_transport_liquid.GetFluxIn(idx);
         tffi += fv_transport_vapor.GetFluxIn(idx);
         tffo += fv_transport_liquid.GetFluxOut(idx);
         tffo += fv_transport_vapor.GetFluxOut(idx);

         mlfi += fv_transport_liquid.GetFluxIn(idx,idx_magmatic_mass);
         mlfo += fv_transport_liquid.GetFluxOut(idx,idx_magmatic_mass);
         mvfi += fv_transport_vapor.GetFluxIn(idx,idx_magmatic_mass);
         mvfo += fv_transport_vapor.GetFluxOut(idx,idx_magmatic_mass);

         if (temperature() < upper_shell_T && temperature() > lower_shell_T)
         {
         slfi += fv_transport_liquid.GetFluxIn(idx,idx_magmatic_mass);
         slfo += fv_transport_liquid.GetFluxOut(idx,idx_magmatic_mass);
         svfi += fv_transport_vapor.GetFluxIn(idx,idx_magmatic_mass);
         svfo += fv_transport_vapor.GetFluxOut(idx,idx_magmatic_mass);
         }

         (*fvit)->Store(tffi_key, tffi );
         (*fvit)->Store(tffo_key, tffo );
         (*fvit)->Store(mlfi_key, mlfi );
         (*fvit)->Store(mlfo_key, mlfo );
         (*fvit)->Store(mvfi_key, mvfi );
         (*fvit)->Store(mvfo_key, mvfo );
         (*fvit)->Store(slfi_key, slfi );
         (*fvit)->Store(slfo_key, slfo );
         (*fvit)->Store(svfi_key, svfi );
         (*fvit)->Store(svfo_key, svfo );
       }

       }

   }

/** adding further variables for FV calculations */
template<size_t dim>
void TwoPhaseTransportPHX<dim>::AddAdvectionVariable( const char* new_lhs_liquid, const char* new_rhs_liquid,
                                                        const char* new_lhs_vapor,  const char* new_rhs_vapor )
  {
    fv_transport_liquid.AddAdvectionVariable(new_lhs_liquid,new_rhs_liquid);
    fv_transport_vapor.AddAdvectionVariable(new_lhs_vapor,new_rhs_vapor);
  }


template class TwoPhaseTransportPHX<1U>;
template class TwoPhaseTransportPHX<2U>;
template class TwoPhaseTransportPHX<3U>;
 
} // end namespace csmp

