// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef EXPLICIT_NODE_CENTERED_FV_TRANSPORT_H
#define EXPLICIT_NODE_CENTERED_FV_TRANSPORT_H

#include "finiteVolumeAuxiliaryFunctions.h"
#include "NodeCenteredFiniteVolumeTransport.h"
#include "ExplicitStencilProcessor.h"
#include "ExplicitNodeCenteredFiniteVolumeTransport.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "STL_utilities.h"
#include "Point.h"
#include "Exception.h"
#if defined(_OPENMP )
#include "omp.h"
#endif


namespace csmp {

template<uint32_t dim,template<uint32_t> class STP>
class ExplicitNodeCenteredFiniteVolumeTransport : public NodeCenteredFiniteVolumeTransport<dim> {
  public:
    /// single phase solute advection-only constructor for a subregion identified as a group
    ExplicitNodeCenteredFiniteVolumeTransport( const char* region,
                                               Model<dim>&,
                                               const char* porosity,
                                               const char* advected_prop,
                                               const char* transp_velocity,
                                               const char* nodal_source,
                                               bool second_order_accuracy);
    
   
    /// single phase solute advection-dispersion constructor for a subregion identified as a group
    ExplicitNodeCenteredFiniteVolumeTransport( const char* region,
                                               Model<dim>&,
                                               const char* porosity,
                                               const char* diffusivity,
                                               const char* advected_prop,
                                               const char* transp_velocity,
                                               const char* nodal_source,
                                               bool second_order_accuracy);

    ExplicitNodeCenteredFiniteVolumeTransport( const char* region,
                                               Model<dim>&,
                                               const char* porosity,
                                               const char* advected_prop_lhs,
                                               const char* advected_prop_rhs,
                                               const char* transp_velocity,
                                               const char* nodal_source,
                                               bool second_order_accurate,
                                               bool second_order_in_time,
                                               const char* elmt_thickness_attribute=NULL );
                                       
    virtual ~ExplicitNodeCenteredFiniteVolumeTransport();
    
    /// single-phase transport (@todo NEEDS TO BECOME A VIRTUAL FUNCTION TEMPLATE)
    virtual double AdvectVariable( double time_interval,
                                   double cfl_multiplication_factor, ///< ideally 0.1
                                   bool apply_flux_balance_correction, ///< usually desirable
                                   bool update_pore_volumes );         ///< normally not needed

    /// single-phase passive advection, does NOT return courant increment, single timestep calculation
    /// no checks are made for courant condition.  Assumes external checks.
    virtual void AdvectVariableSingleStep( double time_increment,
                                           bool apply_flux_balance_correction,
                                           bool update_pore_volumes);

    /// results are stored back wherever flags are not DIRICH
    double OutputResults(const PropertyDatabase<dim>&,
                         const csmp::Index& adv_key,
                         bool  show_range,
                         uint32_t var_comp_nr=0u ) const;

  protected:
    virtual void AccumulateFluxUpwindProducts();
    /// as in base class but for a result vector
    virtual void AssignFluxBoundaryConditions( uint32_t var_comp_nr=0 );
    
    /// passive advection: gives transported variable at end of time_interval
    virtual void ComposeSolution( double time_interval, uint32_t var_comp_nr=0u );
    
    STP<dim>             stencil_;
    std::vector<double>  RESULT;

#if defined(_OPENMP )
    std::vector<STP<dim>* >  thread_stencil_processor_;
#endif

  protected:  
    /// explicit conservative
    void AdvectVariable1stOrder( double time_interval, bool output_result_range );
    
  private:
    /// explicit conservative
    void AdvectAndDiffuseVariable1stOrder( double time_interval, bool output_result_range );

    /// 2nd-order explicit conservative, using slope limiting
    void AdvectVariable2ndOrder( double time_increment, bool output_result_range );

    /// 2nd-order explicit conservative, using slope limiting
    void AdvectAndDiffuseVariable2ndOrder( double time_increment, bool output_result_range );
};


/**
 
@class ExplicitNodeCenteredFiniteVolumeTransport  ExplicitNodeCenteredFiniteVolumeTransport 
"generic_node_centered_finite_volumes\ExplicitNodeCenteredFiniteVolumeTransport.h"
@author S.K. Matthaei
@date 2002

@section motivation Motivation

The ExplicitNodeCenteredFiniteVolumeTransport module is provided to allow users
to perform passive advection and advection of a single non-wetting 
phase in two-phase flow simulations. It represents the fastest but
CFL dependent transport scheme for arbitrary FE meshes including 
different FE types and mixes of volume and surface elements.  

The scheme is also very tolerant in that no specific boundary conditions 
need to be applied. The value of the advected variable at inflow 
boundaries is not modified and governs what gets transported into
the model domain.  

 
@section design Design Intent

The scheme has been designed to capitalize on the new generic polytype
CVFE technology that is part of CSMP.  
 
 
@section applicability Applicability

Use the scheme if the transport velocity field is rapidly changing or
invalidated by the transport itself. In this case the scheme is faster
that the implicit scheme embodied in the FiniteVolumeTransport module
which starts to pay off at about CFL=20 or greater.  

Do not use this scheme when the CFL condition based transport increment
leads global advection steps that are smaller than about 1/1,000 of the 
model length in the direction of the flow. In this case the computations
would be prohibitively slow.
 
 
@section structure Structure

The explicit advection scheme consists out of the three modules:  
 
ExplicitNodeCenteredFiniteVolumeTransport
ExplicitFiniteVolumeAlgorithm
ExplicitFiniteVolumeProcessor
  

The user will only interact with the ExplicitNodeCenteredFiniteVolumeTransport module,
but the ExplicitFiniteVolumeAlgorithm is used to manage the accumulation
and output of the solution in each transport step, and the 
ExplicitFiniteVolumeProcessor class facilitates the generation of the 
finite-volume stencils from each finite element, so that the accumulation
can proceed on an element by element basis and no FV mesh needs to be 
stored.
 
 
@section collaboration Collaboration

The scheme collaborates with the finite elements of the current model.
 

@section consequences Consequences

This implementation of the scheme cannot be restricted to Regions, but 
always needs to be applied to the entire model, i.e. Model.  

Currently it also only works for 3D models. Surface computations can
still be done if the mesh consists entirely of surface elements in
3D space. 

 
@section implementation Implementation

The flux-balances are computed for all finite volumes including those 
which are located on the boundaries of the Region or Model to 
which the transport scheme is applied. They are compensated for also
on all domain boundaries.  

 
@section examples Application Examples

This is an example of how the second order accurate scheme is applied:

@code
 ExplicitNodeCenteredFiniteVolumeTransport  explicit_advector( model, "porosity", "diffusivity", 
                                                  "concentration", "velocity", second_order );

 cout <<"\n\nadvectVariableExplicit: Configuring TRANSPORT simulation: ";
 if ( second_order ) cout <<" IMPES: SECOND ORDER SCHEME."<< endl;
 else                cout <<" IMPES: FIRST ORDER SCHEME."<< endl;
 cout <<"\nThe grid Courant number is "<< explicit_advector.CourantIncrement( sg ) << endl;
 cout <<"\nEnter advection time: ";
 double time_interval;
 cin >> time_interval;

 cout <<"\n\tMeasuring the time required to solve the advection problem."<< endl;
 clock_t ticks = clock();
 explicit_advector.AdvectVariable( model, time_interval );
 ticks = clock() - ticks;
 cout <<"\n\n\tCPU clock ticks used for advection step: "<< ticks << endl;
@endcode

*/


/// single-phase passive advection (group-restricted)
template<uint32_t dim,template<uint32_t> class STP>
ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::ExplicitNodeCenteredFiniteVolumeTransport(  
                                                                            const char* group,
                                                                            Model<dim>& sg, 
                                                                            const char* porosity,
                                                                            const char* advected_prop_t1,
                                                                            const char* transp_velocity,
                                                                            const char* nodal_source,
                                                                            bool second_order_accurate)
 : NodeCenteredFiniteVolumeTransport<dim>( group, sg, porosity, 
                                           advected_prop_t1, transp_velocity, 
                                           nodal_source, second_order_accurate, false),

   stencil_(this->adv1_key_, this->vel_key_),
   RESULT(this->gref_.Nodes())
 {
#if defined(_OPENMP )
    this->thread_stencil_processor_.reserve(omp_get_max_threads());
    for (size_t tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->thread_stencil_processor_.push_back(new STP<dim>(this->adv1_key_, this->vel_key_));
#endif
    if (this->Verbose())
        std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim <<">(constructor): Constructed successfully."<< std::endl;
 } // end constructor (solute transport-only)


/// single-phase passive advection (group-restricted)
template<uint32_t dim,template<uint32_t> class STP>
ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::ExplicitNodeCenteredFiniteVolumeTransport(  
                                                                            const char* group,
                                                                            Model<dim>& sg, 
                                                                            const char* porosity,
                                                                            const char* advected_prop_lhs,
                                                                            const char* advected_prop_rhs,
                                                                            const char* transp_velocity,
                                                                            const char* nodal_source,
                                                                            bool second_order_accurate,
                                                                            bool second_order_in_time,
                                                                            const char* elmt_thickness_attribute)
 : NodeCenteredFiniteVolumeTransport<dim>( group, sg, porosity, 
                                           advected_prop_lhs, transp_velocity,
                                           nodal_source, second_order_accurate,
                                           second_order_in_time, elmt_thickness_attribute ), // introduced thickness attribute.
  // for ExplicitMassBasedTransport add thickness here if necessary
  // also needs to be added in FacetArea calculation

   stencil_(this->adv1_key_, sg.Database().StorageKey(advected_prop_rhs), this->vel_key_),
   RESULT(this->gref_.Nodes())
 {
#if defined(_OPENMP )
    this->thread_stencil_processor_.reserve(omp_get_max_threads());
    for (size_t tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->thread_stencil_processor_.push_back( new STP<dim>(this->adv1_key_, sg.Database().StorageKey(advected_prop_rhs), this->vel_key_));
#endif
    if (this->Verbose())
        std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim <<">(constructor): Constructed successfully."<< std::endl;
    
 } // end constructor (solute transport-only)


/**
 
Constructor for the single-phase passive advection scheme for slightly
compressible flow without fluid density and viscosity variations. The
ExplicitNodeCenteredFiniteVolumeTransport can be constructed for first- and second
order calculations.  

Diffusion is only taken into account if second-order is chosen, because
the first-order scheme already is so diffusive that adding an extra
diffusion term makes little sense.   

@section arguments Input Arguments 

The ExplicitNodeCenteredFiniteVolumeTransport gains access via reference to the 
model to which it is applied to. The user needs to specify the piecewise
constant element variables porosity, diffusivity and transport velocity
as must be defined in the CSP_variables.txt file. The transported 
variable must be nodal since this explicit transport scheme is based
on node-centered finite volumes.  

@section implementation Implementation

The constructor creates minor storage array for repetitively used
variables like the finite (pore) volumes, the flux balance for each
FV cell and a vector<double> of FV_parameters called FV_stencil_data that holds
the FV sector volumes, facet areas and facet normal fluxes which are 
always updated before the first advection step, when the CFL criterion
is computed. This vector<double> will become redundant when the generic FV scheme
has been optimized for speed.  

@section application Application

Use ExplicitNodeCenteredFiniteVolumeTransport as constructed with this particular
constructor for the passive advection of tracers in a rapidly evolving
flow field that necessitates frequent updates of the velocity fields.  

@section messages Messages 

The memory that is consumed during the construction process is reported.  
*/
template<uint32_t dim,template<uint32_t> class STP>
ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::ExplicitNodeCenteredFiniteVolumeTransport(  
                                                                            const char* group,
                                                                            Model<dim>& sg,
                                                                            const char* porosity,
                                                                            const char* diffusivity,
                                                                            const char* advected_prop_t1,
                                                                            const char* transp_velocity,
                                                                            const char* nodal_source,
                                                                            bool second_order_accurate)
 : NodeCenteredFiniteVolumeTransport<dim>( group, sg, porosity, diffusivity, 
                                           advected_prop_t1, transp_velocity,
                                           nodal_source, second_order_accurate, false ),
   stencil_(this->adv1_key_, this->vel_key_, this->diff_key_),
   RESULT(this->gref_.Nodes())
 {
 
#if defined(_OPENMP )
    this->thread_stencil_processor_.reserve(omp_get_max_threads());
    for (size_t tid = 0 ; tid < omp_get_max_threads();tid ++)
        this->thread_stencil_processor_.push_back(new STP<dim>(this->adv1_key_, this->vel_key_, this->diff_key_));
#endif

    if (this->Verbose()){
        std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<";
        std::cout << dim <<">(constructor): Constructed successfully."<< std::endl;
    }
    
 } // end constructor (solute transport)







/** Frees all the allocated memory.
*/
template<uint32_t dim,template<uint32_t> class STP>
ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::~ExplicitNodeCenteredFiniteVolumeTransport()
 {
#if defined(_OPENMP )
    for ( size_t i = 0 ; i < omp_get_max_threads(); i++)
         delete(thread_stencil_processor_[i]);
#endif
 }


/*
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AssignFluxBoundaryConditions(const size_t var_comp_nr)
{
    double        inflow, flux_balance;
    const double  zero(0.);
    
    if (this->adv1_key_.type ==SCALAR){
        // loop over the boundary cells and adjust fluxes
        for ( typename std::vector<Node<dim>*>::const_iterator
              nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ ){
            // if the flux balance cannot be evaluated because the node sits at the model boundary
            if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) )
            {
                // inflow or outflow compensation
                if ( inflow != zero ) RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * -inflow;
            }
            // if we have a boundary flux balance because we are dealing with interior boundaries
            else if ( flux_balance != zero ) // counter balancing
                RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * -flux_balance;
        }
    }
    else{
        // loop over the boundary cells and adjust fluxes
        if (this->adv1_key_.type == ARRAY){
            for ( typename std::vector<Node<dim>*>::const_iterator
                  nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ ){
                // if the flux balance cannot be evaluated because the node sits at the model boundary
                ArrayVariable av;
                (*nit)->Read( this->adv1_key_ ,av);

                if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) ){
                    // inflow or outflow compensation
                    if ( inflow != zero ) RESULT[ (*nit)->Idx() ] += av(var_comp_nr) * -inflow;
                }
                // if we have a boundary flux balance because we are dealing with interior boundaries
                else if ( flux_balance != zero ) // counter balancing
                    RESULT[ (*nit)->Idx() ] += av(var_comp_nr) * -flux_balance;
            }
        } else if (this->adv1_key_.type == FLAGGEDARRAY){
            for ( typename std::vector<Node<dim>*>::const_iterator
                  nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ ){
                // if the flux balance cannot be evaluated because the node sits at the model boundary
                FlaggedArrayVariable fav;
                (*nit)->Read( this->adv1_key_ ,fav);

                if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) ){
                    // inflow or outflow compensation
                    if ( inflow != zero ) RESULT[ (*nit)->Idx() ] += fav(var_comp_nr) * -inflow;
                }
                // if we have a boundary flux balance because we are dealing with interior boundaries
                else if ( flux_balance != zero ) // counter balancing
                    RESULT[ (*nit)->Idx() ] += fav(var_comp_nr) * -flux_balance;
            }
        }
    }
    
} // end AssignFluxBoundaryConditions
*/


/**
         Luat's fix 4/6/2020
*/
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AssignFluxBoundaryConditions( uint32_t  )
  {
      double        inflow, flux_balance;
      const double  zero(0.);

      if (this->adv1_key_.type != SCALAR)
      {
          throw "advective is not scalar variable";
      }
      
      // loop over the boundary cells and adjust fluxes
      for ( auto nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); ++nit)
      {
          this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) ;
          if(inflow != zero)
          {
              RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * -inflow;
          }
          else
          {
              RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * flux_balance;
          }
      }
  } // end AssignFluxBoundaryConditions


   
/**
 
Single-phase transport in a non-divergent flow. The method is either
first- or second-order accurate depending on the boolean flag that 
was set when the 'ExplicitNodeCenteredFiniteVolumeTransport' object was constructed.

@section arguments Input Arguments 

The method acts on the supplied Model object and advects the 
transport variable in the current velocity field for the time interval
that the user specifies as second method argument.  

@return The advected transport variable is written back to the Model
object.  

@section implementation Implementation

The CFL_FACTOR is a multiplier that represents the maximum factor by
which the diameter of the idealized spherical FV can be multiplied
while the transport scheme is stable.  

@section messages Messages 

For each completed advection increment the method writes a dot on the
screen as a progress monitor.
 
@todo (3) Needs timestep reduction and DIRICHLET inflow condition 

*/
template<uint32_t dim,template<uint32_t> class STP>
double ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectVariable( double time_interval,
                                                                             double cfl_multiplication_factor,
                                                                             bool apply_flux_balance_correction,
                                                                             bool update_pore_volumes )
 {
    this->gref_.RenumberNodes();

    // 1. update the sector volumes and FVPOREVOL if so required (as in heat transport for instance)
    // ---------------------------------------------------------------------------------------------
    if ( update_pore_volumes ) {
         const bool multiply_with_thickness_attribute = (this->thi_key_ != csmp::Index() ) ? true : false;
         this->InitializeSectorPoreVolumeData( multiply_with_thickness_attribute );
         this->InitializeArraysForFirstOrderMethod();
      }

    const double courant_increment( this->AnisotropicCourantIncrement() );

    if ( cfl_multiplication_factor > 0.5 ) {
      if (this->Verbose())
        {
            std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim <<">::AdvectVariable: ";
            std::cout <<" In explicit mode, CFL multiplication factor must be less or equal to one(1/2). resetting..."<< std::endl;
        }   
      cfl_multiplication_factor = 0.5;
      }
    if ( !apply_flux_balance_correction ) {
      if (this->Verbose())
        {
            std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim <<">::AdvectVariable: ";
            std::cout <<" Flux balance correction is applied automatically."<< std::endl;
        }
      }

    // 2. compute the CFL condition to identify value for overstepping
    // --------------------------------------------------------------
    double time(0.), time_increment = std::min( courant_increment * cfl_multiplication_factor, time_interval ); 
    bool output_result_range(false);
    bool finish(false);
      
    // 3. compute solution
    // -------------------
    if (this->Verbose()) std::cout <<"\n\n\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim;
    if (this->Verbose()) std::cout <<">::AdvectVariable: Advecting transport variable";
    while ( time <= time_interval  and !finish )
      {
         if ( (time_interval - time) <= time_increment ) {
              time_increment = time_interval - time;
              finish = true;
              if (this->Verbose())
                output_result_range = true;
           }

         if ( !this->SecondOrderInSpace() )
             if ( this->diff_key_ != csmp::Index() )
                 // explicit 1st-order solution of advection-diffusion equation
                 AdvectAndDiffuseVariable1stOrder( time_increment, output_result_range );
             else
               {
                 // explicit 1st-order solution of advection equation
                 AdvectVariable1stOrder( time_increment, output_result_range );
               }    

         else {
             // explicit 2nd-order solution
             if ( this->diff_key_ != csmp::Index() )
                 // explicit 2nd-order solution of advection-diffusion equation
                 AdvectAndDiffuseVariable2ndOrder( time_increment, output_result_range );
             else
                 // explicit 2nd-order solution of advection equation
                 AdvectVariable2ndOrder( time_increment, output_result_range );
           }

         time += time_increment;
      }
                                           
    return courant_increment;

 } // end AdvectVariable

/**
  @author: Julian E. Mindel (23-09-2013)

This method calculates a single timestep of transport, assumming nothing about the time interval
it is "fed".  It will advect using that time interval, assuming the user has externally determined
that this is the correct time interval length.  In contrast with AdvectVariable(, this method does
not break the time interval into sub-parts and guarantee stability. This also means that the courant
increment is not checked for in this method (and should be, externally, of course).

 */
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectVariableSingleStep( double time_increment,
                                                                                   bool apply_flux_balance_correction,
                                                                                   bool update_pore_volumes )
 {

    bool output_result_range(false);
#ifdef DEBUG_ExplicitNodeCenteredFiniteVolumeTransport
    output_result_range=true;
#endif
    this->gref_.RenumberNodes();

    // 1. update the sector volumes and FVPOREVOL if so required (as in heat transport for instance)
    // ---------------------------------------------------------------------------------------------
    if ( update_pore_volumes ) {
        const bool multiply_with_thickness_attribute = (this->thi_key_ != csmp::Index() ) ? true : false;
        this->InitializeSectorPoreVolumeData( multiply_with_thickness_attribute );
        this->InitializeArraysForFirstOrderMethod();
    }

    if ( !apply_flux_balance_correction ) {
      if (this->Verbose())
        {
            std::cout <<"\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim <<">::AdvectVariable: ";
            std::cout <<" Flux balance correction is applied automatically."<< std::endl;
        }
    }

    // 2. compute solution
    // -------------------
    if (this->Verbose()) std::cout <<"\n\n\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim;
    if (this->Verbose()) std::cout <<">::AdvectVariableSingleStep: Advecting transport variable";

    if ( !this->SecondOrderInSpace() )
        if ( this->diff_key_ != csmp::Index() )
            // explicit 1st-order solution of advection-diffusion equation
            AdvectAndDiffuseVariable1stOrder( time_increment, output_result_range );
        else
        {
            // explicit 1st-order solution of advection equation
            AdvectVariable1stOrder( time_increment, output_result_range );
        }

    else {
        // explicit 2nd-order solution
        if ( this->diff_key_ != csmp::Index() )
            // explicit 2nd-order solution of advection-diffusion equation
            AdvectAndDiffuseVariable2ndOrder( time_increment, output_result_range );
        else
            // explicit 2nd-order solution of advection equation
            AdvectVariable2ndOrder( time_increment, output_result_range );
    }

} // end AdvectVariableSingleStep

   






/**

Explicit first-order advection algorithm of a passive conservative tracer.
AdvectVariableFirstOrder() uses the upstream 
values of the advected variable to compose the solution in the result
vector<double> managed by the ExplicitFiniteVolumeAlgorithm which performs 
a range check and writes the RESULT back to the Model.  

@section arguments Input Arguments 

Reference to the current model and the ExplicitFiniteVolumeAlgorithm as
well as the time_increment for the transport step. This increment must 
fulfill the CFL criterion.  

@section implementation Implementation

The method uses the ExplicitFiniteVolumeProcessor object in the 
accumulation of the finite-volume stencils, when it loops over the 
elements.  

@section application Application

Since it costs almost the same to compute the second-order accurate 
explicit solution, there is not much rational to use this method.  

*/
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectVariable1stOrder( 
                                                            double time_increment,
                                                            bool output_result_range )
{
    for ( uint32_t ncom{0u}; ncom< this->var_ncomponents_; ncom++){
        if ( this->Verbose() ) std::cout<<" Advecting component (ENCFVT): "<< ncom <<std::endl;
        std::fill( RESULT.begin(), RESULT.end(), static_cast<double>(0.) );
        std::vector<FV_Parameter>::const_iterator  fvt = this->STENCIL_DATA.begin();

        for ( auto eit=this->gref_.CellsBegin();
              eit!=this->gref_.CellsEnd(); eit++, fvt++ )
        {
            stencil_.eidx_ = (*eit)->Idx();

            // getting all necessary data for the construction of the result vector from the element
            stencil_.InitializeFirstOrder( (*fvt), *(*eit) ,this->adv1_key_.type,ncom);

            // starting accumulation with in and out fluxes which are temporarily stored in result vector
            AccumulateFluxUpwindProducts();

        } // end of accumulation

        // compensate for the inflow and the outflow boundaries
        AssignFluxBoundaryConditions(ncom);

        // 4. Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
        //     - a single loop over all nodes is required
        //     - this steps also considers in and outflow of finite volume cells
        //     - (distributed) sources and sinks will be considered in the future
        ComposeSolution( time_increment ,ncom);


        // 5. Saving the computed new saturations at the FV centers
        OutputResults( this->pref_, this->adv1_key_, output_result_range ,ncom);
    }


} // end AdvectVariable1stOrder








/**
 
Explicit second-order advection algorithm of a passive conservative tracer.
AdvectVariableSecondOrder() uses the slope-limited values of the advected 
variable which were interpolated to the FV facet integration points in order 
to compose the solution into the result vector<double> as managed by the 
ExplicitFiniteVolumeAlgorithm. The latter performs a range check before 
the RESULT are written back to the Model.  

@section arguments Input Arguments 

Reference to the current model and the ExplicitFiniteVolumeAlgorithm as
well as the time_increment for the transport step. This increment must 
fulfill the CFL criterion.  

@section implementation Implementation

The method uses the ExplicitFiniteVolumeProcessor object in the 
accumulation of the finite-volume stencils, when it loops over the 
elements.  

@section application Application

This method is speediest way to get a second-order accurate solution
to the transport equation. However, the CFL constraint applies. Thus
one should use the implicit approach rather than this method, when the
velocity field is not rapidly changing.  
 */
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectVariable2ndOrder( double time_increment,
                                                                                 bool output_result_range )
 {
    std::fill( RESULT.begin(), RESULT.end(), static_cast<double>(0.) );
    
    // 0. diffusion is taken into account if the diffusion key in the base class is initialized
    //const bool with_diffusion( (this->dif_key_ == csmp::Index()) ? false : true );

    // 1. finding smin/smax of previous solution in the neighborhood of each FV
    this->MinMaxAdvectedProperty();

    std::vector<FV_Parameter>::const_iterator  fvt(this->STENCIL_DATA.begin());

    if(!this->with_lsmgrad_limiter_)
        for ( auto eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++, fvt++ )
             // 2.0 initializing the stencil array index (needed for the limiter function)
             stencil_.AccumulateExplicitAdvectionSolution2( this->SMINMAX,
                                                           (*fvt), *(*eit), RESULT );
    else{
        this->grad_advprop_limiter_->CalculateGenericNodalGradient();
        this->grad_advprop_limiter_->CalculateSlopeLimiter( this->gref_, this->SMINMAX );
        for ( auto eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++, fvt++ )
             // 2.0 initializing the stencil array index (needed for the limiter function)
             stencil_.AccumulateExplicitAdvectionSolution2( this->SMINMAX,
                                                           (*fvt), *(*eit), RESULT, this->mass_center_key_,
                                                           this->grad_advprop_key_,this->grad_advprop_limiter_key_ );
    }


    // 2.5 compensate for the inflow and the outflow boundaries
    AssignFluxBoundaryConditions();

    // 3.0 Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
    //     - a single loop over all nodes is required
    //     - this steps also considers in and outflow of finite volume cells
    //     - (distributed) sources and sinks will be considered in the future
    ComposeSolution( time_increment );

    // 3.1 Saving the computed new saturations at the FV centers
    OutputResults( this->pref_, this->adv1_key_,output_result_range );
 	 
 	// indicating that something has been done
    //std::cout <<".";
    //std::cout.flush();

 } // end AdvectVariable2ndOrder	                                                                     




template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectAndDiffuseVariable1stOrder( double time_increment,
                                                                                            bool output_result_range )
 {
    for ( uint32_t ncom{0u}; ncom<this->var_ncomponents_; ncom++ ){
        if (this->Verbose() ) std::cout<<" Advecting and diffusing component (ENCFVT): "<<ncom<<std::endl;
        std::fill( RESULT.begin(), RESULT.end(), static_cast<double>(0.) );
        // 0. diffusion is taken into account if the diffusion key in the base class is initialized
        //const bool with_diffusion( (this->dif_key_ == csmp::Index()) ? false : true );
        std::vector<FV_Parameter>::const_iterator  fvt(this->STENCIL_DATA.begin());

        for ( auto eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++, fvt++ )
            // 1.0 initializing the stencil array index (needed for the limiter function)
            stencil_.AccumulateExplicitAdvectionDiffusionSolution1( (*fvt), *(*eit), RESULT, this->adv1_key_.type, ncom );

        // 2.0 compensate for the inflow and the outflow boundaries
        AssignFluxBoundaryConditions(ncom);

        // 3.0 Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
        //     - a single loop over all nodes is required
        //     - this steps also considers in and outflow of finite volume cells
        //     - (distributed) sources and sinks will be considered in the future
        ComposeSolution( time_increment, ncom );

        // 3.1 Saving the computed new saturations at the FV centers
        OutputResults( this->pref_, this->adv1_key_,output_result_range, ncom );
    }

 } // end AdvectAndDiffuseVariable1stOrder


template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AdvectAndDiffuseVariable2ndOrder( double time_increment,
	                                                                                         bool output_result_range )
 {
    std::fill( RESULT.begin(), RESULT.end(), static_cast<double>(0.) );
    
    // 0. diffusion is taken into account if the diffusion key in the base class is initialized
    //const bool with_diffusion( (this->dif_key_ == csmp::Index()) ? false : true );

    // 1. finding smin/smax of previous solution in the neighborhood of each FV
    this->MinMaxAdvectedProperty();

    std::vector<FV_Parameter>::const_iterator  fvt(this->STENCIL_DATA.begin());

    if(!this->with_lsmgrad_limiter_)
        for ( typename std::vector<Element<dim>*>::const_iterator
              eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++, fvt++ )
             // 1.1 initializing the stencil array index (needed for the limiter function)
             stencil_.AccumulateExplicitAdvectionDiffusionSolution2( this->SMINMAX,
                                                                    (*fvt), *(*eit), RESULT );
    else{
        this->grad_advprop_limiter_->CalculateGenericNodalGradient();
        this->grad_advprop_limiter_->CalculateSlopeLimiter( this->gref_,this->SMINMAX );
        for ( typename std::vector<Element<dim>*>::const_iterator
              eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++, fvt++ )
             // 1.1 initializing the stencil array index (needed for the limiter function)
             stencil_.AccumulateExplicitAdvectionDiffusionSolution2( this->SMINMAX,
                                                                    (*fvt), *(*eit), RESULT, this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);
    }

    // 2.0 compensate for the inflow and the outflow boundaries
    AssignFluxBoundaryConditions();

    // 3.0 Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
    //     - a single loop over all nodes is required
    //     - this steps also considers in and outflow of finite volume cells
    //     - (distributed) sources and sinks will be considered in the future
    ComposeSolution( time_increment );

    // 3.1 Saving the computed new saturations at the FV centers
    OutputResults( this->pref_, this->adv1_key_,output_result_range );
 	 
    // indicating that something has been done
    //std::cout <<".";
    //cout.flush();

 } // end AdvectAndDiffuseVariable2ndOrder



/**
 
Computes the values of the transported variable at the new time-level 
(passive advection case). Fluid sources or sinks as well as in or out-
fluxes at Model or region boundaries as well as poroelastic sources 
are also taken care off by consideration of FLUX_BALANCE vector.  

@section arguments Input Arguments 

The method needs access to the current transport model and two vectors
which store the pore volume of each finite volume cell and the flux
balance of that cell, respectively. The method further requires a
csmp::Index accessor to the transport variable in order to read 
current nodal values. Also, the current time-increment for
the advection step needs to be specified in order to compose the 
solution. This time increment must satisfy the CFL condition, because 
this is an explicit transport scheme.  

@section implementation Implementation


@section application Application

For the passive advection of tracers.  
*/
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::ComposeSolution( double time_interval, uint32_t var_comp_nr )
{
    if (this->adv1_key_.type == SCALAR){
        for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
        {
            // reading the transported variable value
            double solution = this->gref_.N(nidx)->Read( this->adv1_key_ );
            // adding potential (volumetric) source or sink terms due to a divergence of the flow (+ dt sum_j^e 1/3 V_e q_j)
            if ( this->FLUX_BALANCE[nidx] != static_cast<double>(0.) )
                RESULT[nidx] += solution * -this->FLUX_BALANCE[nidx];
            // subtracting the flux time-interval product
            RESULT[nidx] = solution - (time_interval / this->FVPOREVOL[nidx]) * RESULT[nidx];
            // adding externally assigned source or sink terms
            RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );

        }
    }
    else if (this->adv1_key_.type == ARRAY){
        ArrayVariable av;
        if (this->src_key_.type == SCALAR)
        {
          for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
          {
              // reading the transported variable value
              this->gref_.N(nidx)->Read( this->adv1_key_ ,av);
              double solution = av(var_comp_nr);
              // adding potential (volumetric) source or sink terms due to a divergence of the flow (+ dt sum_j^e 1/3 V_e q_j)
              if ( this->FLUX_BALANCE[nidx] != static_cast<double>(0.) )
                  RESULT[nidx] += solution * -this->FLUX_BALANCE[nidx];
              // subtracting the flux time-interval product
              RESULT[nidx] = solution - (time_interval / this->FVPOREVOL[nidx]) * RESULT[nidx];
              // adding externally assigned source or sink terms
              RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );
          }
        }
        else if (this->src_key_.type == ARRAY)
        {
          ArrayVariable src_av;
          double source;
          for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
          {
              // reading the transported variable value
              this->gref_.N(nidx)->Read( this->adv1_key_ ,av);
              double solution = av(var_comp_nr);
              this->gref_.N(nidx)->Read( this->src_key_ , src_av);
              source = src_av[var_comp_nr];  
              // adding potential (volumetric) source or sink terms due to a divergence of the flow (+ dt sum_j^e 1/3 V_e q_j)
              if ( this->FLUX_BALANCE[nidx] != static_cast<double>(0.) )
                  RESULT[nidx] += solution * -this->FLUX_BALANCE[nidx];
              // subtracting the flux time-interval product
              RESULT[nidx] = solution - (time_interval / this->FVPOREVOL[nidx]) * RESULT[nidx];
              // adding externally assigned source or sink terms
              RESULT[nidx] += time_interval * source;
          }
        }
        else
        {
        std::cout <<"ExplicitNodeCenteredFiniteVolumeTransport<dim>::ComposeSolution(): source variable type not implemented"<< std::endl;
        }
    }
    else if (this->adv1_key_.type == FLAGGEDARRAY){
        FlaggedArrayVariable fav;
        for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
        {
            // reading the transported variable value
            this->gref_.N(nidx)->Read( this->adv1_key_ ,fav);
            double solution = fav(var_comp_nr);
            // adding potential (volumetric) source or sink terms due to a divergence of the flow (+ dt sum_j^e 1/3 V_e q_j)
            if ( this->FLUX_BALANCE[nidx] != static_cast<double>(0.) )
                RESULT[nidx] += solution * -this->FLUX_BALANCE[nidx];
            // subtracting the flux time-interval product
            RESULT[nidx] = solution - (time_interval / this->FVPOREVOL[nidx]) * RESULT[nidx];
            // adding externally assigned source or sink terms
            RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );

        }
    }

} // end ComposeSolution (passive advection case)

/**
    Accumulates facet fluxes (volume * saturation) coming into the control volumes into the result vector.
*/
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::AccumulateFluxUpwindProducts()
{
   // for all finite-volume facets
   for ( uint32_t i{0u}; i<this->gref_.E(stencil_.eidx_)->FV()->Facets(); i++ )
      {
      // identifying the finite volumes to which the flux will be distributed
      this->gref_.E(stencil_.eidx_)->FV()->FacetEdgeNodes( i, stencil_.inside_node_, stencil_.outside_node_ );

      // for the "inside" node
      if ( stencil_.facet_flux_[i] < 0. ) {
         // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
         //    (outside node = upstream)
         // -------------------------------------------------------------------------------------------------
         RESULT[ this->gref_.E(stencil_.eidx_)->N(stencil_.inside_node_)->Idx() ]  +=
               stencil_.facet_flux_[i] * stencil_.psi1_[stencil_.outside_node_];
         // 2. outgoing fluxes are subtracted
         // ---------------------------------
         RESULT[ this->gref_.E(stencil_.eidx_)->N(stencil_.outside_node_)->Idx() ] -=
               stencil_.facet_flux_[i] * stencil_.psi1_[stencil_.outside_node_];
         } 
      else {
         RESULT[ this->gref_.E(stencil_.eidx_)->N(stencil_.inside_node_)->Idx() ]  +=
               stencil_.facet_flux_[i] * stencil_.psi1_[stencil_.inside_node_]; // out
         RESULT[ this->gref_.E(stencil_.eidx_)->N(stencil_.outside_node_)->Idx() ] -=
               stencil_.facet_flux_[i] * stencil_.psi1_[stencil_.inside_node_];
         }
      }

} // end AccumulateFluxUpwindProducts




/** Where the condition flag is not DIRICH, results are mapped from solution
vector<double> back to finite volume cells after performing a range cheque.

@section messages Messages

If the result is out of range at more than 10 nodes, an error message
is returned. If errors occur at more than 2 per cent of the nodes, an
'out_of_range' exception is thrown.

*/
template<uint32_t dim,template<uint32_t> class STP>
double ExplicitNodeCenteredFiniteVolumeTransport<dim,STP>::OutputResults( const PropertyDatabase<dim>& p,
                                                                          const csmp::Index& adv_key,
                                                                          bool show_range,
                                                                          uint32_t var_comp_nr ) const
{
    double  rmin, rmax,
              amin = RESULT[0],
              amax = RESULT[0],
              difference_to_last_output(0.);
    size_t    error_counter(0);

    p.RangeOf( p.Name(adv_key), rmin, rmax );
    if (this->adv1_key_.type==SCALAR){
        for ( size_t i{0U}; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );
            if ( this->gref_.N(i)->Status( adv_key ) != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                double sc = this->gref_.N(i)->Read( adv_key );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc) );
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin )
                    this->gref_.N(i)->Store( adv_key, makeScalar( this->gref_.N(i)->Status(adv_key), RESULT[i] ) );
                else {
                    std::cerr <<"ExplicitNodeCenteredFiniteVolumeTransport<dim>::OutputResults: value: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << std::endl;
                    if ( RESULT[i] > rmax )
                        this->gref_.N(i)->Store( adv_key, makeScalar( this->gref_.N(i)->Status(adv_key), rmax) );
                    else if ( RESULT[i] < rmin )
                        this->gref_.N(i)->Store( adv_key, makeScalar( this->gref_.N(i)->Status(adv_key), rmin) );
                    error_counter++;
                }
            }
        }
    }
    else if (this->adv1_key_.type==ARRAY){
        ArrayVariable av;
        for ( size_t i{0U}; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );

            // array variables have only one flag, hence they are read similar to Scalar variables
            if ( this->gref_.N(i)->Status( adv_key ) != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                this->gref_.N(i)->Read( adv_key, av );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-av(var_comp_nr) ));
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin ){
                    av(var_comp_nr)=RESULT[i];
                    this->gref_.N(i)->Store( adv_key, av );
                }
                else {
                    std::cout <<"ExplicitNodeCenteredFiniteVolumeTransport<dim>::OutputResults: array comp.number: "<<var_comp_nr<< " Result[i]:"<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << std::endl;
                    if ( RESULT[i] > rmax ){
                        av(var_comp_nr)=rmax;
                        this->gref_.N(i)->Store( adv_key, av );
                    }
                    else if ( RESULT[i] < rmin ){
                        av(var_comp_nr)=rmin;
                        this->gref_.N(i)->Store( adv_key, av );
                    }
                    error_counter++;
                }
            }
        }
    }
    else if (this->adv1_key_.type==FLAGGEDARRAY){
        FlaggedArrayVariable fav;
        for ( size_t i{0U}; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );
            VARIABLE_FLAG flag=this->gref_.N(i)->Status( adv_key ,var_comp_nr);
            // Flagged array variables as well as vectors and tensors need a special status call with no return value (it seems)
            if ( flag != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                this->gref_.N(i)->Read( adv_key, fav );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-fav(var_comp_nr) ));
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin ){
                    fav.Component(var_comp_nr,RESULT[i]);
                    this->gref_.N(i)->Store( adv_key, fav );
                }
                else {
                    std::cout <<"ExplicitNodeCenteredFiniteVolumeTransport<dim>::OutputResults: Flagged array comp.number: "<<var_comp_nr<< " Result[i]:"<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << std::endl;
                    if ( RESULT[i] > rmax ){
                        fav.Component(var_comp_nr,rmax);
                        this->gref_.N(i)->Store( adv_key, fav );
                    }
                    else if ( RESULT[i] < rmin ){
                        fav.Component(var_comp_nr,rmin);
                        this->gref_.N(i)->Store( adv_key, fav );
                    }
                    error_counter++;
                }
            }
        }
    }
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if (error_counter>0){
        std::cout <<"\n "<<error_counter<<" nodes have a transport variable range outside that of the PropertyDatabase: "<< rmin <<"-"<< rmax << std::endl;
        std::cout<< " Range is from: "<<amin<<" to: "<<amax<< std::endl;
    }
    // reporting problems
    if ( error_counter >= csmp_error.MaximumNumberOfErrors() ) {
        throw csmp::Exception( ERROR, "ExplicitNodeCenteredFiniteVolumeTransport::OutputResults",
                               "Output property was out of range, legal (min/max) was stored instead");
    }
    if ( error_counter > (RESULT.size() / csmp_error.MaximumNumberOfErrors()) ){
        std::cerr<<"errors: "<<error_counter<<" max: " << csmp_error.MaximumNumberOfErrors()<<" result size:"<<RESULT.size()<<std::endl;
        throw std::out_of_range("ExplicitNodeCenteredFiniteVolumeTransport::OutputResults: Advected variable out of range");
    }

    if ( show_range ) {
         std::cout <<"\n\nExplicitNodeCenteredFiniteVolumeTransport<"<< dim;
         std::cout <<">::OutputResults: Variable range after advection: ";
         this->gref_.MinMaxOf(this->adv1_key_, amin, amax); // takes into account, that DIRICH nodes were not changed
         std::cout << amin <<" to "<< amax << std::endl << std::endl;
      }

    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults






} // end namespace csmp



#endif
