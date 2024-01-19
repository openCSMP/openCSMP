#include "ExplicitMassBasedTransport.h"
#include "MassBasedStencilProcessor.h"

using namespace std;

namespace csmp{

/// single-phase passive advection (group-restricted)
/// with different lhs and rhs
template<uint32_t dim,template<uint32_t> class STP>
ExplicitMassBasedTransport<dim,STP>::ExplicitMassBasedTransport(const char* group,
                                                                Model<dim>& sg,
                                                                const char* porosity,
                                                                const char* advected_prop_lhs,
                                                                const char* advected_prop_rhs,
                                                                const char* transp_velocity,
                                                                const char* nodal_source,
                                                                bool second_order_in_space,
                                                                bool second_order_in_time,
                                                                const char *thickness)
    : ExplicitNodeCenteredFiniteVolumeTransport<dim, STP>( group, sg, porosity,
                                                           advected_prop_lhs, advected_prop_rhs, transp_velocity,
                                                           nodal_source, second_order_in_space, second_order_in_time,thickness),
      ad_rhs_key_(sg.Database ().StorageKey (advected_prop_rhs))
{
    if(this->Verbose())
        cout <<"ExplicitMassBasedTransport<"<< dim <<">(constructor): Constructed successfully."<< endl;
    //! important
    /// if lower dimensional elements are included in transport calculation
    /// 2. multiply FacetArea by thickness in InitializeFVData of NodeCenteredFiniteVolumeTransport
    ///
} // end constructor (solute transport-only)

template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AssignFluxBoundaryConditions( uint32_t var_comp_nr )
{
    /// unneccessary for compressible flow(??) (there needs to be a better comment here! julian, july 2014)
    double        inflow, flux_balance;
    const double  zero(0.);

    if (this->Verbose())
        cout<<" Correcting fluxes in ExplicitMassBasedTransport"<<endl;
    // loop over the boundary FV cells/Nodes and adjust fluxes

    if (this->adv1_key_.type == SCALAR)

        for ( typename std::vector<Node<dim>*>::const_iterator
              nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
        {
            // if the flux balance cannot be evaluated because the node sits at the model boundary
//            if ((*nit)->Status(this->adv1_key_)!=DIRICH) {
                if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) )
                {
                    // inflow or outflow compensation
                    if ( inflow != zero ) this->RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * -inflow;
                }
                // if we have a boundary flux balance because we are dealing with interior boundaries
                else if ( flux_balance != zero ) // counter balancing
                    this->RESULT[ (*nit)->Idx() ] += (*nit)->Read( this->adv1_key_ ) * -flux_balance;
//            }
        }
    else if (this->adv1_key_.type == ARRAY)
        for ( auto nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
        {
//            if ((*nit)->Status(this->adv1_key_)!=DIRICH) {
                ArrayVariable av;
                (*nit)->Read( this->adv1_key_ ,av);
                // if the flux balance cannot be evaluated because the node sits at the model boundary
                if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) )
                {
                    // inflow or outflow compensation
                    if ( inflow != zero ) this->RESULT[ (*nit)->Idx() ] += av(var_comp_nr) * -inflow;
                }
                // if we have a boundary flux balance because we are dealing with interior boundaries
                else if ( flux_balance != zero ) // counter balancing
                    this->RESULT[ (*nit)->Idx() ] += av(var_comp_nr) * -flux_balance;
//            }
        }
    else if (this->adv1_key_.type == FLAGGEDARRAY)

        for ( typename std::vector<Node<dim>*>::const_iterator
              nit=this->gref_.PerimeterNodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
        {
//            if ((*nit)->Status(this->adv1_key_)!=DIRICH) {
                FlaggedArrayVariable fav;
                (*nit)->Read( this->adv1_key_ ,fav);
                // if the flux balance cannot be evaluated because the node sits at the model boundary
                if ( !this->FluxThroughBoundaryFiniteVolume( (*nit), inflow, flux_balance ) )
                {
                    // inflow or outflow compensation
                    if ( inflow != zero ) this->RESULT[ (*nit)->Idx() ] += fav(var_comp_nr) * -inflow;
                }
                // if we have a boundary flux balance because we are dealing with interior boundaries
                else if ( flux_balance != zero ) // counter balancing
                    this->RESULT[ (*nit)->Idx() ] += fav(var_comp_nr) * -flux_balance;
//            }
        }
} // end AssignFluxBoundaryConditions

template<uint32_t dim,template<uint32_t> class STP>
double  ExplicitMassBasedTransport<dim, STP>::AnisotropicCourantIncrement()
{
    this->UpdateProjectedVelocitiesAndFluxBalances();

    VectorVariable<dim>         vc;
    const double              zero(0.);
    double                    courant_increment(8640000.); // 100 days
    double                    velocity;
    
    // loop over the elements finding their transsect length in the direction of flow
    for ( typename vector<Element<dim>*>::const_iterator
          eit=this->gref_.CellsBegin(); eit!=this->gref_.CellsEnd(); eit++ )
    {
        (*eit)->Read( this->vel_key_, vc );
        velocity = vc.Length();

        vc       /= velocity; // normalize vc to avoid round-off error during geometrical projection
        double ediameter((*eit)->LengthInDirection( vc ) * (*eit)->Read( this->phi_key_ ) );

        // guarding against degenerate cases
        if ( velocity > zero and ediameter > zero ) {
            courant_increment = std::min( courant_increment, ediameter / velocity );
        }
        
    }

    if (this->Verbose()) cout <<"\nExplicitMassBasedTransport<dim>::AnisotropicCourantIncrement: "<< courant_increment;
    if ( courant_increment <= 1.0e-3 ) cout <<" \n (CFL constraint is very tight)."<<endl;

    return courant_increment;

} // end AnisotropicCourantIncrement
/**
    Accumulates facet fluxes (volume * saturation) coming into the control volumes into the result vector.
*/
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AccumulateFluxUpwindProducts()
{
    // for all finite-volume facets
    for ( auto i{0U}; i<this->gref_.E(this->stencil_.eidx_)->FV()->Facets(); i++ )
    {
        // identifying the finite volumes to which the flux will be distributed
        this->gref_.E(this->stencil_.eidx_)->FV()->FacetEdgeNodes( i, this->stencil_.inside_node_, this->stencil_.outside_node_ );

        // for the "inside" node
        if ( this->stencil_.facet_flux_[i] < 0. ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            this->RESULT[ this->gref_.E(this->stencil_.eidx_)->N(this->stencil_.inside_node_)->Idx() ]  +=
                    this->stencil_.facet_flux_[i] * this->stencil_.psi2_[this->stencil_.outside_node_];
            // 2. outgoing fluxes are subtracted
            // ---------------------------------
            this->RESULT[ this->gref_.E(this->stencil_.eidx_)->N(this->stencil_.outside_node_)->Idx() ] -=
                    this->stencil_.facet_flux_[i] * this->stencil_.psi2_[this->stencil_.outside_node_];
        }
        else {
            this->RESULT[ this->gref_.E(this->stencil_.eidx_)->N(this->stencil_.inside_node_)->Idx() ]  +=
                    this->stencil_.facet_flux_[i] * this->stencil_.psi2_[this->stencil_.inside_node_]; // out
            this->RESULT[ this->gref_.E(this->stencil_.eidx_)->N(this->stencil_.outside_node_)->Idx() ] -=
                    this->stencil_.facet_flux_[i] * this->stencil_.psi2_[this->stencil_.inside_node_];
        }

    }
}

template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AccumulateFluxUpwindProductsOMP(vector<double>& RESULT)
{
#if defined(_OPENMP )
    size_t tid = omp_get_thread_num();
    Element<dim>* ep = this->gref_.E(this->thread_stencil_processor_[tid]->eidx_);

    // for all finite-volume facets
    for ( auto i{0U}; i<ep->FV()->Facets(); i++ )
    {
        // identifying the finite volumes to which the flux will be distributed
        ep->FV()->FacetEdgeNodes( i, this->thread_stencil_processor_[tid]->inside_node_, this->thread_stencil_processor_[tid]->outside_node_ );

        // for the "inside" node
        if ( this->thread_stencil_processor_[tid]->facet_flux_[i] < 0. ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            RESULT[ ep->N(this->thread_stencil_processor_[tid]->inside_node_)->Idx() ]  +=
                    this->thread_stencil_processor_[tid]->facet_flux_[i] * this->thread_stencil_processor_[tid]->psi2_[this->thread_stencil_processor_[tid]->outside_node_];
            // 2. outgoing fluxes are subtracted
            // ---------------------------------
            RESULT[ ep->N(this->thread_stencil_processor_[tid]->outside_node_)->Idx() ] -=
                    this->thread_stencil_processor_[tid]->facet_flux_[i] * this->thread_stencil_processor_[tid]->psi2_[this->thread_stencil_processor_[tid]->outside_node_];
        }
        else {
            RESULT[ ep->N(this->thread_stencil_processor_[tid]->inside_node_)->Idx() ]  +=
                    this->thread_stencil_processor_[tid]->facet_flux_[i] * this->thread_stencil_processor_[tid]->psi2_[this->thread_stencil_processor_[tid]->inside_node_]; // out
            RESULT[ ep->N(this->thread_stencil_processor_[tid]->outside_node_)->Idx() ] -=
                    this->thread_stencil_processor_[tid]->facet_flux_[i] * this->thread_stencil_processor_[tid]->psi2_[this->thread_stencil_processor_[tid]->inside_node_];
        }

    }
#endif
} // end AccumulateFluxUpwindProductsOMP

template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::ComposeSolution(double time_interval, uint32_t var_comp_nr )
{
    if (this->adv1_key_.type == SCALAR)
        for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
        {
            // reading the transported variable value
            double solution = this->gref_.N(nidx)->Read( this->adv1_key_ );
            // subtracting the flux time-interval product
            this->RESULT[nidx] = solution - (time_interval / this->FVPOREVOL[nidx]) * this->RESULT[nidx];
            // adding externally assigned source or sink terms
            this->RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );
        }
    else if (this->adv1_key_.type == ARRAY)
    {
        for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
        {
            // reading the transported variable value
            ArrayVariable solution;
            this->gref_.N(nidx)->Read( this->adv1_key_ ,solution);
            // subtracting the flux time-interval product
            this->RESULT[nidx] = solution(var_comp_nr) - (time_interval / this->FVPOREVOL[nidx]) * this->RESULT[nidx];
            // adding externally assigned source or sink terms
            this->RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );
        }
    }
    else if (this->adv1_key_.type == FLAGGEDARRAY){
        for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ )
        {
            // reading the transported variable value
            FlaggedArrayVariable solution;
            this->gref_.N(nidx)->Read( this->adv1_key_ ,solution);
            // subtracting the flux time-interval product
            this->RESULT[nidx] = solution(var_comp_nr) - (time_interval / this->FVPOREVOL[nidx]) * this->RESULT[nidx];
            // adding externally assigned source or sink terms
            this->RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );
        }
    }
} // end ComposeSolution (passive advection case)

template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AdvectVariable( double time_interval)
{
    //TODO if porevolume changes -> updateporevolumes
    this->gref_.RenumberNodes();

    //if (true) // update pore volumes
    //  {
    //  this->InitializeSectorPoreVolumeData (true);
    //  this->InitializeArraysForFirstOrderMethod();
    //  }

    this->UpdateProjectedVelocitiesAndFluxBalances();
    bool output_result_range(false);
    // explicit 1st-order solution of advection equation
    this->AdvectVariable1stOrder( time_interval, output_result_range );

} // end AdvectVariable

template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AdvectVariable1stOrder(
        double time_increment, bool output_result_range )
{
    for ( uint32_t ncom{0u}; ncom < this->var_ncomponents_;ncom++){
        if (this->Verbose()) cout<<" Advecting component (mass based): "<<ncom<<" time increment used: "<<time_increment<<endl;
        std::fill( this->RESULT.begin(), this->RESULT.end(), static_cast<double>(0.) );

#if !defined(_OPENMP)
        /// IMPORTANT NOTE: If you every change/improve this loop, make sure you
        /// input the equivalent changes in the openmp section below. Otherwise
        /// your changes/improvements will only be evident in serial simulations.
        std::vector<FV_Parameter>::const_iterator  fvt = this->STENCIL_DATA.begin();
        for ( typename std::vector<Element<dim>*>::const_iterator
              eit=this->gref_.CellsBegin();
              eit!=this->gref_.CellsEnd(); eit++, fvt++ )
        {
            this->stencil_.eidx_ = (*eit)->Idx();

            // getting all necessary data for the construction of the result vector from the element
            this->stencil_.InitializeFirstOrder( (*fvt), *(*eit) ,this->adv1_key_.type,ncom);
            // starting accumulation with in and out fluxes which are temporarily stored in result vector
            AccumulateFluxUpwindProducts();

        } // end of accumulation
#else
        vector<vector<double> > RESULT(omp_get_max_threads());
        for (size_t tid = 0 ; tid < omp_get_max_threads();tid++){
            RESULT[tid].resize(this->gref_.Nodes());
            std::fill( RESULT[tid].begin(), RESULT[tid].end(), static_cast<double>(0.) );
        }

#pragma omp parallel
        {

            Element<dim>* ep;
            size_t tid = omp_get_thread_num();
            /// OPENMP note: given the likely small size of ncom vs the number of elements,
            /// only the inner loop was parallelized as it was assumed to be the most efficient way.
            /// Perhaps an outer dynamic scheduled for loop instead of the inner one should be tried
            /// Julian - 02/09/2015

#pragma omp for
            for ( int32_t e = 0 ; e < this->gref_.Cells(); e++ )
            {
                ep = this->gref_.E(e);
                this->thread_stencil_processor_[tid]->eidx_ = ep->Idx();

                //----------------------------------------------------
                //change the element stencil to one for this thread, temporarily.
                FiniteElement* fe_tmp=ep->FE();
                // change pointer here
                ep->Assign(this->femgrs_[tid].E(ep->FE_Type()));
                //----------------------------------------------------

                //reassign fv stencil here (to put it back later, at the end of the loop.
                const FiniteVolumeStencil<dim>* tmp_fvstencil= ep->FV();
                ep->Assign(this->fvmgrs_[tid].Stencil( ep->FE_Type()));

                // getting all necessary data for the construction of the result vector from the element
                this->thread_stencil_processor_[tid]->InitializeFirstOrder( this->STENCIL_DATA[e], *ep ,this->adv1_key_.type,ncom);
                // starting accumulation with in and out fluxes which are temporarily stored in result vector

                AccumulateFluxUpwindProductsOMP(RESULT[tid]);
                //put the stencil back
                ep->Assign(tmp_fvstencil);

                // put the FEM back here
                ep->Assign(fe_tmp);

            } // end of accumulation
        }
        for  (size_t tid = 0 ; tid < omp_get_max_threads();tid++){
            for (auto i = 0 ; i < this->gref_.Nodes();i++)
                this->RESULT[i]+=RESULT[tid][i];
        }
#endif

        // compensate for the inflow and the outflow boundaries
        this->AssignFluxBoundaryConditions(ncom);

        // 4. Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
        //     - a single loop over all nodes is required
        //     - this steps also considers in and outflow of finite volume cells
        //     - (distributed) sources and sinks will be considered in the future
        this->ComposeSolution( time_increment,ncom);


        // 5. Saving the computed new saturations at the FV centers
        this->OutputResults( this->pref_, this->adv1_key_, output_result_range,ncom);
    }


} // end AdvectVariable1stOrder

/**
  @author: Julian E. Mindel (23-09-2013)

This method calculates a single timestep of transport, assumming nothing about the time interval
it is "fed".  It will advect using that time interval, assuming the user has externally determined
that this is the correct time interval length.  In contrast with AdvectVariable(, this method does
not break the time interval into sub-parts and guarantee stability. This also means that the courant
increment is not checked for in this method (and should be, externally, of course).

Note, that in this case, AdvectVariable as created before the date in this documentation entry did not
contain subcycling of time intervals.  This statement is valid only for ExplicitMassBasedTransport.

 */
template<uint32_t dim,template<uint32_t> class STP>
void ExplicitMassBasedTransport<dim,STP>::AdvectVariableSingleStep( double time_increment,
                                                                    bool apply_flux_balance_correction,
                                                                    bool update_pore_volumes )
{
    this->gref_.RenumberNodes();
    this->UpdateProjectedVelocitiesAndFluxBalances();
    bool output_result_range(false);
#ifdef DEBUG_ExplicitNodeCenteredFiniteVolumeTransport
    output_result_range=true;
#endif

    // 1. update the sector volumes and FVPOREVOL if so required (as in heat transport for instance)
    // ---------------------------------------------------------------------------------------------
    if ( update_pore_volumes ) {
        const bool multiply_with_thickness_attribute = (this->thi_key_ != csmp::Index() ) ? true : false;
        this->InitializeSectorPoreVolumeData( multiply_with_thickness_attribute );
        this->InitializeArraysForFirstOrderMethod();
    }

    if ( !apply_flux_balance_correction && this->Verbose()) {
        std::cout <<"ExplicitMassBasedTransport<"<< dim <<">::AdvectVariable: ";
        std::cout <<" Flux balance correction is applied automatically."<< std::endl;
    }

    // 2. compute solution
    // -------------------
    if (this->Verbose()) std::cout <<"\n\n\ExplicitMassBasedTransport<"<< dim;
    if (this->Verbose()) std::cout <<">::AdvectVariableSingleStep: Advecting transport variable";
        if ( !this->SecondOrderInSpace() )
        if ( this->diff_key_ != csmp::Index() )
            // explicit 1st-order solution of advection-diffusion equation
            throw csmp::Exception( ERROR, "ExplicitMassBasedTransport<dim>::AdvectVariableSingleStep",
                                         " No diffusion is possible in ExplicitMassBasedTransport" );
        else
        {
            // explicit 1st-order solution of advection equation
            this->AdvectVariable1stOrder( time_increment, output_result_range );
        }

    else {
        // explicit 2nd-order solution
        throw csmp::Exception( ERROR, "ExplicitMassBasedTransport<dim>::AdvectVariableSingleStep",
                                     " No 2nd order is possible in ExplicitMassBasedTransport" );
    }

} // end AdvectVariableSingleStep


template class ExplicitMassBasedTransport<1U,MassBasedStencilProcessor>;
template class ExplicitMassBasedTransport<2U,MassBasedStencilProcessor>;
template class ExplicitMassBasedTransport<3U,MassBasedStencilProcessor>;

} // end namespace csmp
