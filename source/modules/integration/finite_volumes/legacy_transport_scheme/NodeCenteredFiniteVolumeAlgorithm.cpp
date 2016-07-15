#include "NodeCenteredFiniteVolumeAlgorithm.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "StencilProcessor.h"
#include "ExplicitStencilProcessor.h"
#include "DenseMatrix.h"
#include "Element.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Solver.h"
#include "GaussJordan_Solver.h"
#include "TwoPhaseModel.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

/** Initializes LHS, XVEC, and RHS arrays.
*/
template<size_t dim>
NodeCenteredFiniteVolumeAlgorithm<dim>::NodeCenteredFiniteVolumeAlgorithm( Region<dim>& sg )
 : gref_(sg),
   RESULT( sg.Nodes(), 0. ),
   DN(dim,3), DNT(3,dim),
   MAX_NODES_GAUSS_SOLVER(50U),
   solver_(0),
   firstCall_(true),
   verbose_(true)
 {
    RHS.resize( gref_.Nodes(), 0. );
    LHS.Resize( gref_.Nodes() );

    if ( RHS.size() <= MAX_NODES_GAUSS_SOLVER ) solver_ = new GaussJordan_Solver();
    else {

   solver_ = new CSMP_DEFAULT_LINEAR_SOLVER();
#ifdef CSMP_WITH_SAMG_SOLVER
   GetSolverSettings().Set_iout1(1);
   GetSolverSettings().Set_iout2(0);
#endif

   }
      
 } // end ctor


 
 
template<size_t dim>
NodeCenteredFiniteVolumeAlgorithm<dim>::~NodeCenteredFiniteVolumeAlgorithm()
 {
    delete solver_;
 }





/** Removes potential leftover entries from previous solution matrix or
re-initializes the righthand vector<double64> to zero.  

@section arguments Input Arguments 

The size of the computational problem = number of FVs.
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::ResetLHS( size_t nodes )
 {
    LHS.Erase();
    LHS.Resize( nodes );
 }



template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::ResetRHS( size_t nodes )
 {
    // if the RHS was used before
    if ( !RHS.empty() ) {
          RHS.resize( nodes );
          vector<double64>(RHS).swap(RHS);
          fill( RHS.begin(), RHS.end(), static_cast<double64>(0.) );
       }
    // if for some reason there is no RHS vector
    else
    fill( RESULT.begin(), RESULT.end(), static_cast<double64>(0.) );
 }

template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::ResetRow( size_t nid )
{
    LHS.ZeroRow(nid);
    RHS[nid]=0.0;
}

template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::ZeroRHS( size_t nid )
{
    RHS[nid]=0.0;
}

/// returns the settings object if SAMG is used to solve the matrix equations
template<size_t dim>
CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS& NodeCenteredFiniteVolumeAlgorithm<dim>::GetSolverSettings()
{
    return *dynamic_cast<CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS*>((solver_)->GetSolverSettings());
}

/** returns the solver object if SAMG is used to solve the matrix equations
      It is important to note that this, just like the old method above (GetSolverSettings), asks for an instance of
      SAMG_Solver that may not exist, given that if the problem is smaller than 50 nodes, a Gauss_Jordan solver is created by
      default.  I believe this needs re-consideration.  Rather , I would utilize SAMG_Solver even for small problems, with special
      settings.  -- Julian 24-09.2013
*/
template<size_t dim>
CSMP_DEFAULT_LINEAR_SOLVER* NodeCenteredFiniteVolumeAlgorithm<dim>::GetSolver()
    {
        return dynamic_cast<CSMP_DEFAULT_LINEAR_SOLVER*>(solver_);
    }

/**
 
The transport velocity values projected onto the facet normals and multiplied
with the facet areas are multiplied with the upstream values of the 
transported variable at each facet integration point and accumulated into 
the solution vector<double64> as the upstream weighted first-order fluxes.  

@section arguments Input Arguments 

To assign the higher order fluxes, the method needs access to the 
current finite element and the associated FV stencil processor.  
 */
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateFluxUpwindSaturationProducts( const StencilProcessor<dim>& es )
{
   // for all finite-volume facets
  for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
    {
       // identifying the finite volumes to which the flux will be distributed
       gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

       // for the "inside" node
       if ( es.facet_flux_[i] < 0. ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            RESULT[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ]  += 
                                                                    es.facet_flux_[i] * es.psi1_[es.outside_node_];
            // 2. outgoing fluxes are subtracted
            // ---------------------------------
            RESULT[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ] -= 
                                                                    es.facet_flux_[i] * es.psi1_[es.outside_node_];
         } 
       else {
            RESULT[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ]  += 
                                                              es.facet_flux_[i] * es.psi1_[es.inside_node_]; // out
            RESULT[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ] -= 
                                                              es.facet_flux_[i] * es.psi1_[es.inside_node_];
         }
    }

 } // end AccumulateFluxUpwindSaturationProducts






/**
 
The transport velocity values projected onto the facet normals and multiplied
with the facet areas are multiplied with the (limited) 'isat' variable values 
at each facet integration point to define the higher order fluxes.  

@section arguments Input Arguments 

To assign the higher order fluxes, the method needs access to the 
current finite element and the associated FV stencil processor.  
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateHigherOrderFluxSaturationProducts( const StencilProcessor<dim>& es )
{
   // for all finite-volume facets
   for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
     {
         // identifying the finite volumes to which the flux will be distributed
         gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

         // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
         // -------------------------------------------------------------------------------------------------
         RESULT[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ]  += es.facet_flux_[i] * es.ipsi1_[i];
         
         // 2. outgoing fluxes are subtracted
         // ---------------------------------
         RESULT[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ] -= es.facet_flux_[i] * es.ipsi1_[i];
     }

 } // end AccumulateHigherOrderFluxSaturationProducts 	                                       

/** Uses the finite element method to add diffusion to the FV scheme.
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateIntegral_DNT_op_DN_dV_LHS( 
                                                                           const StencilProcessor<dim>& es )
 {

    for ( size_t i=0U; i<gref_.E(es.eidx_)->FE()->IntegrationPoints(); i++ )
      {
         // getting global intpol. function derivative matrix and determinant of
         // byproduct Jacobian matrix (B is already in global coordinates)
         double64 detJ = (*gref_.E(es.eidx_)).dN_AtIntegrationPoint( DN, i, SCALAR );
         DN.Transposed( DNT );
         DN  *= es.diff_coeff_;
         DNT *= DN;
         DNT *= (*gref_.E(es.eidx_)).WeightAtIntegrationPoint(i) * detJ; 

         // assigning the matrix contribution to the solution matrix
         for ( size_t j=0U; j<DNT.Rows(); j++ )
           for ( size_t k=0U; k<DNT.Cols(); k++ )
             LHS.Add( gref_.E(es.eidx_)->N(j)->Idx(),
                      gref_.E(es.eidx_)->N(k)->Idx(), DNT(j,k) );
      }


 } // end AccumulateIntegral_DNT_op_DN_dV_LHS





template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateIntegral_DN_op_dS_LHS( const StencilProcessor<dim>& es )
 {

    gref_.E(es.eidx_)->dN_AtBaryCenter(DN);
    DN  *= -es.diff_coeff_;

    double64 flux(0.0);
    Point<dim>  n;
    for ( size_t iFacet=0U; iFacet<gref_.E(es.eidx_)->FV_Stencil()->Facets(); iFacet++ )
    {
        gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( iFacet, es.inside_node_, es.outside_node_ );

        n = gref_.E(es.eidx_)->FacetNormal(iFacet);

        for ( size_t j=0U; j<DN.Cols(); j++ ){

            flux = 0.0;
            for ( size_t k=0U; k< DN.Rows(); k++ )
                flux += DN(k,j)*n[k];
            flux *= gref_.E(es.eidx_)->FacetArea(iFacet);

            LHS.Add( gref_.E(es.eidx_)->N(es.inside_node_)->Idx(),
                     gref_.E(es.eidx_)->N(j)->Idx(),  flux );

            LHS.Add( gref_.E(es.eidx_)->N(es.outside_node_)->Idx(),
                     gref_.E(es.eidx_)->N(j)->Idx(), -flux );

        }

    }


 } // end AccumulateIntegral_DN_op_dS_LHS






/** Source term per finite volume sector is added to the diagonal of the
matrix.  
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateSectorSourceTermsInLHS( const StencilProcessor<dim>& es ) 
 {
    for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Sectors(); i++ )
      {
         // identifying the finite volumes to which the flux will be distributed
         const size_t j(gref_.E(es.eidx_)->N(i)->Idx());
         
         // diffusion / nodal sources or sink terms are add to the diagonal of lhs
         LHS.Add( j, j, es.src_[i] );
      }
      
 } // end AccumulateSectorSourceTermsInLHS











/**

Accumulates the volume over the time-increment and the outgoing facet-
integrated fluxes into the diagonal of the solution matrix.  

An upstream coupling is created by accumulating the fluxes that come
into the finite-volumes sectors as off-diagonal elements. For each
finite-volume facet, there are the nodes i (inside) and j(outside).  

@section arguments Input Arguments 

The current Element and uptodate FiniteVolumeProcessor are used to 
retrieve the variables of interest.  
 */
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateLHS( const StencilProcessor<dim>& es,
                                                                double64 time_multiplier )
{
  const double64 zero(0.);
  
  // putting contributions to pore volume into the matrix diagonal
  for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Sectors(); i++ )
    //                                 phi * sector-volume
    LHS.Add( gref_.E(es.eidx_)->N(i)->Idx(),    
             gref_.E(es.eidx_)->N(i)->Idx(), 
             es.sector_pore_volume_[i] / time_multiplier );
    
   // for all finite-volume facets
  for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
    {
       // identifying the finite volumes to which the flux will be distributed
       gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

       if ( es.facet_flux_[i] < zero ) {
            // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
            //    (outside node = upstream)
            // -------------------------------------------------------------------------------------------------
            LHS.Add( gref_.E(es.eidx_)->N(es.inside_node_)->Idx(), 
                     gref_.E(es.eidx_)->N(es.outside_node_)->Idx(),  
                     es.facet_flux_[i] );  // incoming flux
            
            // 2. outgoing fluxes are added to the matrix diagonal
            // ---------------------------------------------------
            LHS.Add( gref_.E(es.eidx_)->N(es.outside_node_)->Idx(), 
                     gref_.E(es.eidx_)->N(es.outside_node_)->Idx(),
                    -es.facet_flux_[i] ); // outgoing flux
         } 
       else {
            LHS.Add( gref_.E(es.eidx_)->N(es.inside_node_)->Idx(), 
                     gref_.E(es.eidx_)->N(es.inside_node_)->Idx(),   
                     es.facet_flux_[i] );  // outgoing flux
            LHS.Add( gref_.E(es.eidx_)->N(es.outside_node_)->Idx(), 
                     gref_.E(es.eidx_)->N(es.inside_node_)->Idx(), 
                    -es.facet_flux_[i] );  // incoming flux
         }
    }
} // end AccumulateLHS








template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateMatrix_NonlinearNewtonRaphson( const StencilProcessor<dim>& es,
                                                                                       std::vector<double64>& SAT0,
                                                                                       double64 time_multiplier)
 {
    // putting contributions to pore volume into the matrix diagonal
    for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Sectors(); i++ ){
      //  phi * sector-volume
      LHS.Add( gref_.E(es.eidx_)->N(i)->Idx(),
               gref_.E(es.eidx_)->N(i)->Idx(),
               (es.sector_pore_volume_[i] / time_multiplier) );
    }

    // for all finite-volume facets
    for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
      {
         // identifying the finite volumes to which the flux will be distributed
         gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

         LHS.Add( gref_.E(es.eidx_)->N(es.inside_node_)->Idx(),
                  gref_.E(es.eidx_)->N(es.upstream_node_[i])->Idx(),
                  es.facet_flux_[i] );  // if (es.facet_flux_[i]>0)  - incoming flux, else outgoing flux

         LHS.Add( gref_.E(es.eidx_)->N(es.outside_node_)->Idx(),
                  gref_.E(es.eidx_)->N(es.upstream_node_[i])->Idx(),
                  -es.facet_flux_[i] ); // if (es.facet_flux_[i]>0)  - outgoing flux, else incoming flux
      }

 } // end AccumulateMatrix_NolinearNewtonRaphson












template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateResidual_NonlinearNewtonRaphson( const StencilProcessor<dim>& es,
                                                                                        std::vector<double64>& SAT0,
                                                                                        double64 time_multiplier)
 {
    Element<dim>* e(gref_.E(es.eidx_));
    size_t NNodes(e->Nodes());
    // fill righthandside with the  prop_t0 * pore_vol/time_increment  products
    for ( size_t i=0U; i<NNodes; i++ )
      RHS[ e->N(i)->Idx() ] +=
            - ((es.psi1_[i] - SAT0[ e->N(i)->Idx() ]) * es.sector_pore_volume_[i]) / time_multiplier;

    for ( size_t i=0U; i<NNodes; i++ )
        RHS[ e->N(i)->Idx() ] +=
            -es.facet_flux_rhs_[i]
            ;

 } // end AccumulateResidual_NolinearNewtonRaphson








template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateMatrixAtBoundary_NonlinearNewtonRaphson( const StencilProcessor<dim>& es, std::vector<double64>& SAT0, size_t nid, size_t pnid, double64 time_multiplier,const csmp::Index& adv1_key)
 {

    if(gref_.N(pnid)->Status( adv1_key ) != DIRICH ){
        // putting contributions to pore volume into the matrix diagonal
        //  phi * sector-volume
        LHS.Add( gref_.E(es.eidx_)->N(pnid)->Idx(),
                 gref_.E(es.eidx_)->N(pnid)->Idx(),
                 (es.sector_pore_volume_[pnid] / time_multiplier) );

        for ( size_t k=0U; k<gref_.E(es.eidx_)->FV_Stencil()->FacetsPerSector(pnid); k++ ) {
              const size_t i(gref_.E(es.eidx_)->FV_Stencil()->FacetSurroundingSector(pnid,k));
              // if the sector node is the inside node then an incoming flux will create a positive source term
              gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

              if ( pnid == es.inside_node_ ){
                  LHS.Add( gref_.E(es.eidx_)->N(es.inside_node_)->Idx(),
                           gref_.E(es.eidx_)->N(es.upstream_node_[i])->Idx(),
                           es.facet_flux_[i] );  // if (es.facet_flux_[i]>0)  - incoming flux, else outgoing flux
              }

              if ( pnid == es.outside_node_ ){
                  LHS.Add( gref_.E(es.eidx_)->N(es.outside_node_)->Idx(),
                           gref_.E(es.eidx_)->N(es.upstream_node_[i])->Idx(),
                           -es.facet_flux_[i] ); // if (es.facet_flux_[i]>0)  - outgoing flux, else incoming flux
            }
        }

    }else{

        LHS.Assign( gref_.E(es.eidx_)->N(pnid)->Idx(),
                    gref_.E(es.eidx_)->N(pnid)->Idx(),
                    1.0 );

    }

 } // end AccumulateMatrixAtBoundary_NolinearNewtonRaphson




template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateResidualAtBoundary_NonlinearNewtonRaphson( const StencilProcessor<dim>& es, std::vector<double64>& SAT0, size_t nid, size_t pnid, double64 time_multiplier, const csmp::Index& adv1_key)
 {
    if(gref_.N(pnid)->Status( adv1_key ) != DIRICH ){

      RHS[ gref_.E(es.eidx_)->N(pnid)->Idx() ] +=
            -((es.psi1_[pnid] - SAT0[ gref_.E(es.eidx_)->N(pnid)->Idx() ]) * es.sector_pore_volume_[pnid]) / time_multiplier
            ;
      RHS[ gref_.E(es.eidx_)->N(pnid)->Idx() ] +=
           -es.facet_flux_rhs_[pnid]
           ;
   }else{

        RHS[ gref_.E(es.eidx_)->N(pnid)->Idx() ]=-(es.psi1_[pnid] - SAT0[ gref_.E(es.eidx_)->N(pnid)->Idx() ]);

    }

 } // end AccumulateResidual_NolinearNewtonRaphson


template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::CompensateInflowOutFlowBoundaries(size_t nid,
                                                                       const double64& inflow,
                                                                       const double64& flux_balance)
 {
    if(gref_.IsPerimeterNode(nid))
        // inflow compensation (both in and outflow compensations are necessary in explicit scheme)
        RHS[ nid ]-= (-inflow);
    else
        // if we have a flux balance because this is an internal boundary
        RHS[ nid ] -= (-flux_balance);

 } // end CompensateInflowOutFlowBoundaries




/**
 
Righthand side for Backward-Euler FD time stepping  

    RHS = (FV-volume * sat_t0) / time-increment

@section arguments Input Arguments 

The current Element and up-to-date FiniteVolumeProcessor are used to
retrieve the variables of interest. The saturations from the current
timelevel are supplied by the vector<double64> SAT0. The time-multiplier
is needed to divide the finite-volume sector volume by.  

*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateRHS( const StencilProcessor<dim>& es,
                                                            vector<double64>& SAT0,
                                                            double64 time_multiplier )
 {
     // fill righthandside with the  prop_t0 * pore_vol/time_increment  products
     for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Sectors(); i++ )
       //                                           phi * sector volume
       RHS[ gref_.E(es.eidx_)->N(i)->Idx() ] += 
           (SAT0[ gref_.E(es.eidx_)->N(i)->Idx() ] * es.sector_pore_volume_[i]) / time_multiplier;

 } // end AccumulateRHS



// simpler version for first-order method in time only, O.K.
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateRHS( const StencilProcessor<dim>& es,
                                                            double64 time_multiplier )
 {
     // fill righthandside with the  prop_t0 * pore_vol/time_increment  products
     for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Sectors(); i++ )
       RHS[ gref_.E(es.eidx_)->N(i)->Idx() ] += (es.psi1_[i] * es.sector_pore_volume_[i]) / time_multiplier;

 } // end AccumulateRHS


/**
 
Higher-order fluxes are calculated from the slope-limited saturations at 
the current and previous time levels and the corresponding volume fluxes
integrated over the facets dividing the finite-volume sectors.  

The contributions from the two time levels are weighted using the theta
limiter value computed for each facet. When theta=1 the scheme is fully
implicit, i.e. only the contributions for the new time level are 
considered. Theta limiting attempts however the equally weight the 
contributions on the basis of their magnitude. If this is achieved, the
scheme is second-order accurate in time and numerical diffusion is
minimized.  

The computed higher-order fluxes are subtracted from the righthand vector.
  
@section arguments Input Arguments 

To compute the higher-order fluxes, the method requires access to
the current Element and up-to-date finite volume processor object.  

The method requires the facet-fluxes and the slope-limited interpolated 
values of the transported variable at the facet integration points. The
temporal limiter value, theta, is supplied as a 7th method argument.  

@section implementation Implementation

@section application Application

In higher-order theta-limited transport scheme.  
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateHigherOrderRHS( const StencilProcessor<dim>& es )
 {
     // summing flux saturation products over the finite volume cell
     for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ ) {
	        // computing higher order flux
	        const double64  hflux(es.facet_flux_[i] * es.ipsi1_[i]);
          
          // subtracting higher-order solution from righthand side
          // incoming fluxes, first cell (outward pointing normal)
            RHS[ gref_.E(es.eidx_)->N(gref_.E(es.eidx_)->FV_Stencil()->InsideNode(i))->Idx() ]  -= hflux;
          // incoming fluxes, second cell (inward pointing normal)
            RHS[ gref_.E(es.eidx_)->N(gref_.E(es.eidx_)->FV_Stencil()->OutsideNode(i))->Idx() ] += hflux;
       }

 } // end AccumulateHigherOrderRHS (version for first-order-accurate method in time)

template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AccumulateHigherOrderRHS( 
                                                      const StencilProcessor<dim>& es,
                                                      const vector<vector<double64> >& FACETFLUXES0,
                                                      const vector<vector<double64> >& LTDSATS0 )
 {
     // summing flux saturation products over the finite volume cell
     // theta  n +  1/2 term          higher-order flux contribution
     for ( size_t i=0; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ ) {
  	      // computing higher order flux
  	      double64 hflux  = es.theta_[i]     * es.facet_flux_[i]         * es.ipsi1_[i];
  	      hflux += (1. - es.theta_[i]) * FACETFLUXES0[es.eidx_][i] * LTDSATS0[es.eidx_][i]; 

  	      // identifying the finite volumes to which the flux will be distributed
            // subtracting higher-order solution from righthand side
            // incoming fluxes, first cell (outward pointing normal)
   	      RHS[ gref_.E(es.eidx_)->N(gref_.E(es.eidx_)->FV_Stencil()->InsideNode(i))->Idx() ]  -= hflux;    
            // incoming fluxes, second cell (inward pointing normal)
   	      RHS[ gref_.E(es.eidx_)->N(gref_.E(es.eidx_)->FV_Stencil()->OutsideNode(i))->Idx() ] += hflux;
       }

 } // end AccumulateHigherOrderRHS (version for second-order-accurate method in time)

/**
 
Add first-order Backward Euler solution to RHS (using upstream saturations).
This integrates upstream sats over the segments by summing (flux * saturation) 
products over faces.  

@section arguments Input Arguments 

References to the current element, finite-volume stencil processor,
and the fluxes integrated over the facets inside the current element.  
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AddToRHS( const StencilProcessor<dim>& es )
 {
     const double64 zero(0.);
     double64       uvar; // upstream value of advected variable
 
     for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
       {
	       // identifying the finite volumes to which the flux will be distributed
	       gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );
	       
	       // identifying upstream value of advected variable
	       if ( es.facet_flux_[i] < zero ) uvar = es.psi1_[es.outside_node_];
           else                            uvar = es.psi1_[es.inside_node_];
           
           // assigning upstream fluxes
	  	   RHS[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ]  += es.facet_flux_[i] * uvar;
  		   RHS[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ] -= es.facet_flux_[i] * uvar;
	   }

 } // end AddToRHS



template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AddToRHS( const StencilProcessor<dim>& es, vector<double64>& SAT0  )
 {
     const double64 zero(0.);
     double64       uvar; // upstream value of advected variable

     for ( size_t i=0U; i<gref_.E(es.eidx_)->FV_Stencil()->Facets(); i++ )
       {
           // identifying the finite volumes to which the flux will be distributed
           gref_.E(es.eidx_)->FV_Stencil()->FacetEdgeNodes( i, es.inside_node_, es.outside_node_ );

           // identifying upstream value of advected variable
           if ( es.facet_flux_[i] < zero ) uvar = SAT0[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ];
           else                            uvar = SAT0[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ];

           // assigning upstream fluxes
           RHS[ gref_.E(es.eidx_)->N(es.inside_node_)->Idx() ]  += es.facet_flux_[i] * uvar;
           RHS[ gref_.E(es.eidx_)->N(es.outside_node_)->Idx() ] -= es.facet_flux_[i] * uvar;
       }

 } // end AddToRHS




/** including potential point sources or sinks in the model.
 */
/*
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AddDivergenceOfFluxes( 
                                                   const vector<double64>& SAT0, 
                                                   const vector<double64>& pore_volume,
                                                   const vector<double64>& divsrc,
                                                   double64 dt, 
                                                   bool righthandside_only )
 {
     const double64 zero(0.);
 
     for ( size_t nidx=0U; nidx<divsrc.size(); nidx++ )
         // any nodal source or sink will lead to a divergence of the fluxes
         if ( divsrc[nidx] != zero ) { 
	          // limit the source terms according to equation 61, p. 31 (original manuscript),
	          // but max criterion is not applied
	          double64 omega_i = omega( dt, pore_volume[nidx], divsrc[nidx] );

		      // balance the contribution between left and righthand sides
		      if ( !righthandside_only ) {
		           LHS.Add( nidx, nidx, omega_i * divsrc[nidx] );
		           RHS[ nidx ] -= (static_cast<double64>(1.) - omega_i) * divsrc[nidx] * SAT0[nidx];
		        }
		      else RHS[ nidx ] -= divsrc[nidx] * SAT0[nidx];
		      
		      // verify that lefthandside has something in it
		      if ( static_cast<double64>(0.) >= LHS( nidx, nidx ) ) {
		           cout <<"\NodeCenteredFiniteVolumeAlgorithm::AddDivergenceOfFluxes: LHS diagonal at i=j=";
		           cout << nidx + 1 <<" is close to zero: ";
		           cout << LHS( nidx, nidx ) <<", augmenting value by 1.0e-15."<< endl;
		           LHS.Add( nidx, nidx, 1.0e-15 );
		        }
	       }

 } // end AddDivergenceOfFluxes (point sources or sinks)


// adds contribution to lefthandside no matter what !
// CURRENTLY USED
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::AddDivergenceOfFluxes( const vector<double64>& divsrc )
 {
     for ( size_t nidx=0; nidx<divsrc.size(); nidx++ )
       // any nodal source or sink will lead to a divergence of the fluxes
       if ( divsrc[nidx] != static_cast<double64>(0.) ) LHS.Add( nidx, nidx, divsrc[nidx] );

 } // end AddDivergenceOfFluxes (point sources or sinks)

*/


/**
 
If the fluxbalance is negative, this means that there is outflow from a 
truncated boundary cell that only has inward pointing normals.

This function compensates this negative balance.

@warning THIS FUNCTION ACTUALLY ADDS CONTRIBUTIONS DUE TO NON_CONSERVATIVE OUTFLUXES To RHS
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::SubtractNonConservativeOutFluxesFromRHS( 
                                                                       TwoPhaseModel<dim>&  relperm,  
                                                                       const vector<double64>& FLUX_BALANCE,
                                                                       const csmp::Index& ad_key )
 {
     // pair: local parent node number (0...nodes-1), & global element ID (n...elements)
     double64  flux_balance; 
     double64  psi;

     for ( size_t nidx=0; nidx<FLUX_BALANCE.size(); nidx++ ) 
       // at the model boundaries
       if ( gref_.N(nidx)->AtBoundary() != NOT ) {
            // compute the fractional flow using the first available parent element
            size_t eidx = gref_.N(nidx)->Parent(0u)->Idx();
            
            // setting up the relative permeability model for saturation value at finite volume
            relperm.Initialize( *gref_.E( eidx ) );
            //relperm.InitializeForBaryCenter( e ); // densities and viscosities as when mobility was calculated
            relperm.SaturationWettingPhase( 1. - (psi=gref_.N(nidx)->Read( ad_key )) );
            relperm.EffectiveSaturation();        
            
            // absolute value of mismatch denotes amount that goes out of model
            flux_balance = fabs(FLUX_BALANCE[nidx]) * relperm.AdvectionMultiplier();
            
            // if the fluxbalance is negative, this means that there is inflow, so that 
            // there should be a compensating negative term in the RHS
            if ( flux_balance < 0. ) RHS[ nidx ] -= flux_balance * psi;
         }
    
 } // end SubtractNonConservativeOutFluxesFromRHS







/** The matrices and vectors are copied to condensed row storage and are then
used by the AMG to find a solution.  
*/
template<size_t dim>
void NodeCenteredFiniteVolumeAlgorithm<dim>::SolveMatrixEquation()
 {

#ifdef CSMP_WITH_SAMG_SOLVER
    if (this->Verbose()){
       cout <<"\n\nNodeCenteredFiniteVolumeAlgorithm<dim>::SolveMatrixEquation: calling solver instance: ";
       cout << GetSolverSettings().GetSolverInstance() << endl;
       cout <<"\tSAMG settings:";
       cout <<"\n\t\tiswit  = " << GetSolverSettings().Get_iswit();
       cout <<"\n\t\titypu  = " << GetSolverSettings().Get_ifirst();
       cout <<"\n\t\tlevelx = " << GetSolverSettings().Get_levelx();
       cout << endl;
    }
#endif

     solver_->Solve( LHS, RHS, RESULT );

 } // end SolveMatrixEquation










/** 

Where the condition flag is not DIRICH, results are mapped from solution
vector<double64> back to finite volume cells after performing a range cheque.  

@section messages Messages 

If the result is out of range at more than 10 nodes, an error message
is returned. If errors occur at more than 2 per cent of the nodes, an
'out_of_range' exception is thrown.  

*/
template<size_t dim>
double64 NodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults( const PropertyDatabase<dim>& p, 
                                                                const csmp::Index& adv_key, 
                                                                bool show_range ,
                                                                const size_t var_comp_nr ) const
 {
    double64  rmin, rmax,
              amin = RESULT[0],
              amax = RESULT[0],
              difference_to_last_output(0.);
    size_t    error_counter(0);
    ScalarVariable   sc;
    
    p.RangeOf( p.Name(adv_key), rmin, rmax );
    if (adv_key.type==SCALAR){
        for ( size_t i=0; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );
            if ( gref_.N(i)->Status( adv_key ) != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                gref_.N(i)->Read( adv_key, sc );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc.Value()) );
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin )
                    gref_.N(i)->Store( adv_key, sc=RESULT[i] );
                else {
                    cout <<"\nNodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults: value: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                    if ( RESULT[i] > rmax )
                        gref_.N(i)->Store( adv_key, sc=rmax );
                    else if ( RESULT[i] < rmin )
                        gref_.N(i)->Store( adv_key, sc=rmin );
                    error_counter++;
                }
            }
        }
    }
    else if (adv_key.type==ARRAY){
        ArrayVariable av;
        for ( size_t i=0; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );

            // array variables have only one flag, hence they are read similar to Scalar variables
            if ( gref_.N(i)->Status( adv_key ) != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                gref_.N(i)->Read( adv_key, av );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-av(var_comp_nr) ));
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin ){
                    av(var_comp_nr)=RESULT[i];
                    gref_.N(i)->Store( adv_key, av );
                }
                else {
                    cout <<"\nNodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults: value: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                    if ( RESULT[i] > rmax ){
                        av(var_comp_nr)=rmax;
                        gref_.N(i)->Store( adv_key, av );
                    }
                    else if ( RESULT[i] < rmin ){
                        av(var_comp_nr)=rmin;
                        gref_.N(i)->Store( adv_key, av );
                    }
                    error_counter++;
                }
            }
        }
    }
    else if (adv_key.type==FLAGGEDARRAY){
        FlaggedArrayVariable fav;
        for ( size_t i=0; i<RESULT.size(); i++ ) {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );
            VARIABLE_FLAG flag(ANY);
            gref_.N(i)->Status( adv_key ,var_comp_nr,flag);
            // Flagged array variables as well as vectors and tensors need a special status call with no return value (it seems)
            if ( flag != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                gref_.N(i)->Read( adv_key, fav );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-fav(var_comp_nr) ));
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin ){
                    fav(var_comp_nr)=RESULT[i];
                    gref_.N(i)->Store( adv_key, fav );
                }
                else {
                    cout <<"\nNodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults: value: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                    if ( RESULT[i] > rmax ){
                        fav(var_comp_nr)=rmax;
                        gref_.N(i)->Store( adv_key, fav );
                    }
                    else if ( RESULT[i] < rmin ){
                        fav(var_comp_nr)=rmin;
                        gref_.N(i)->Store( adv_key, fav );
                    }
                    error_counter++;
                }
            }
        }
    }
    // reporting problems
    if ( error_counter > 20U ) {
          throw csmp::Exception( ERROR, "NodeCenteredFiniteVolumeAlgorithm::OutputResults", 
                                     "Output property was out of range, legal (min/max) was stored instead");
      }
    if ( error_counter > (gref_.Nodes() / 20U) ) 
      throw out_of_range("NodeCenteredFiniteVolumeAlgorithm::OutputResults: Advected variable out of range");

    if ( show_range && this->Verbose()) {
         cout <<"\n\nNodeCenteredFiniteVolumeAlgorithm<"<< dim;
         cout <<">::OutputResults: Variable range after advection: ";
         cout << fixed << setprecision(5) << amin <<" to "<< amax << endl << endl;
      }
      
    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults














/**

Updates saturation water as well! - by default phase1 = wetting phase
multiphase version - assumes range to be between 0 and 1

*/
template<size_t dim>
double64 NodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults( const PropertyDatabase<dim>& p, 
                                                                size_t result_phase,
                                                                const csmp::Index& adv1_key, 
                                                                const csmp::Index& adv2_key,
                                                                bool show_range ) const
 {
    assert( result_phase >  0 );
    assert( result_phase <= 2 );
 
     double64       rmin, rmax,
                    amin = RESULT[0],
                    amax = RESULT[0],
                    difference_to_last_output(0.);
    size_t          error_counter(0);
    ScalarVariable  sc;
    
    if ( result_phase == 1U ) {
      p.RangeOf( p.Name(adv1_key), rmin, rmax );
	    for ( size_t i=0; i<RESULT.size(); i++ ) 
	      {
	         // recording output range
	         amin = std::min( amin, RESULT[i] );
	         amax = std::max( amax, RESULT[i] );
		       if ( gref_.N(i)->Status( adv1_key ) != DIRICH ) {
  		         // reading the pre-existing value and calculating the maximum change per node
  		         gref_.N(i)->Read( adv1_key, sc );
  		         difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc.Value()) );
  		         // result checking
  		         if ( RESULT[i] <= rmax && RESULT[i] >= rmin ) {
  		              gref_.N(i)->Store( adv1_key, sc=RESULT[i] );
  		              gref_.N(i)->Store( adv2_key, sc=1.-RESULT[i] );
  		           }
  		         else {
  		              cout <<"\nNodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults: value: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
  		              if ( RESULT[i] > rmax ) {                           
  		                   gref_.N(i)->Store( adv1_key, sc=rmax );
  		                   gref_.N(i)->Store( adv2_key, sc=1.-rmax );
  		                }
  		              else if ( RESULT[i] < rmin ) {                          
  		                   gref_.N(i)->Store( adv1_key, sc=rmin );
  		                   gref_.N(i)->Store( adv2_key, sc=1.-rmin );
  		                }
  		              error_counter++;
  		           }
  		       }
	      }
      }

   else if ( result_phase == 2U ) {
     p.RangeOf( p.Name(adv2_key), rmin, rmax );
	   for ( size_t i=0; i<RESULT.size(); i++ ) 
	      {
          // recording output range
          amin = std::min( amin, RESULT[i] );
          amax = std::max( amax, RESULT[i] );
          if ( gref_.N(i)->Status( adv2_key ) != DIRICH ) {
             // reading the pre-existing value and calculating the maximum change per node
             gref_.N(i)->Read( adv2_key, sc );
             difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc.Value()) );
             // result checking
             if ( RESULT[i] <= rmax && RESULT[i] >= rmin ) {
                  gref_.N(i)->Store( adv2_key, sc=RESULT[i] );
                  gref_.N(i)->Store( adv1_key, sc=1.-RESULT[i] );
             }
             else {
                 cout <<"\nvalue: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                 if ( RESULT[i] > rmax ) {
                      gref_.N(i)->Store( adv2_key, sc=rmax );
                      gref_.N(i)->Store( adv1_key, sc=1.-rmax );
                 }
                 else if ( RESULT[i] < rmin ) {
                      gref_.N(i)->Store( adv2_key, sc=rmin );
                      gref_.N(i)->Store( adv1_key, sc=1.-rmin );
                     }
                   error_counter++;
                 }
  		       }
	      }
   }
      
    // reporting problems
    if ( error_counter > 20U ) {
          throw csmp::Exception( ERROR, "NodeCenteredFiniteVolumeAlgorithm::OutputResults", 
                                "Output property was out of range, legal (min/max) was stored instead");
      }
    if ( error_counter > (gref_.Nodes() / 20U) ) 
      throw out_of_range("NodeCenteredFiniteVolumeAlgorithm::OutputResults: Advected variable out of range");

    if ( show_range ) {  
         cout <<"\n\nNodeCenteredFiniteVolumeAlgorithm<"<< dim;
         cout <<">::OutputResults: Variable range after advection: ";
         cout << fixed << setprecision(5) << amin <<" to "<< amax << endl << endl;
      }
      
    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults




/**

   @todo JM move these kinds of testing methods into the subclasses that are used for the testing:
     JM addition
*/
template<size_t dim>
double64 NodeCenteredFiniteVolumeAlgorithm<dim>::OutputResultsWithL2NormRes( const PropertyDatabase<dim>& p,
                                                                             const csmp::Index& adv_key,
                                                                             bool show_range ) const
 {
    double64             rmin, rmax,
                         amin = RESULT[0],
                         amax = RESULT[0],
                         difference_to_last_output(0.);
    size_t               error_counter(0);
    ScalarVariable   sc;

    p.RangeOf( p.Name(adv_key), rmin, rmax );

    for ( size_t i=0; i<RESULT.size(); i++ ) {
         // recording output range
         amin = std::min( amin, RESULT[i] );
         amax = std::max( amax, RESULT[i] );
           if ( gref_.N(i)->Status( adv_key ) != DIRICH ) {
             // reading the pre-existing value and calculating the maximum change per node
             gref_.N(i)->Read( adv_key, sc );
             difference_to_last_output += (RESULT[i]-sc.Value())*(RESULT[i]-sc.Value());
             // result checking
             if ( RESULT[i] <= rmax && RESULT[i] >= rmin )
               gref_.N(i)->Store( adv_key, sc=RESULT[i] );
             else {
                  cout <<"\nvalue: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                  if ( RESULT[i] > rmax )
                    gref_.N(i)->Store( adv_key, sc=rmax );
                  else if ( RESULT[i] < rmin )
                    gref_.N(i)->Store( adv_key, sc=rmin );
                  error_counter++;
               }
            }
      }

    // reporting problems
    if ( error_counter > 20U ) {
          throw csmp::Exception( ERROR, "NodeCenteredFiniteVolumeAlgorithm::OutputResultsWithL2NormRes",
                                     "Output property was out of range, legal (min/max) was stored instead");
      }
    if ( error_counter > (gref_.Nodes() / 20U) )
      throw out_of_range("NodeCenteredFiniteVolumeAlgorithm::OutputResultsWithL2NormRes: Advected variable out of range");

    if ( show_range ) {
         cout <<"\n\nNodeCenteredFiniteVolumeAlgorithm<"<< dim;
         cout <<">::OutputResultsWithL2NormRes: Variable range after advection: ";
         cout << scientific << setprecision(5) << amin <<" to "<< amax << endl << endl;
      }

    return sqrt(difference_to_last_output) / std::max( amax - amin, 1.0e-20 );

 } // end OutputResultsWithL2NormRes




  /**

  Save results for saturation in order to use it on the next nonlinear
  iteration step. Updates saturation water as well! - by default phase1 = wetting phase
  multiphase version - assumes range to be between 0 and 1
  
  JM 2013
  @todo Dear SKM, please use svn blame to find out who has created certain functions and contact them directly rather than leaving "blameful" comments on the code itself for everyone to see!!!
  this function was not created by JM 2013.

  */

 template<size_t dim>
 int32 NodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults_NonlinearNewtonRaphson( const PropertyDatabase<dim>& p,
                                                                      size_t result_phase,
                                                                      const csmp::Index& adv1_key,
                                                                      const csmp::Index& adv2_key,
                                                                      bool show_range,
                                                                      std::vector<double64>& DS,
                                                                      std::vector<double64>& SN)
   {
       assert( result_phase >  0 );
       assert( result_phase <= 2 );
       double64       rmin, rmax;
       ScalarVariable  sc;
       if ( result_phase == 1U ) {
           p.RangeOf( p.Name(adv1_key), rmin, rmax );
           for ( size_t i=0; i<RESULT.size(); i++ )
           {
               // recording output range
               gref_.N(i)->Read( adv1_key, sc );
               SN[i]=sc.Value();
               DS[i]=RESULT[i];
           }
       }
       else if ( result_phase == 2U ) {
           p.RangeOf( p.Name(adv2_key), rmin, rmax );
           for ( size_t i=0; i<RESULT.size(); i++ )
           {
               // recording output range
               gref_.N(i)->Read( adv2_key, sc );
               SN[i]=sc.Value();
               DS[i]=RESULT[i];
           }
       }
       return 0;

  } // end OutputResults_NonlinearNewtonRaphson




 /**
    @todo JM move these kinds of testing methods into the subclasses that are used for the testing:

    i.e., class NodeCenteredFiniteVolumeAlgorithm_Test : public NodeCenteredFiniteVolumeAlgorithm { ...

    and do not just sloppily copy everything without cleaning up the parts that are no longer
    used or have changed.

    @todo SMK: Dear SKM (from JM):  This method is needed by the Newton Raphson procedure.
    Please contact me directly if these things are important.  It is NOT a testing method.
    I do my best to not perform any "sloppy" copying.
 */
template<size_t dim>
double64 NodeCenteredFiniteVolumeAlgorithm<dim>::OutputResults_NonlinearNewtonRaphson( const PropertyDatabase<dim>& p,
                                                                 size_t result_phase,
                                                                 const csmp::Index& adv1_key,
                                                                 const csmp::Index& adv2_key,
                                                                 bool show_range )
  {
      assert( result_phase >  0 );
      assert( result_phase <= 2 );

      double64  rmin, rmax,
                amin = numeric_limits<double64>::max(),
                amax = numeric_limits<double64>::min(),
                difference_to_last_output(0.);
      ScalarVariable  sc;

      if ( result_phase == 1U ) {
          p.RangeOf( p.Name(adv1_key), rmin, rmax );
          for ( size_t i=0; i<RESULT.size(); i++ )
          {
              // recording output range
              gref_.N(i)->Read( adv1_key, sc );
              amin = std::min( amin, RESULT[i] );
              amax = std::max( amax, RESULT[i] );
              if ( gref_.N(i)->Status( adv1_key ) != DIRICH ) {
                  // reading the pre-existing value and calculating the maximum change per node
                  difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc.Value()) );
                  // result checking
                  if ( RESULT[i] <= rmax && RESULT[i] >= rmin ) {
                      gref_.N(i)->Store( adv1_key, sc=RESULT[i] );
                      gref_.N(i)->Store( adv2_key, sc=1.-RESULT[i] );
                  }
                  else {
                      if ( RESULT[i] > rmax ) {
                          gref_.N(i)->Store( adv1_key, sc=rmax );
                          gref_.N(i)->Store( adv2_key, sc=1.-rmax );
                      }
                      else if ( RESULT[i] < rmin ) {
                          gref_.N(i)->Store( adv1_key, sc=rmin );
                          gref_.N(i)->Store( adv2_key, sc=1.-rmin );
                      }
                  }
              }
          }
      }
      else if ( result_phase == 2U ) {
          p.RangeOf( p.Name(adv2_key), rmin, rmax );
          for ( size_t i=0; i<RESULT.size(); i++ )
          {
              // recording output range
              gref_.N(i)->Read( adv2_key, sc );
              amin = std::min( amin, RESULT[i] );
              amax = std::max( amax, RESULT[i] );
              if ( gref_.N(i)->Status( adv2_key ) != DIRICH ) {
                  // reading the pre-existing value and calculating the maximum change per node
                  difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc.Value()) );
                  // result checking
                  if ( RESULT[i] <= rmax && RESULT[i] >= rmin ) {
                      gref_.N(i)->Store( adv2_key, sc=RESULT[i] );
                      gref_.N(i)->Store( adv1_key, sc=1.-RESULT[i] );
                  }
                  else {
                      if ( RESULT[i] > rmax ) {
                          gref_.N(i)->Store( adv2_key, sc=rmax );
                          gref_.N(i)->Store( adv1_key, sc=1.-rmax );
                      }
                      else if ( RESULT[i] < rmin ) {
                          gref_.N(i)->Store( adv2_key, sc=rmin );
                          gref_.N(i)->Store( adv1_key, sc=1.-rmin );
                      }
                  }
              }
          }
      }
      if ( show_range ) {
          cout <<"\n\nNodeCenteredFiniteVolumeAlgorithm<"<< typeid(double64).name() <<","<< dim;
          cout <<">::RefreshResults: Variable range after advection: ";
          cout << scientific << setprecision(5) << amin <<" to "<< amax << endl << endl;
      }

      return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults_NonlinearNewtonRaphson



template class NodeCenteredFiniteVolumeAlgorithm<1U>; // see default argument
template class NodeCenteredFiniteVolumeAlgorithm<2U>;
template class NodeCenteredFiniteVolumeAlgorithm<3U>;

} // end namespace csmp
