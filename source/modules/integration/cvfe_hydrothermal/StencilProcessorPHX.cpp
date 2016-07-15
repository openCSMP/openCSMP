#include "StencilProcessorPHX.h"
#include "finiteVolumeAuxiliaryFunctions.h"

using namespace std;

namespace csmp {


/** custom destructor */
template<size_t dim>
StencilProcessorPHX<dim>::StencilProcessorPHX( )
 {
   std::vector<double64> g;
   g.resize(dim,0.0);
   if (dim==1) g[0] = -9.80665;
   else g[1] = -9.80665;
   VectorVariable<dim> gvv( g );
   gravity_ = gvv;   
 } // end StencilProcessorPHX


/** default destructor */
template<size_t dim>
StencilProcessorPHX<dim>::~StencilProcessorPHX()
 {
 } // end destructor


/** Calculating fluxes out of the control volumes, including gravity */
template<size_t dim>
void StencilProcessorPHX<dim>::DetermineFluxOutWithGravity( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              DenseMatrix<DM_MIN>& upwind,
                                              std::vector<std::vector<double64> >& facet_flux,
                                              std::vector<double64>& flux_out,
                                              csmp::Index rhs_key,
                                              csmp::Index rho_key,
                                              csmp::Index k_key)
  {

   eidx_ = e.Idx();
   e.Read( k_key, perm_ );

   for ( size_t i=0U; i<e.FV_Stencil()->Facets(); i++ )
      {
 
        grav_ = param.FacetNormalProjection(i,gravity_);
        e.FV_Stencil()->FacetEdgeNodes( i, inside_node_, outside_node_ );
        
    	if (upwind(inside_node_,outside_node_) == 1)
    	  {
            e.N(inside_node_)->Read( rho_key, density_ );
            g_component_ = grav_*density_()*perm_();
    	    vel_ =  param.FacetNormalVelocity(i) + g_component_;
    	    if (vel_ < 0.0)
    	       vel_ = 0.0;
            e.N(inside_node_)->Read( rhs_key, rhs_property_ );
    	  }
        else if (upwind(inside_node_,outside_node_) == 2)
    	  {
            e.N(outside_node_)->Read( rho_key, density_ );
            g_component_ = grav_*density_()*perm_();
    	    vel_ =  param.FacetNormalVelocity(i) + g_component_;
    	    if ( vel_ > 0.0 )
    	      vel_ = 0.0;
            e.N(outside_node_)->Read( rhs_key, rhs_property_ );
    	  }
        else vel_ =0.0;

        if ( vel_ != 0.0 )
          {
             facet_flux[eidx_][i] = vel_ * param.FacetArea(i) * rhs_property_();
             if (facet_flux[eidx_][i]<0)
                flux_out[ e.N(outside_node_)->Idx() ] -= facet_flux[eidx_][i];
             else
                flux_out[ e.N( inside_node_)->Idx() ] += facet_flux[eidx_][i];
          }
       }
       
  } // end DetermineFluxOutWithGravity


/** Calculating fluxes out of the control volumes without a gravity component */
template<size_t dim>
void StencilProcessorPHX<dim>::DetermineFluxOutWithoutGravity( const FV_Parameter& param,
                                              const Element<dim>& e,
                                              DenseMatrix<DM_MIN>& upwind,
                                              std::vector<std::vector<double64> >& facet_flux,
                                              std::vector<double64>& flux_out,
                                              csmp::Index rhs_key)
  {

   eidx_ = e.Idx();
   for ( size_t i=0U; i<e.FV_Stencil()->Facets(); i++ )
      {
 
       e.FV_Stencil()->FacetEdgeNodes( i, inside_node_, outside_node_ );
       vel_ =  param.FacetNormalVelocity(i);
        
    	if (upwind(inside_node_,outside_node_) == 1)
    	  {
    	    if (vel_ < 0.0) vel_ = 0.0;
              e.N( inside_node_ )->Read( rhs_key, rhs_property_ );
    	  }
        else if (upwind(inside_node_,outside_node_) == 2)
    	  {
    	    if ( vel_ > 0.0 ) vel_ = 0.0;
              e.N( outside_node_ )->Read( rhs_key, rhs_property_ );
    	  }
        else vel_ = 0.0;

        if ( vel_ != 0.0 )
          {
             facet_flux[eidx_][i] = vel_ * param.FacetArea(i) * rhs_property_();             
             if (facet_flux[eidx_][i]<0)
                flux_out[ e.N( outside_node_ )->Idx() ] -= facet_flux[eidx_][i];
             else
                flux_out[ e.N( inside_node_ )->Idx() ] += facet_flux[eidx_][i];
          }
   }
         
  } // end DetermineFluxOutWithoutGravity  

/** Calculating fluxes into the control volumes */
template<size_t dim>
void StencilProcessorPHX<dim>::DetermineFluxIn( const Element<dim>& e,
                                              std::vector<std::vector<double64> >& facet_flux,
                                              std::vector<double64>& flux_in,
                                              std::vector<double64>& mass_balance )
  {
  
   eidx_ = e.Idx();

   for ( size_t i=0U; i<e.FV_Stencil()->Facets(); i++ )
      {
 
        e.FV_Stencil()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        if (facet_flux[eidx_][i]<0)
          flux_in [ e.N( inside_node_ )->Idx() ] += facet_flux[eidx_][i]*mass_balance[ e.N( outside_node_ )->Idx() ];
        else
          flux_in [ e.N( outside_node_ )->Idx() ] -= facet_flux[eidx_][i]*mass_balance[ e.N( inside_node_ )->Idx() ];
		      
      }
  } // end DetermineFluxIn

/** Calculating fluxes out of the control volumes without pre-defined upwind nodes */
template<size_t dim>
void StencilProcessorPHX<dim>::DetermineFluxOut(const FV_Parameter& param,
                                              const Element<dim>& e,
                                              std::vector<std::vector<double64> >& facet_flux,
                                              std::vector<double64>& flux_out,
                                              csmp::Index rhs_key )
  {

  const double64 zero(0.);
  ScalarVariable rhs_property;
  eidx_ = e.Idx();
  
     for ( size_t i=0U; i<e.FV_Stencil()->Facets(); i++ )
       {
 
          e.FV_Stencil()->FacetEdgeNodes( i, inside_node_, outside_node_ );

          if ( param.FacetNormalVelocity(i) != zero ) {

               const double64 n_x_A(param.FacetNormalVelocity(i) * param.FacetArea(i));
 	           // identifying the upstream node and initialize variables for it
               const size_t  nidx( (param.FacetNormalVelocity(i) < zero) ? outside_node_ : inside_node_ );

             if (rhs_key.place == NODE) e.N(nidx)->Read( rhs_key, rhs_property );
             else if (rhs_key.place == ELEMENT) e.Read( rhs_key, rhs_property );
             facet_flux[eidx_][i] = n_x_A * rhs_property();
             
             if (facet_flux[eidx_][i]<0)
                flux_out[ e.N(outside_node_)->Idx() ] -= facet_flux[eidx_][i];
             else
                flux_out[ e.N( inside_node_)->Idx() ] += facet_flux[eidx_][i];

          }
       }

} // end DetermineFluxOut

template class StencilProcessorPHX<1U>;
template class StencilProcessorPHX<2U>;
template class StencilProcessorPHX<3U>;


} // namespace csmp



















