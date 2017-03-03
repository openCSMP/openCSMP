#include "UpwindControlVisitor.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<size_t dim>
    UpwindControlVisitor<dim>::UpwindControlVisitor( Model<dim>& model, 
                                  ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid,
                                  ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor,
                                  const char* permeability,
                                  const char* porosity,
                                  std::vector<std::string>& densities,
                                  std::vector<std::string>& relperm_vis,
                                  std::vector<std::string>& saturations,
                                  std::vector<std::string>& cfl_variables,
                                  std::vector<std::string>& velocities,
                                  std::vector<std::string>& pore_velocities)
  : fv_transport_vapor( fv_vapor ),
    fv_transport_liquid( fv_liquid ),
    phases(densities.size()),
    gravity(-9.80665), // scalar acts to increase the pressure
    pot_crit(1.e-20),
    xyz(dim),
    recalculate(false),
    flipping(false),
    grav(true),
    with_velocity(false),
    largest_time_step(60.*60.*24.*365.),
    VERTICAL_AXIS( (dim==1u) ? 0u : 1u ),
    cfl_scaling(0.1),
    cfl_with_pore_velocity(true)
 { 

    if (dim==3U) xyz = 2;

    this->ApplicationLevel(REGION);
    this->ApplicationTarget(ELEMENT);

    k_key = model.Database().StorageKey(permeability);
    phi_key = model.Database().StorageKey(porosity);
    KgradP_key = model.Database().StorageKey("KgradP");
    tot_pore_vel_key = model.Database().StorageKey("pore velocity");
    tot_vel_key = model.Database().StorageKey("velocity");
    p_key = model.Database().StorageKey("fluid pressure");
    sh_key = model.Database().StorageKey("saturation halite");

    facet_pore_velocity.resize(phases);
    pore_vel_key.resize( phases );
    pore_phase_velocity.resize(phases);
    vel_key.resize( phases );
    phase_velocity.resize(phases);

    facet_normal.resize(dim);

    if ( k_key.type != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                     permeability, " must be a scalar property." );

    if ( phi_key.type != SCALAR )
      throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                     porosity, " must be a scalar property." );

    if ( densities.size() != phases )
       throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                      " densities-vector has wrong size: " );

    if ( relperm_vis.size() != phases )
       throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                      " relperm-visc-vector has wrong size: " );

    if ( saturations.size() != phases )
       throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)",
                      " saturations-vector has wrong size: " );

    rho_key.resize( phases );
    rho.resize( phases ); 
    cfl_key.resize( phases );
    cfl.resize( phases ); 

    relperm_visc_key.resize( phases );
    relperm_visc.resize( phases ); 

    S_key.resize( phases );
    S.resize( phases );

    size_t i;
    for ( i = 0; i < phases; i++)
       {
        rho_key[i]  = model.Database().StorageKey(densities[i].c_str());
        relperm_visc_key[i]  = model.Database().StorageKey(relperm_vis[i].c_str());
        S_key[i]  = model.Database().StorageKey(saturations[i].c_str());
        cfl_key[i]  = model.Database().StorageKey(cfl_variables[i].c_str());
        vel_key[i]  = model.Database().StorageKey(velocities[i].c_str());
        pore_vel_key[i]  = model.Database().StorageKey(pore_velocities[i].c_str());

        if ( rho_key[i].place != NODE || rho_key[i].type != SCALAR )
          throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                      densities[i].c_str(), " must be a scalar property placed on the nodes." );

        if ( relperm_visc_key[i].place != NODE || relperm_visc_key[i].type != SCALAR )
          throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)", 
                      relperm_vis[i].c_str(), " must be a scalar property placed on the nodes." );

        if ( S_key[i].place != NODE || S_key[i].type != SCALAR )
          throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor::(constructor)",
                      saturations[i].c_str(), " must be a scalar property placed on the nodes." );

       }

     typename std::vector<Element<dim>* >::const_iterator it;

    Upwinder.resize( phases );
    
    for ( i = 0; i < Upwinder.size(); i++)
      {
      Upwinder.at(i).resize( model.Region("Model").Elements() );
      for ( it = model.Region("Model").ElementsBegin(); it < model.Region("Model").ElementsEnd(); it++)
        {
          Upwinder.at(i).at( (*(*it)).Idx() ).Resize((*(*it)).Nodes(),(*(*it)).Nodes());
          Upwinder.at(i).at(( *(*it)).Idx() ) = 0.0;
        }
      }
    
   uc_key.resize( phases );
   uc_key[0]  = model.Database().StorageKey("upwind control liquid");
   uc_key[1]  = model.Database().StorageKey("upwind control vapor");

  }

/** default destructor */
template<size_t dim>
UpwindControlVisitor<dim>::~UpwindControlVisitor() 
 {}


/** visit function for region */
template<size_t dim>
void UpwindControlVisitor<dim>::Visit(Region<dim>* n)
{


}

/** visit function for elements */
template<size_t dim>
void UpwindControlVisitor<dim>::Visit(Element<dim>* n)   
  {

	facets = n->FV()->Facets();

    n->Read( k_key, k );
    n->Read( KgradP_key, KgradP );
    n->NodePropertyVector( sh_key , sh );

    for (size_t i = 0; i < phases; i++)
      {
         n->NodePropertyVector( rho_key[i], rho[i] );
         n->NodePropertyVector( relperm_visc_key[i], relperm_visc[i] );
         n->NodePropertyVector( S_key[i], S[i] );
      }
      
    if (with_velocity)
      {
      n->Read( phi_key, phi );
      for (size_t i = 0; i < phases; i++)
        {
         cfl[i] = largest_time_step;      
         facet_pore_velocity[i].resize(facets);
        }
      }

     DetermineUpwindNodes( *n );

	for (size_t p = 0; p < phases; ++p) {
     uc_scal() = 0;
     for ( size_t i=0U; i<n->FV()->Facets(); i++ )
       {
        n->FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );
        uc_scal() += Upwinder[p][n->Idx()](inside_node_,outside_node_);
       }
     n->Store( uc_key[p], uc_scal );
     n->Store( cfl_key[p], cfl[p] );
    }
    
   
   if (with_velocity)
     {
     for (size_t p = 0; p < phases; ++p)
      {
        n->Store( vel_key[p], phase_velocity[p] );
        n->Store( pore_vel_key[p], pore_phase_velocity[p] );     
      }
      n->Store( tot_pore_vel_key, total_pore_velocity );
      n->Store( tot_vel_key, total_velocity );
      }     

}
     

/** define upwind nodes for CVFEM scheme - and calculate velocitites if requested */
template<size_t dim>
void UpwindControlVisitor<dim>::DetermineUpwindNodes( Element<dim>& e )
{

   velocity  = 0.0;
   decision1 = 0.0;
   decision2 = 0.0;
   normal_component = 0.0;
   g = 0.0;
   total_pore_velocity = 0.0;
   total_velocity = 0.0;
   
  for (size_t p = 0; p < phases; ++p) {
  
    norm     = 0.0;
    mobility = 0.0;
    density  = 0.0;
    sat = 0.0;

    // loop over facets of the element
    for ( size_t i=0U; i<facets; i++ )
       {
           facet_cfl = true;
           e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

         if (p==0)
            velocity = fv_transport_liquid.GetFacetNormalVelocity( e.Idx(), i );
         else
            velocity = fv_transport_vapor.GetFacetNormalVelocity( e.Idx(), i );
         decision1 = decision2 = velocity;
         
         if (grav)
    	    {
    	      normal_component = fv_transport_vapor.GetFacetNormalComponent( e.Idx(), i, xyz-1 );
    	      g = gravity*normal_component;
    		  decision1  += k()*rho[p][inside_node_]()*g;
    		  decision2  += k()*rho[p][outside_node_]()*g;
            }
   
         if (decision1 < 0) decision1 = 0.0;
         if (decision2 > 0) decision2 = 0.0;

    	  if (rho[p][inside_node_]() <= 0)
    		 decision1  = 0.0;
    	  if (rho[p][outside_node_]() <= 0)
    	     decision2  = 0.0;

    	 decision1 *= relperm_visc[p][inside_node_]();
    	 decision2 *= relperm_visc[p][outside_node_]();

    	  if (decision1 > 0 && decision2 < 0)
    	     {
    		   if (abs(decision1) <= abs(decision2)) decision1 = 0.0;
    		   if (abs(decision1) >= abs(decision2)) decision2 = 0.0;
    	     }

    	  if (decision1 > 0 && decision2 >= 0)
    	     {
    		   velocity = decision1;
    		   if (!flipping)
    		     {
    		      if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 1)
    		        recalculate = true;
    		      Upwinder[p][e.Idx()](inside_node_,outside_node_) = 1;
    		      Upwinder[p][e.Idx()](outside_node_,inside_node_) = 2;
    		     }
               else
                 if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 10)
                   {
    		         if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 1)
    		           {
    		            recalculate = true;
    		            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 10;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 10;
    		            velocity = 0.;
    		           }
    		         else
    		           {
                        Upwinder[p][e.Idx()](inside_node_,outside_node_) = 1;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 2;
    		           }
    		       }
    		   abs_velocity = abs(velocity);
    		   density  += rho[p][inside_node_]()*abs_velocity;
    		   mobility += relperm_visc[p][inside_node_]()*abs_velocity;
    		   norm     += abs_velocity;
            sat = S[p][inside_node_]();
            if (e.N(outside_node_)->Status(p_key) == DIRICH)
              facet_cfl = false;
             }
          else if (decision1 <= 0 && decision2 < 0)
    	     {
    		   velocity = decision2;
    		   if (!flipping)
    		     {
    		      if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 2)
    		         recalculate = true;
    		      Upwinder[p][e.Idx()](inside_node_,outside_node_) = 2;
    		      Upwinder[p][e.Idx()](outside_node_,inside_node_) = 1;
    		     }
               else
                 if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 10)
                   {
    		         if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 2)
    		           {
    		            recalculate = true;
    		            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 10;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 10;
    		            velocity = 0.;
    		           }
    		         else
    		           {
    		            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 2;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 1;
    		           }
    		       }
    		   abs_velocity = abs(velocity);
    		   density  += rho[p][outside_node_]()*abs_velocity;
    		   mobility += relperm_visc[p][outside_node_]()*abs_velocity;
    		   norm     += abs_velocity;
            sat = S[p][outside_node_]();
            if (e.N(inside_node_)->Status(p_key) == DIRICH)
              facet_cfl = false;
             }
    	   else
    	     {
    	       velocity = 0.;
    		   if (!flipping)
    		     {
    		      if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 0)
    		         recalculate = true;
    		      Upwinder[p][e.Idx()](inside_node_,outside_node_) = 0;
    		      Upwinder[p][e.Idx()](outside_node_,inside_node_) = 0;
    		     }
               else
                 if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 10)
                   {
    		         if (Upwinder[p][e.Idx()](inside_node_,outside_node_) != 0)
    		           {
    		            recalculate = true;
    		            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 10;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 10;
    		           }
    		         else
    		           {
    		            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 0;
    		            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 0;
    		           }
                }
               facet_cfl = false;
    	      }

          if (sh[inside_node_]() >= 1.0)
            {
            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 2;
            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 1;
            }
          if (sh[outside_node_]() >= 1.0)
            {
            Upwinder[p][e.Idx()](outside_node_,inside_node_) = 2;
            Upwinder[p][e.Idx()](inside_node_,outside_node_) = 1;
            }


     if (with_velocity)
       {
//         if (phi() != 0.) facet_pore_velocity[p][i] = velocity / phi();
         if (sat*phi() != 0.) facet_pore_velocity[p][i] = velocity / (sat*phi());
         else facet_pore_velocity[p][i] = 0.;
         distance = e.N(inside_node_)->Coordinate().DistanceTo(
                                e.N(outside_node_)->Coordinate());

         if (facet_pore_velocity[p][i] != 0. && facet_cfl)
           {
           if (cfl_with_pore_velocity)
              cfl[p]() = std::min(cfl[p](),distance*cfl_scaling/fabs(facet_pore_velocity[p][i]));
           else
              cfl[p]() = std::min(cfl[p](),distance*cfl_scaling/(fabs(facet_pore_velocity[p][i])*phi()));
           }
       }

   }
   
   if (with_velocity)
     {
       phase_velocity[p] = KgradP;
       if ( grav && norm != 0.) phase_velocity[p](VERTICAL_AXIS) += gravity * density * k() / norm;
       phase_velocity[p] *= mobility;
       if (norm != 0.) phase_velocity[p] /= norm;
       if (phi()!= 0.) pore_phase_velocity[p] = phase_velocity[p]/phi();
       
       total_velocity      += phase_velocity[p];
       total_pore_velocity += pore_phase_velocity[p];

     }
   }
    
} // end DetermineUpwindNodes


/** reset booleans */
template<size_t dim>
void UpwindControlVisitor<dim>::Reset()
    {
      recalculate   = false;
      flipping      = false;
      with_velocity = false;
    } // end Reset

/** set boolean recalculate */
template<size_t dim>
bool UpwindControlVisitor<dim>::Recalculate()
    {
      return recalculate;
    } // end Recalculate

/** set boolean flip */
template<size_t dim>
void UpwindControlVisitor<dim>::Flipping(bool flip)
    {
      flipping = flip;
    } // end Flipping

/** access to upwind matrix of specified phase (density) and element */
template<size_t dim>
DenseMatrix<DM_MIN>  UpwindControlVisitor<dim>::UpwindMatrix(csmp::Index rho_index, size_t eidx)
    {
      for (size_t i = 0; i < rho_key.size(); i++)
         if (rho_index == rho_key[i])
            return Upwinder[i][eidx];
      
      // SKM FIX FOR XCODE
      throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor<dim>::UpwindMatrix:", "density index never identified." );
      return DenseMatrix<DM_MIN>();
            
    } // end UpwindMatrix

//// access to vector of upwind matrices of specified phase (density); TODO: refactor: creates DenseMatrix as temporary!
template<size_t dim>
std::vector<DenseMatrix<DM_MIN> >& UpwindControlVisitor<dim>::UpwindMatrices(csmp::Index rho_index)
    {
      for (size_t i = 0; i < rho_key.size(); i++)
         if (rho_index == rho_key[i])
            return Upwinder[i];
            
      // SKM FIX FOR XCODE
      throw csmp::Exception( CSMP_ERROR, "UpwindControlVisitor<dim>::UpwindMatrices:", "density index never identified." );
      return Upwinder[0];

    } // end UpwindMatrices
  
  
  
/** set recalculate boolean */
template<size_t dim>
void UpwindControlVisitor<dim>::Recalculate(bool recalc)
    {
      recalculate = recalc;
    } // end Recalculate

/** activate or deactivate gravity component */
template<size_t dim>
void UpwindControlVisitor<dim>::Gravity(bool with_gravity)
    {
      grav = with_gravity;
    } // end Gravity

/** set boolean for velocity calculations */
template<size_t dim>
void UpwindControlVisitor<dim>::WithVelocity(bool velo)
    {
      with_velocity = velo;
      // change JPW
      flipping      = true;
    } // end Gravity

/** set maximum size of time step */
template<size_t dim>
void UpwindControlVisitor<dim>::SetLargestTimeStep(double64 timestep)
    {
      largest_time_step = timestep;
    } // end SetLargestTimeStep
 
/** modify calculation of cfl criterion */
template<size_t dim>
void UpwindControlVisitor<dim>::Adjust_CFL_Criterion(double64 scale_factor, bool take_pore_velocity)
    {
      cfl_scaling = scale_factor;
      cfl_with_pore_velocity = take_pore_velocity;
    } // end Adjust_CFL_Criterion

template class UpwindControlVisitor<1U>;
template class UpwindControlVisitor<2U>;
template class UpwindControlVisitor<3U>;

} // csmp






















