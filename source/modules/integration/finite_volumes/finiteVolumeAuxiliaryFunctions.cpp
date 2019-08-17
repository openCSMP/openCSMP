#include "finiteVolumeAuxiliaryFunctions.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"


using namespace std;

namespace csmp {


/**

Computes the spatially limited variable value from Xi, U_f the value of the
property at the segment and the upstream value of the advected property
U_c.

@return U_tilde_f at segment.

@section application Application

In the method which calculates the isotropic limiter values.
*/
double64 NVD_Function( double64 xi, double64 U_f, double64 U_c )
 {
    constexpr double64 one(1.), zero(0.);
    if ( U_c > one || U_c < zero ) return U_c;
    double64  U_tilde_f = fmax( zero, U_f );
    return fmin( fmin(xi * U_c, U_tilde_f), one );
 }


/**

When looping over the faces, this function computes temporal limiter values
at each face(=segment) in this formulation).

@section arguments Input Arguments

Input are values of the advected property psi at
 u=upwind, c=current, d=downwind FV cells. Notation:

- vol_upstr, vol_dwstr are volumes of control volumes i and j sharing the face
- hf_t0, hf_t1:  	    value of higher order flux (flux * slope-limited advected
                        quantity) across face at time t and t+1
- dt:					time increment value dt
- vol_upstr, vol_dwstr: pore volume represented by the finite volume
- cval_t0, cval_t1: 	value of advected quantity of current FV cell at t and t+dt
- dval_t0, dval_t1: 	value of advected quantity of downstream FV cell at t and t+dt

- method is implemented exactly like in the paper

@return Method returns the temporally limited segment values of the advected
property.

*/
double64 thetaLimiter( double64 hf_t0, double64 hf_t1,
                        double64 vol_upstr, double64 vol_dwstr,
                        double64 psi_hat_c_t1_minus_psi_hat_c_t0,
                        double64 psi_hat_d_t1_minus_psi_hat_d_t0,
                        double64 dt )
 {
    // gf(t+1/2) * dt (Pain et al., eqn 48)
    double64 gf = hf_t0 - hf_t1;
    gf = dt * ( gf < 0.0 ? -1.0 : 1.0) * std::max( std::fabs(gf), 1.0e-30 );

    // for the current CV (Pain et al., 1/eqn 47)
    double64 ufc = (psi_hat_c_t1_minus_psi_hat_c_t0 * vol_upstr) / gf;

    // for the downwind CV (Pain et al., 1/eqn 47)
    double64 ufd = (psi_hat_d_t1_minus_psi_hat_d_t0 * vol_dwstr) / gf;

    const double64 omega(0.5);
    // limiter formula (Pain et al., eqn 46) omega=0.5=Crank-Nicholson, O.K.
    return std::max( 0.5, 1. - omega * std::min(std::fabs(ufc),std::fabs(ufd)) );

 } // end thetaLimiter



/**

Computes limiter values in the case when diffusive fluxes are incorporated
into the higher order transport scheme.

@section arguments Input Arguments

The calculation utilizes the first-order and higher-order gradients of
in the projected variable as projected on the current facet normal in
order to find the limiter value.

@return The method returns the value of the diffusion limiter for the current
face.

@section implementation Implementation

See Pain et al. "Num. transp. meth. for radiation and multi-phase fluid
flow modelling", p. 14-15, eqn. 29 ff.

@section application Application

To avoid oscillations that may arise if diffusive fluxes are not limited
in the higher-order transport scheme.

*/
double64 diffusionLimiter( double64 grad_psi_O1, double64 grad_psi_O2 )
 {
    const double64 gamma1(2.), gamma2(0.5);
    const double64 sf((grad_psi_O1 < 0.) ? -1. : 1.);

    return sf * std::max( gamma2 * sf * grad_psi_O1, std::min( sf * grad_psi_O2, gamma1 * sf * grad_psi_O1 ) );

 } // end DiffusionLimiter


/**

Calculates the isotropic limiter value using the MINMOD scheme
for a slope limiter xi with default value 2.

@section arguments Input Arguments

SMINMAX stores min and max values of the advected property in the
neighborhood of the upstream node, including its own value.

@return U_tilde_f at segment.

@section application Application

In the second-order transport methods of the ExplicitStencilProcessor classes

*/
double64  limitProperty( double64 psi_hat_c, double64 psi_hat_d,
                         double64 psi_facet,
                         const std::pair<double64,double64>& SMINMAX_upstr,
                         double64 xi )
 {
    // isotropic multi-dimensional limiting of facet saturations        min                   max
    const double64 psi_hat_u = (psi_hat_d > psi_hat_c) ? SMINMAX_upstr.first : SMINMAX_upstr.second;

    // finding min, max values of advected property in neighborhood of upstream node
    double64 psi_diff = psi_hat_d - psi_hat_u;
    // avoiding division by zero

    if ( std::fabs(psi_diff) < std::numeric_limits<double64>::epsilon() ) psi_diff = 1.0e-20 * ( psi_diff < 0.0 ? -1.0 : 1.0);

    // computing U_f, U_c factors for slope limiting (Pain et al. paper)
    // 'psi_dash_f' = fsn (advected property values extrapolated to facet integration point)
    const double64 U_f = (psi_facet  - psi_hat_u) / psi_diff;
    const double64 U_c = (psi_hat_c  - psi_hat_u) / psi_diff;

    // computing normalized (always positive) slope-limiter value at face (eqn 38, Pain et al. 2001)
    const double64 U_tilde_f = NVD_Function( xi, U_f, U_c );

    // finding the slope limited value of the advected property at the segment according to eqn. 36
    const double64 psi_hat = lerp(U_tilde_f, psi_hat_u, psi_hat_d);

    return psi_hat;
 }




/**
    Limiter MINMOD?
    
    @author propably Roman Manasipov?
*/
double64  limitProperty( double64 psi_hat_c, double64 psi_hat_u,
                         double64 min, double64 max )
 {
    double64 r(1.);
   
    if ( psi_hat_c > psi_hat_u ) {
        r  = max;
        r -= psi_hat_u;
        r /= ( psi_hat_c - psi_hat_u );
    }
    else if ( psi_hat_c < psi_hat_u ) {
        r  = min;
        r -= psi_hat_u;
        r /= ( psi_hat_c - psi_hat_u );
    }

    return r;
 }




template <size_t dim>
void limitProperty_LSMGRAD( const Element<dim>& e,
                            const csmp::Index& mass_center_key,
                            const csmp::Index& grad_sn_key,
                            const csmp::Index& grad_sn_limiter_key,
                            size_t inside_node, size_t outside_node, size_t iFacet,
                            const double64 sn_inside_node, const double64 sn_outside_node,
                            double64& limited_sn_inside_node, double64& limited_sn_outside_node )
 {

    ScalarVariable limiter_sn_inside_node, limiter_sn_outside_node;
    VectorVariable<dim> grad_sn_inside_node(PLAIN,0.0), grad_sn_outside_node(PLAIN,0.0);
    VectorVariable<dim> distance_inside_node(PLAIN,0.0), distance_outside_node(PLAIN,0.0);
    VectorVariable<dim> mass_center_inside_node(PLAIN,0.0), mass_center_outside_node(PLAIN,0.0);

    // Calculate the coordinates of the Facet Integration Point
    const Point<dim> local_c_point(e.FV()->FacetIntegrationPoint( iFacet, 0U ));
    std::vector<double64> temp(e.Nodes()); //has the local interp. function values
    std::vector<double64> global_c(dim),local_c(local_c_point.Coordinates());

    if( e.IsLineElement() ){
        e.FE()->Nr( local_c[0], temp );
    }else if( e.IsSurfaceElement()){
        e.FE()->Nrs( local_c[0], local_c[1], temp );
    }else{
        e.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
    }

    e.CoordinateMatrix();

    global_c.assign( dim, 0.0);

    // transform local c's to global c's
    for (size_t m = 0; m<e.Nodes(); m++)
        for (size_t n = 0; n<dim; n++){
          global_c[n] += e.FE()->XY(m,n)*temp[m];
    }


    e.N(inside_node)->Read( mass_center_key, mass_center_inside_node);
    e.N(inside_node)->Read( grad_sn_key, grad_sn_inside_node);
    e.N(inside_node)->Read( grad_sn_limiter_key, limiter_sn_inside_node);

    for(size_t l=0U;l<dim;l++)
        distance_inside_node.Component(l,global_c[l]-mass_center_inside_node[l]);

    double64 sn_linear_increment_inside_node=0.0;
    for(size_t l=0U;l<dim;l++)
        sn_linear_increment_inside_node+=grad_sn_inside_node[l]*distance_inside_node[l];

    limited_sn_inside_node = sn_inside_node + limiter_sn_inside_node()*sn_linear_increment_inside_node;


    e.N(outside_node)->Read( mass_center_key,mass_center_outside_node);
    e.N(outside_node)->Read( grad_sn_key, grad_sn_outside_node);
    e.N(outside_node)->Read( grad_sn_limiter_key, limiter_sn_outside_node);

    for(size_t l=0U;l<dim;l++)
        distance_outside_node.Component(l,global_c[l]-mass_center_outside_node[l]);

    double64 sn_linear_increment_outside_node=0.0;
    for(size_t l=0U;l<dim;l++)
        sn_linear_increment_outside_node+=grad_sn_outside_node[l]*distance_outside_node[l];

    limited_sn_outside_node = sn_outside_node  + limiter_sn_outside_node()*sn_linear_increment_outside_node;


 }
 
 
 
/**
    Limiter for gradient based method.
    
    @author propably Roman Manasipov?
*/
template <size_t dim>
double64 limitProperty_LSMGRAD( const Element<dim>& e,
                                const csmp::Index& mass_center_key,
                                const csmp::Index& grad_sn_key,
                                const csmp::Index& grad_sn_limiter_key,
                                size_t upstream_node, size_t iFacet,
                                const double64 sn_upstream_node)
 {

    ScalarVariable limiter_sn_upstream_node;
    VectorVariable<dim> grad_sn_upstream_node(PLAIN,0.0);
    VectorVariable<dim> distance_upstream_node(PLAIN,0.0);
    VectorVariable<dim> mass_center_upstream_node(PLAIN,0.0);

    // Calculate the coordinates of the Facet Integration Point
    const Point<dim> local_c_point(e.FV()->FacetIntegrationPoint( iFacet, 0U ));
    std::vector<double64> temp(e.Nodes()); //has the local interp. function values
    std::vector<double64> global_c(dim),local_c(local_c_point.Coordinates());

    if( e.IsLineElement() ){
        e.FE()->Nr( local_c[0], temp );
    }else if( e.IsSurfaceElement()){
        e.FE()->Nrs( local_c[0], local_c[1], temp );
    }else{
        e.FE()->Nrst( local_c[0], local_c[1], local_c[2], temp );
    }

    e.CoordinateMatrix();

    global_c.assign( dim, 0.0);

    // transform local c's to global c's
    for (size_t m = 0; m<e.Nodes(); m++)
        for (size_t n = 0; n<dim; n++){
          global_c[n] += e.FE()->XY(m,n)*temp[m];
    }


    e.N(upstream_node)->Read( mass_center_key, mass_center_upstream_node);
    e.N(upstream_node)->Read( grad_sn_key, grad_sn_upstream_node);
    e.N(upstream_node)->Read( grad_sn_limiter_key, limiter_sn_upstream_node);

    for(size_t l=0U;l<dim;l++)
        distance_upstream_node.Component(l,global_c[l]-mass_center_upstream_node[l]);

    double64 sn_linear_increment_upstream_node=0.0;
    for(size_t l=0U;l<dim;l++)
        sn_linear_increment_upstream_node+=grad_sn_upstream_node[l]*distance_upstream_node[l];

    return sn_upstream_node + limiter_sn_upstream_node()*sn_linear_increment_upstream_node;

 }




/**

length of edge 01s as described by a vector
get the distance between control volume centres which share this face
(this distance is equivalent to the length of corresponding edge)
get the edge length to calculate the saturation gradient

*/
double64  diffusionVelocity( const Region<1>& sg,
                              const Node<1U>* const nd,
                              const set<size_t>& ngraph_entry,
                              const csmp::Index& advected_var_key )
 {
    double64 diff_flux, max_diff_flux(0.);
    
    for ( set<size_t>::const_iterator it=ngraph_entry.begin(); it!=ngraph_entry.end(); it++ ) 
      {
	     diff_flux = sg.N( *it )->x() - nd->x();
	     diff_flux *= sg.N( *it )->Read( advected_var_key ) - nd->Read( advected_var_key ) ;

         max_diff_flux = std::max( max_diff_flux, fabs(diff_flux) );
      }
    
    return max_diff_flux;
    
 } // end diffusionVelocity



double64  diffusionVelocity( const Region<2>& sg,
                              const Node<2U>* const nd,
                              const set<size_t>& ngraph_entry,
                              const csmp::Index& advected_var_key )
 {
    double64 edge[2], diff_flux, max_diff_flux(0.);
    
    for ( set<size_t>::const_iterator  it=ngraph_entry.begin(); it!=ngraph_entry.end(); it++ ) 
      {
	     // length of edge 01s as described by a vector
	     // get the distance between control volume centres which share this face
	     // (this distance is equivalent to the length of corresponding Triangle edge)
	     edge[0] = sg.N( *it )->x() - nd->x();
	     edge[1] = sg.N( *it )->y() - nd->y();
	     // get the edge length to calculate the saturation gradient
	     diff_flux  = hypot(edge[0],edge[1]);
	     diff_flux *= sg.N( *it )->Read( advected_var_key ) - nd->Read( advected_var_key ) ;

         max_diff_flux = std::max( max_diff_flux, fabs(diff_flux) );
      }
    
    return max_diff_flux;
    
 } // end diffusionVelocity



double64  diffusionVelocity( const Region<3>& sg,
                              const Node<3U>* const nd,
                              const set<size_t>& ngraph_entry,
                              const csmp::Index& advected_var_key )
 {
    double64 edge[3], diff_flux, max_diff_flux(0.);
    
    for ( set<size_t>::const_iterator  it=ngraph_entry.begin(); it!=ngraph_entry.end(); it++ ) 
      {
	     // length of edge 01s as described by a vector
	     // get the distance between control volume centres which share this face
	     // (this distance is equivalent to the length of corresponding Triangle edge)
	     edge[0] = sg.N( *it )->x() - nd->x();
	     edge[1] = sg.N( *it )->y() - nd->y();
	     edge[2] = sg.N( *it )->z() - nd->z();
	     // get the edge length to calculate the saturation gradient
	     diff_flux  = sqrt(edge[0]*edge[0]+edge[1]*edge[1]+edge[2]*edge[2]);
	     diff_flux *= sg.N( *it )->Read( advected_var_key ) - nd->Read( advected_var_key ) ;

         max_diff_flux = std::max( max_diff_flux, fabs(diff_flux) );
      }
    
    return max_diff_flux;
    
 } // end diffusionVelocity

template<size_t dim>
double64 fluxThroughFiniteVolume( Node<dim> const& node, Index const& velocityKey )
{
  double64 totalFlux(0.);
  VectorVariable<dim>  vel;

  // for all those sectors of the FE_FV-stencils which contribute to finite volume
  for ( size_t t(0); t < node.Parents(); ++t ) {
    const size_t nid(node.ParentNodeNumber(t));
    double64  flux(0.);
    // for all facets surrounding the finite volume at the boundary
    for ( size_t i(0); i < node.Parent(t)->FV()->FacetsPerSector(nid); ++i )
    {
      size_t iFacet( node.Parent(t)->FV()->FacetSurroundingSector(nid,i) );
      // fluxes are determined for the sectors inside and outside of the advection region
      node.Parent(t)->Read( velocityKey, vel );
      if ( nid == node.Parent(t)->FV()->InsideNode( iFacet ) )
        flux += ( *node.Parent(t) ).FacetArea(iFacet) * // facet normal velocity
        ( *node.Parent(t) ).ProjectionOnFacetNormal( iFacet, velocityKey );
      else
        flux -= ( *node.Parent(t) ).FacetArea(iFacet) * // facet normal velocity
        ( *node.Parent(t) ).ProjectionOnFacetNormal( iFacet, velocityKey );
    }
    totalFlux += flux; 
  }
  return totalFlux;
}


template<size_t dim, typename ForwardIt>
double64 fluxThroughFiniteVolumes( ForwardIt nodesBegin, ForwardIt nodesEnd, Index const& velocityKey )
{
  double64 totalFlux(0.);

  while( nodesBegin != nodesEnd )
    totalFlux += fluxThroughFiniteVolume( *(*nodesBegin++), velocityKey );

  return totalFlux;
}




/**

Calculates the critical tranport distance in the direction of flow
delta x from the volume of the finite-volume cell under the assumption
that the FV cell is spherical. For this case this distance is equivalent
to the diameter of the sphere d = 2 (3 V / 4 pi)^1/3 where V is the
FV volume.

For skewed elements this method may fail, but it should give a rather
conservative estimate for the case where the element is elongated in
the direction of the flow.

@section arguments Input Arguments

The method needs the volume of the current FV cell as input.

@return The method returns the diameter of the hypothetically spherical
FV cell.
 */
double64 delta_X_FromFV_Volume( double64 FV_volume, size_t dim )
 {
    // cross-section length from volume of a sphere
    if ( dim == 3U ) return 2. * std::pow( (3. * FV_volume) / (4. * PI), 1./3. );
    // diameter of a circle from its area
    if ( dim == 2U ) return 2. * std::sqrt(FV_volume / PI);
    // length (1D case)
    return FV_volume;

 } // end ExplicitNodeCenteredFiniteVolumeTransport




/**

Member function to balance LHS and RHS contributions to the system of equations.
Ascertains that when source term, src, gets negative, the return value omega
goes to zero.

@section arguments Input Arguments

Third argument is the source term that shall be distributed between
left and right.
 */
double64 omega( double64 dt, double64 pore_vol, double64 src )
 {
    // limit the source terms according to equation 61, p. 31 (original manuscript),
    // but max criterion is not applied
    double64 product     = dt * (src / pore_vol);
    double64 e_prod      = std::exp( -product );
    double64 denominator = product * (static_cast<double64>(1.) - e_prod);

    // very important to avoid division by zero
    if ( denominator == static_cast<double64>(0.) ) return static_cast<double64>(1.);

    return (product - static_cast<double64>(1.) + e_prod) / denominator;

 } // end omega







template double64 fluxThroughFiniteVolume( Node<3> const&, Index const& );
template double64 fluxThroughFiniteVolume( Node<2> const&, Index const& );
template double64 fluxThroughFiniteVolume( Node<1> const&, Index const& );



template void limitProperty_LSMGRAD( const Element<1U>&,
                                     const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                     size_t,size_t, size_t,
                                     const double64, const double64,
                                     double64&, double64&);

template void limitProperty_LSMGRAD( const Element<2U>&,
                                     const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                     size_t,size_t, size_t,
                                     const double64, const double64,
                                     double64&, double64&);

template void limitProperty_LSMGRAD( const Element<3U>&,
                                     const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                     size_t,size_t, size_t,
                                     const double64, const double64,
                                     double64&, double64&);

template double64 limitProperty_LSMGRAD( const Element<1U>&,
                                const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                size_t, size_t, const double64);

template double64 limitProperty_LSMGRAD( const Element<2U>&,
                                const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                size_t, size_t, const double64);

template double64 limitProperty_LSMGRAD( const Element<3U>&,
                                const csmp::Index&, const csmp::Index&, const csmp::Index&,
                                size_t, size_t, const double64);



/**
     We store the volume of the finite volumes and their pore volumes
 
     - finite volume
     - finite volume (effective) pore volume
     - sector volume (stored at sector integration point)
     - sector weight: sector pore volume / finite volume pore volume = weighting factor
 
     Those parameters are also computed for the "halo elements" i.e., those elements
     that contribute sectors to the finite volumes that belong to the target region
     but are themselves not part of it.
*/
template<size_t dim>
void initializeFiniteVolumeProperties( Model<dim>& model, Region<dim>& gref, bool initialize_flux )
 {
    // input variables
    const csmp::Index phi_key = model.Database().StorageKey("porosity");
    const csmp::Index thi_key = model.Database().StorageKey("thickness");
    const csmp::Index vt_key  = model.Database().StorageKey("total velocity");
    // output variables
    const csmp::Index fv_key  = model.Database().StorageKey("finite volume");
    const csmp::Index pv_key  = model.Database().StorageKey("FV pore volume");
    const csmp::Index sv_key  = model.Database().StorageKey("sector volume");
    const csmp::Index spv_key = model.Database().StorageKey("sector pore volume");
    const csmp::Index fa_key  = model.Database().StorageKey("facet area");
    const csmp::Index ff_key  = model.Database().StorageKey("facet flux");
    const csmp::Index fn_key  = model.Database().StorageKey("facet normal");
    const csmp::Index fb_key  = model.Database().StorageKey("flux balance");
   
    Point<dim>           nrml;
    VectorVariable<dim>  fnrml, vt;
   
    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "finite volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
    gref.InputPropertyValue( "flux balance", makeScalar(PLAIN,0.), COMPLETE );
   
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // element-based total velocity
         if ( initialize_flux ) (*it)->Read( vt_key, vt );

         // 1. computing sector pore volumes
         // --------------------------------
         // (scaled by the cell thickness attribute=1 for volumetric elements)
         const double64 phi = (*it)->Read( phi_key ) * (*it)->Read( thi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);
              (*it)->Store( i, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
              (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
              // sector volume is added to  pore volume of FV's containing this sector
              double64 finite_volume = (*it)->N(i)->Read( fv_key );
              double64 pore_volume   = (*it)->N(i)->Read( pv_key );
              // sector volume from FV traits
              finite_volume += sector_volume;
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( fv_key, makeScalar(PLAIN,finite_volume) );
              (*it)->N(i)->Store( pv_key, makeScalar(PLAIN,pore_volume) );
           }

         // 2. computing facet normals and areas
         // ------------------------------------
         for ( size_t j=0U; j<facets; ++j ) {
              // computing facet areas
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
              // computing facet normals
              nrml = (*it)->FacetNormal(j);
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, fn_key, fnrml );
           
              // 3. computing total facet fluxes and flux balance
              // ------------------------------------------------
              if ( initialize_flux ) {
                   double64 facet_flux(nrml[0] * vt[0]);
                   if ( dim != 1U ) facet_flux += nrml[1] * vt[1];
                   if ( dim == 3U ) facet_flux += nrml[2] * vt[2];
                   facet_flux *= facet_area;
                   (*it)->Store( j, 0U, ff_key, makeScalar((*it)->Status( j, 0U, ff_key),facet_flux) );
                }
           }
      }

   // 4. initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // -------------------------------------------------------------------------------------------------------------------
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
   
   for ( typename vector<Node<dim>*>::iterator nit=gref.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents());
        for ( size_t i=0U; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // ---------------------------------
             // computing facet normals and areas
             // ---------------------------------
             const size_t facets(eptr->Facets());
             for ( size_t j=0U; j<facets; ++j ) {
                  // computing facet areas
                  const double64 facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
                  // computing facet normals
                  nrml = eptr->FacetNormal(j);
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, fn_key, fnrml );
               }
             // ---------------------------------------
             // computing sector volumes & pore volumes
             // ---------------------------------------
             const double64 porosity = eptr->Read( phi_key ) * eptr->Read( thi_key );
             const size_t sectors(eptr->Sectors());
             for ( size_t j=0U; j<sectors; ++j ) {
                  const double64 sector_volume = eptr->SectorVolume(j);
                  eptr->Store( j, 0U, sv_key, makeScalar( PLAIN, sector_volume ) );
                  eptr->Store( j, 0U, spv_key, makeScalar( PLAIN, sector_volume * porosity ) );
               }
          }
     }
   
   // 5. computing FV flux balances over the complete stencils
   // --------------------------------------------------------
   if ( initialize_flux ) {
        // loop over FV stencils, computing the relevant variable values
        const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
        double64 bmin(1e30), bmax(-1e30);
     
        for ( typename vector<Node<dim>*>::iterator nit=gref.NodesBegin(); nit!=nit_end; ++nit )
          if ( (*nit)->AtBoundary() == NOT )
            {
               const size_t parent_elements((*nit)->Parents());
               double64 flux_balance(0.);
               for ( size_t i=0U; i<parent_elements; ++i ) {
                    const Element<dim>* const eptr = (*nit)->Parent(i);
                    const size_t sector_node = (*nit)->ParentNodeNumber(i);
                    for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
                         const size_t facet  = eptr->FV()->FacetSurroundingSector( sector_node, j );
                         const double64 sign = (sector_node==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
                         const double64 facet_flux = sign * eptr->Read( facet, 0U, ff_key );
                         flux_balance += facet_flux;
                      }
                 }
               (*nit)->Store( fb_key, makeScalar((*nit)->Status(fb_key),flux_balance) );
            
               bmin = std::min( bmin, (*nit)->Read( fb_key ) );
               bmax = std::max( bmax, (*nit)->Read( fb_key ) );
            }
        cout <<"\ninitializeFiniteVolumeProperties: initial flux balance: "<< std::max(fabs(bmin), fabs(bmax)) << endl;
     }
 
 } // end initializeFiniteVolumeProperties

// explicit instantiation of function template in 2 and 3D
template void initializeFiniteVolumeProperties( Model<2U>&, Region<2U>&, bool );
template void initializeFiniteVolumeProperties( Model<3U>&, Region<3U>&, bool );




/**
    Inflow / outflow integration over finte volume sectors.
 
    @note if lower-dimensional elements are part of the mesh, a thickness attribute has to be specified for them,
    this is taken into account when the facet fluxes are first computed.
    If so, thickness is considered by this method.

    @return accumulated influx into the FV sector (added), outflux subtracted.
*/
template<size_t dim>
double64 sectorFlux( const Element<dim>* const eptr, size_t sector, const csmp::Index& flux_key )
 {
   double64 sector_flux(0.);
   for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector); ++j ) {
         const size_t facet  = eptr->FV()->FacetSurroundingSector( sector, j );
         const double64 sign = (sector==eptr->FV()->InsideNode(facet)) ? 1. : -1.;
         const double64 facet_flux = sign * eptr->Read( facet, 0U, flux_key );
         sector_flux += facet_flux;
      }
   return sector_flux;
 
 } // end

template double64 sectorFlux( const Element<2U>* const, size_t, const csmp::Index& );
template double64 sectorFlux( const Element<3U>* const, size_t, const csmp::Index& );




template<size_t dim>
void initializeBasicFiniteVolumeProperties( Model<dim>& model, Region<dim>& gref )
 {
    // input variables
    const csmp::Index phi_key = model.Database().StorageKey("porosity");
    const csmp::Index thi_key = model.Database().StorageKey("thickness");

    // output variables
    const csmp::Index pv_key  = model.Database().StorageKey("FV pore volume");
    const csmp::Index spv_key = model.Database().StorageKey("sector pore volume");
    const csmp::Index fa_key  = model.Database().StorageKey("facet area");
    const csmp::Index fn_key  = model.Database().StorageKey("facet normal");

    Point<dim>           nrml;
    VectorVariable<dim>  fnrml;
 
    // 0. zeroing sector pore volumes for accumulation in element loop
    // ---------------------------------------------------------------
    gref.InputPropertyValue( "FV pore volume", makeScalar(PLAIN,0.), COMPLETE );
 
    // For the interior elements of the region compute relevant variable values
    const typename vector<Element<dim>*>::iterator it_end(gref.ElementsEnd());
    for ( typename vector<Element<dim>*>::iterator it=gref.ElementsBegin(); it!=it_end; ++it )
      {
         const size_t sectors((*it)->Sectors());
         const size_t facets((*it)->Facets());

         // 1. computing sector pore volumes and sector rock compressibilities
         // ------------------------------------------------------------------
         // (scaled by the cell thickness attribute=1 for volumetric elements)
         const double64 phi = (*it)->Read( phi_key ) * (*it)->Read( thi_key );
         for ( size_t i=0U; i<sectors; ++i ) {
              // sector pore volume
              const double64 sector_volume = (*it)->SectorVolume(i);
              (*it)->Store( i, 0U, spv_key, makeScalar( PLAIN, phi * sector_volume ) );
              // sector volume is added to  pore volume of FV's containing this sector
              double64 pore_volume  = (*it)->N(i)->Read( pv_key );
              // sector volume from FV traits
              pore_volume   += phi * sector_volume;
              (*it)->N(i)->Store( pv_key, makeScalar(PLAIN,pore_volume) );
              // accumulating 'total rock compressibility' from the 'compressibility rock' values on the sectors
           }

         // 2. computing facet normals and areas
         // ------------------------------------
         for ( size_t j=0U; j<facets; ++j ) {
              // computing facet areas
              const double64 facet_area = (*it)->FacetArea(j);
              (*it)->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
              // computing facet normals
              nrml = (*it)->FacetNormal(j);
              fnrml(0) = nrml[0];
              if ( dim != 1U ) fnrml(1) = nrml[1];
              if ( dim == 3U ) fnrml(2) = nrml[2];
              (*it)->Store( j, 0U, fn_key, fnrml );
           }
      }

   // 4. initialising facet area, facet normals, sector volume (/pore volume) in the elements surrounding perimeter nodes
   // -------------------------------------------------------------------------------------------------------------------
   // (here the pore volumes do not include the sectors outside the region)
   const typename vector<Node<dim>*>::iterator nit_end(gref.NodesEnd());
 
   for ( typename vector<Node<dim>*>::iterator nit=gref.PerimeterNodesBegin(); nit!=nit_end; ++nit ) {
        const size_t parent_elements((*nit)->Parents());
        for ( size_t i=0U; i<parent_elements; ++i ) {
             Element<dim>* const eptr = (*nit)->Parent(i);
             // ---------------------------------
             // computing facet normals and areas
             // ---------------------------------
             const size_t facets(eptr->Facets());
             for ( size_t j=0U; j<facets; ++j ) {
                  // computing facet areas
                  const double64 facet_area = eptr->FacetArea(j);
                  eptr->Store( j, 0U, fa_key, makeScalar( PLAIN, facet_area ) );
                  // computing facet normals
                  nrml = eptr->FacetNormal(j);
                  fnrml(0) = nrml[0];
                  if ( dim != 1U ) fnrml(1) = nrml[1];
                  if ( dim == 3U ) fnrml(2) = nrml[2];
                  eptr->Store( j, 0U, fn_key, fnrml );
               }
             // ---------------------------------------
             // computing sector volumes & pore volumes
             // ---------------------------------------
             const double64 porosity = eptr->Read( phi_key );
             const size_t sectors(eptr->Sectors());
             for ( size_t j=0U; j<sectors; ++j ) {
                  const double64 sector_volume = eptr->SectorVolume(j);
                  eptr->Store( j, 0U, spv_key, makeScalar( PLAIN, sector_volume * porosity ) );
               }
          }
     }
 
 } // end initializeFiniteVolumeProperties

// explicit instantiation of function template in 2 and 3D
template void initializeBasicFiniteVolumeProperties( Model<2U>&, Region<2U>& );
template void initializeBasicFiniteVolumeProperties( Model<3U>&, Region<3U>& );


} // end namespace csmp
