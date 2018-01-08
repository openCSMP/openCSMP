#include "finiteVolumeAuxiliaryFunctions.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"


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
    double64 U_tilde_f = NVD_Function( xi, U_f, U_c );

    // finding the slope limited value of the advected property at the segment according to eqn. 36
    U_tilde_f = U_tilde_f * psi_hat_d + (1.0 - U_tilde_f) * psi_hat_u;

    return U_tilde_f;
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
get the edge length to calculate the saturation gradient*/
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



} // end namespace csmp
