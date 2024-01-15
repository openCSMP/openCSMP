#ifndef TWO_PHASE_VELOCITY_AND_VOLUME_FLUX_H
#define TWO_PHASE_VELOCITY_AND_VOLUME_FLUX_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "VectorVariable.h"

namespace csmp {

template<uint32_t,template<uint32_t> class> class MathOperatorLHS;
template<uint32_t> class Element;
template<uint32_t> class Node;
template<uint32_t> class Face;
template<uint32_t> class Model;
template<uint32_t> class TwoPhaseModel;


/// flow "velocity" & "volume flux" are output as vector<double> and scalar variables, respecitively
/// upon request these properties are extrapolated to the nodes and averaged between adjacent elements
template<uint32_t dim, template<uint32_t> class CELL>
class TwoPhaseVelocityAndVolumeFlux : public MathOperatorLHS<dim,CELL> {
  public:
    TwoPhaseVelocityAndVolumeFlux( const Model<dim>&,
                           TwoPhaseModel<dim>& satfunc,
                           const char* oper,                // conductivity
                           const char* basic,               // porosity
                           const char* test,                // fluid pressure
                           const char* density_nw,          // density of non-wetting phase
                           const char* density_w,           // density of wetting phase
                           const char* gravity_proj,        // gravity projection
                           const char* gravity_term,        // gravity term
                           bool  with_gravity               = false,    // gravity forces
                           bool  with_capillary_spreading   = false,    // capillary effect
                           bool  phase_velocities           = false,    // phase velocities
                           bool  node_averaging             = false,    // node output
                           const char* thickness            = NULL,
                           const char* prop_multiplier      = NULL,     // operand multiplier
                           const char* velocity             = "velocity",
                           const char* pore_velocity        = "pore velocity",
                           const char* volume_flux          = "volume flux",
                           const char* velocity_nw          = "velocity oil",
                           const char* velocity_w           = "velocity water",
                           const char* nodal_velocity       = "nodal velocity",
                           const char* nodal_pore_velocity  = "nodal pore velocity",
                           const char* nodal_volume_flux    = "nodal volume flux",
                           const char* nodal_velocity_nw    = "nodal velocity oil",
                           const char* nodal_velocity_w     = "nodal velocity water");
                           
    virtual ~TwoPhaseVelocityAndVolumeFlux() {}
    
    void Verbose( bool stdoutput );

    virtual void GetOperands( const CELL<dim>& );

    /// {V} = [grad P]{k}
    virtual void ComputeContribution( const CELL<dim>& );

    virtual void ComputeContribution( const Node<dim>& );

    virtual void WriteOperands( CELL<dim>& );

    
  private:

    void ExtractVelocity( const DenseMatrix<DM_MIN>& INP, uint32_t col, VectorVariable<dim>& );
    void ExtractVolumeFlux( const DenseMatrix<DM_MIN>& INP, uint32_t col, ScalarVariable& );
    void ExtractInterstitialVelocity( const DenseMatrix<DM_MIN>& INP, uint32_t col, VectorVariable<dim>& );
    void TestRangeOfOutputVariables() const;

    // Compute Gravity Term and Total mobility
    void ComputeTotalMobilityRelativeDensityAndGravityTerm( CELL<dim>& );

    bool                            verbose_,
                                    nodal_averaging_,
                                    phase_velocities_,
                                    with_gravity_, with_capillary_, with_multiplier_, multiply_with_cell_thickness_;

    std::vector<bool>               node_output_;

    TwoPhaseModel<dim>&             satFunc_;

    const uint32_t                  components_;
    const uint32_t                  VERTICAL_AXIS_;

    csmp::Index                     velo_key_, ivelo_key_, flux_key_, velo_nw_key_, velo_w_key_,
                                    nvelo_key_, nivelo_key_, nflux_key_, nvelo_nw_key_, nvelo_w_key_,
                                    rho_nw_key_,rho_w_key_,
                                    grav_key_, gravity_term_key_,
                                    thi_key_,
                                    mult_key_;
    csmp::Index                     facet_normal_idx_,              // Facet Normal property index
                                    facet_area_idx_,                // Facet Area property index
                                    sector_volume_idx_;             // Sector Volume property index
    std::pair<double,double>        minmaxV_, minmaxF_;
    DenseMatrix<DM_MIN>             DERIV_, RESULT_;
    std::vector<double>             VELOFLUX_, IVELOFLUX_, VELOFLUX_W_, VELOFLUX_NW_,
                                    IPVF_, NVF_, veloflux_,
                                    IPOL_;

    std::vector<ScalarVariable >    PF_,
                                    mult_vec_,
                                    rho_w_vec_, rho_nw_vec_;

    std::vector<std::list<std::vector<double> > >  temp_veloflux_;

    VectorVariable<dim>             vt_, ivelo_, velo_nw_, velo_w_, gproj_;
    VectorVariable<dim>             facet_n_;
    ScalarVariable                  phi_, flux_, rhor_, rhot_, rho_w_, rho_nw_;
    ScalarVariable                  sc_;
    DenseMatrix<DM_MIN>             DN_;
    Point<dim>                      dsdn_;
    double                          sum_, ac_gravity_, gravTerm_,
                                    rho_fac_, rhot_fac_, rho_w_fac_,
                                    rho_nw_fac_, mult_fac_, cell_thickness_;

};

 
/**
 
@class TwoPhaseVelocityAndVolumeFlux
@author S.K. Matthai
@date 2005

@section motivation Motivation

Post-Processing: Compute flow velocities and their magnitudes from the 
fluid pressure field, the conductivity, and the porosity of the medium.
The flow velocities are computed at the element centers or integration
points but they can be extrapolated to the nodes if this option is 
chosen. The variable names which are hardwired into this math operator 
are:

@code
"velocity"              element     vector
"pore velocity"			element		vector
"volume flux"           element		scalar
"nodal velocity"        node        vector
"nodal pore velocity"	node        vector
"nodal volume flux"		element		scalar
@code
 
Those with a the "nodal" prefix are used only when the option for 
output to the nodes is chosen.  

 
@section implementation Implementation
 
Procedure: In 2 post-processing application cycles of this class 
object the following is performed: 

1. From the nodal "fluid pressure" and the shape function derivatives, 
   the "fluid pressure gradient" at the element center (for elements with 
   linear interpolation functions), or at the integration points is
   computed. For higher order elements these locations are optimal locations, 
   see Barlow, 1977).  

2. The product of the fluid pressure gradient, the conductivity, and the
   porosity (in case of the interstitial velocity) and -1.0, gives the
   flow velocity. The volume flux, as a scalar property, is given by 
   the magnitude of the Darcy velocity vector.  

3. If the properties shall be mapped to the nodes, they are linearly 
   extrapolated to these from the element center or the integration points.
   In the case of linear elements, nodal velocities represent the 
   averages of adjacent elements. These are not weighted by element 
   area !
    
   For quadratic elements, bilinear extrapolation is carried out
   by the specific finite element used. Only linear variations in the velocity
   over the element are supported at this stage. The nodal velocities are 
   averaged among the adjacent elements since they diverge at these locations
   from one-another.  
   
Steps 1 and 2 are accomplished in the first application cycle of the 
post-processing operator. The nodal averaging and following operations
are done in a second cycle.  
*/

 
 
} // end csmp
 
 
#endif








