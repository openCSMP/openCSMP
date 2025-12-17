#ifndef UPWIND_CONTROL_VISITOR_H
#define UPWIND_CONTROL_VISITOR_H

#include "Visitor.h"
#include "Model.h"
#include "ExplicitFiniteVolumeTransportPHX.h"

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++ and strong simplification as compared to the csmp5-version.
*/

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class ExplicitFiniteVolumeTransportPHX;


template<uint32_t dim>
class UpwindControlVisitor : public Visitor<dim> {
  public:
    UpwindControlVisitor( Model<dim>& model, 
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_vapor, // specialized FV transport class for vapor phase
                          ExplicitFiniteVolumeTransportPHX<dim>& fv_liquid, // specialized FV transport class for liquid phase
                          const char* permeability, // element permeability
                          const char* porosity, // element porosity for pore velocity calculations
                          std::vector<std::string>& densities, // vapor and liquid densities
                          std::vector<std::string>& relperm_vis, // relative permeabilities for vapor and liquid
                          std::vector<std::string>& saturations, // vapor and liquid saturations
                          std::vector<std::string>& cfl_variables, // vapor and lquid cfl values
                          std::vector<std::string>& velocities, // vapor and liquid velocities
                          std::vector<std::string>& pore_velocities); // vapor and liquid pore velocities

    virtual ~UpwindControlVisitor();
    
    virtual void Visit(Element<dim>* );
	  virtual void Visit(Region<dim>* );

    void Reset(); // reset all bopleans
    bool Recalculate(); // get boolean recalculate
    void Flipping( bool flip ); // set boolean flip
    void Recalculate( bool recalc ); // ser boolean recalculate 
    void Gravity( bool with_gravity ); // activate or deactivate gravity component
    void WithVelocity( bool with_velocity ); // activate or deactive velocity calculations
    void SetLargestTimeStep( double timestep ); // change maximum time step size

    void Adjust_CFL_Criterion( double scale_factor, bool take_pore_velocity); // modify cfl criterion
    
    DenseMatrix<DM_MIN> UpwindMatrix(csmp::Index rho_index, size_t eidx); // access function to specific upwind matrix
    std::vector<DenseMatrix<DM_MIN> >& UpwindMatrices(csmp::Index rho_index); // access function to specific vector of upwind matrices
    
  private:
    void DetermineUpwindNodes( Element<dim>& e );

    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_vapor;
    ExplicitFiniteVolumeTransportPHX<dim>& fv_transport_liquid;
    
    std::vector<csmp::Index> rho_key, relperm_visc_key, S_key;
    std::vector<csmp::Index> uc_key, cfl_key, pore_vel_key, vel_key;
    csmp::Index k_key, phi_key, KgradP_key, tot_pore_vel_key, tot_vel_key, p_key, sh_key;

    ScalarVariable                              k, uc_scal, phi;
    std::vector<ScalarVariable>                 sh;
    VectorVariable<dim>                         KgradP, total_pore_velocity, total_velocity;
    std::vector<VectorVariable<dim> >           pore_phase_velocity, phase_velocity;
    std::vector<std::vector<ScalarVariable> >   rho, relperm_visc, S;

    std::vector<std::vector<DenseMatrix<DM_MIN> > >  Upwinder;
    DenseMatrix<DM_MIN> facet_normals, projection, grad, inversed_facet_normals;

    std::vector<ScalarVariable> cfl;

    std::vector<std::vector<double> >  facet_pore_velocity;

    std::vector<double>  facet_normal;
    
    uint32_t                  VERTICAL_AXIS;
    double gravity, pot_crit;
    double decision1, decision2;
    double velocity, abs_velocity, norm, mobility, density, sat;
    double normal_component, g;
    double distance, largest_time_step;
    double cfl_scaling;

    uint32_t  xyz;       // 1=x, 2=y, 3=z
    uint32_t  inside_node_, outside_node_;
    
    bool recalculate, flipping, grav, with_velocity, facet_cfl;
    bool cfl_with_pore_velocity;

    uint32_t phases, facets; 

};

   /**
     @class UpwindControlVisitor UpwindControlVisitor.h

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      Specialized visitor to define upwind nodes for CVFEM scheme.

     @section usage Usage
      To be used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 Calculates upwind nodes and store information in a matrix.
	 Calculates velocity and pore velocity vectors for graphical representations only!
          
     @endcode
     
     @section dependencies Dependencies
	 The functions are tailored for use in ExplicitFiniteVolumeTransportPHX.
	 This visitor is needed in TwoPhaseTransportPHX and Upwind-PDE-operators.
     
     @section issues Known issues
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

 } // csmp

#endif
