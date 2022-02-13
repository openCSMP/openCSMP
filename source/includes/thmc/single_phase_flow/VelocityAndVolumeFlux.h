#ifndef VELOCITY_AND_VOLUME_FLUX_H
#define VELOCITY_AND_VOLUME_FLUX_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "VectorVariable.h"

namespace csmp {

template<uint32_t> class MathOperatorLHS;
template<uint32_t> class Element;
template<uint32_t> class Model;

/**
 
@brief Post-processes flow "velocity" & "volume flux" that are output as vector<double> and scalar variables, respectively
upon request these properties are extrapolated to the nodes and averaged between adjacent elements.

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
@endcode
 
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
template<uint32_t dim,class CELL=Element<dim> >
class VelocityAndVolumeFlux : public MathOperatorLHS<dim> {
  public:
    VelocityAndVolumeFlux( const Model<dim>&, 
                           const char* oper,             // conductivity
                           const char* basic,            // porosity
                           const char* test,             // fluid pressure
                           bool  node_averaging=false,   // node output 
                           const char* velocity = "velocity",
                           const char* pore_velocity = "pore velocity",
                           const char* volume_flux = "volume flux",
                           const char* nodal_velocity = "nodal velocity",
                           const char* nodal_pore_velocity = "nodal pore velocity",
                           const char* nodal_volume_flux = "nodal volume flux");

    // including gravity (buoyancy effects)
    VelocityAndVolumeFlux( const Model<dim>&, 
                           const char* oper,             // conductivity
                           const char* basic,            // porosity
                           const char* test,             // fluid pressure 
                           const char* relative_density, // relative fluid density 
                           bool  node_averaging=false,   // node output 
                           const char* velocity = "velocity",
                           const char* pore_velocity = "pore velocity",
                           const char* volume_flux = "volume flux",
                           const char* nodal_velocity = "nodal velocity",
                           const char* nodal_pore_velocity = "nodal pore velocity",
                           const char* nodal_volume_flux = "nodal volume flux");

    /// as previous but offering possibility to multiply operand by a distributed variable
    VelocityAndVolumeFlux( const Model<dim>&, 
                           const char* oper,             // conductivity
                           const char* basic,            // porosity
                           const char* test,             // fluid pressure 
                           const char* relative_density, // relative fluid density 
                           const char* prop_multiplier,  // operand multiplier 
                           bool  node_averaging=false,   // node output 
                           const char* velocity = "velocity",
                           const char* pore_velocity = "pore velocity",
                           const char* volume_flux = "volume flux",
                           const char* nodal_velocity = "nodal velocity",
                           const char* nodal_pore_velocity = "nodal pore velocity",
                           const char* nodal_volume_flux = "nodal volume flux");
    
    void Verbose( bool stdoutput );

    virtual void GetOperands( CELL& );
    virtual void WriteOperands( CELL& );

    /// {V} = [grad P]{k}
    virtual void ComputeContribution( CELL& );
    virtual VelocityAndVolumeFlux<dim,CELL>* clone() const { return new VelocityAndVolumeFlux<dim,CELL> (*this); }
  private:

    void ExtractVelocity( const DenseMatrix<DM_MIN>& INP, uint32_t col, VectorVariable<dim>& );
    void ExtractVolumeFlux( const DenseMatrix<DM_MIN>& INP, uint32_t col, ScalarVariable& );
    void ExtractInterstitialVelocity( const DenseMatrix<DM_MIN>& INP, uint32_t col, VectorVariable<dim>& );
    void TestRangeOfOutputVariables() const;
    bool WithLowerDimensionalElements( const Model<dim>& ) const;

    const uint32_t components_;
    
    csmp::Index  velo_key_, ivelo_key_, nvelo_key_, nivelo_key_,
                 flux_key_, nflux_key_, rhor_key_, mult_key_;
               
    std::pair<double,double>     minmaxV_, minmaxF_;
    DenseMatrix<DM_MIN>          DERIV_, RESULT_;
    std::vector<double>          VELOFLUX_, IVELOFLUX_;
    std::vector<double>          IPVF_, NVF_, veloflux_;
    std::vector<double>          IPOL_;
    std::vector<ScalarVariable > PF_, mult_vec_, rho_vec_;
    VectorVariable<dim>          velo_, ivelo_;
    ScalarVariable               phi_, flux_, rhor_;
    bool                         verbose_, nodal_averaging_, with_gravity_, with_multiplier_;
    double                       sum_, ac_gravity_, rho_fac_, mult_fac_;
    std::vector<bool>            node_output_;
    std::vector<std::list<std::vector<double> > >  temp_veloflux_;
    const uint32_t               VERTICAL_AXIS_;
};
 
} // end csmp
 
 
#endif








