//
//  PVTX_Calculator_H2O_CO2_NaCl.h
//  CSMP-geo-sequestration-simulator
//
//  Created by Stephan Matthai on 2/1/19.
//  core thermodynamics by Thomas Driesner (ETHZ).
//
//  Copyright © 2019 CSMP Originators Group, Stephan Matthai. All rights reserved.
//

#ifndef CSMP_PVTX_CACULATOR_H2O_CO2_NACL_H
#define CSMP_PVTX_CACULATOR_H2O_CO2_NACL_H

#include "VariableSet_CO2GeoSequestration.h"
#include "PhaseStateFinder_H2O_CO2_NaCl.h"
#include "ArrayVariable.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;
template<size_t> class Model;

enum AQUEOUS_PHASE : std::int8_t { XH2O=0, XCO2=1, XNACl_aq=2 }; ///< mass fractions stored
enum CARBONIC_PHASE : std::int8_t { YCO2=0, YH2O=1 };

/// phases that are continuous across element, note extra NEITHER state required by contrast to SYSTEM_STATE; neither corresponds to full & undefined
enum PHASES_CONTINUOUS_ACROSS_ELEMENT : std::int8_t { AQUEOUS, CARBONIC, HALITE, AQUEOUS_HALITE, AQUEOUS_CARBONIC,
                                                      CARBONIC_HALITE, AQUEOUS_CARBONIC_HALITE, NEITHER };

/// returns phases which are present on all nodes of the element so that their physical properties can be interpolated
template<size_t dim>
PHASES_CONTINUOUS_ACROSS_ELEMENT  continuousPhases( const variables::VariableSet_CO2GeoSequestration&, const Element<dim>* );

/// reports the PhaseStateFinder's system state inferred from continuity of the phases across the element
SYSTEM_STATE  elementState( PHASES_CONTINUOUS_ACROSS_ELEMENT );


/**
   @note uses 2 fitted functions to establish relations between aqueous and carbonic phase as function of salinity.
 
   @todo lookup tables not considered yet.
   @todo salt only case fails
*/
template<size_t dim>
class PVTX_Calculator_H2O_CO2_NaCl {
  public:
    explicit PVTX_Calculator_H2O_CO2_NaCl( const variables::VariableSet_CO2GeoSequestration& );
    ~PVTX_Calculator_H2O_CO2_NaCl() {}
  
    /**
        (1) Computes input mass fractions (Xbulk) for Equilibrate() method from saturations, phase compositions, halite volume fraction and fluid densities:
            Inputs: saturations, halite (volume fraction of pore space), mass_NaCl_aq (kg/m3), Pf_, ToC_, rho_aq, rho_carb.
            Outputs: The results are written to the "bulk mass and fractions" array variable with the entries totalmass, Xaq, Xcarb and Xsalt.
     
            @attention use this method before calling Equilibrate() when the model is first initialised.
     */
    void InitialisePVTX_FromFieldData( const Node<dim>* );
  
     /**
        (2) Determines phase stabilities in the triangular compositional diagram with the 3 end-member compositions:
              aqueous - carbonic - salt.
            Computes mass_aq, mass_carb, sw, the mass fractions of the phases YCO2, YH2O, XCO2, XH2O, XNaCl, and stores the results
            in aq_ph_composition_, carb_ph_composition_, and the variable (NaCl_cr = crystalline salt as a volume)
            Proceeds with the computation of the physical properties of the fluids: density, viscosity and compressibility.
     
            @return the phase state that the system is in which is useful for the simulator as a whole and the DES evolution predictions in particular.
    */
    SYSTEM_STATE Equilibrate( Node<dim>* );

    /// element version with properties interpolated to the element barycentre
    // TODO: would it be better if this returns PHASES_CONTINUOUS_ACROSS_ELEMENT
    SYSTEM_STATE Equilibrate( Element<dim>* );

    // accessors for the Element calculation where result properties cannot be stored on the element
  
    /// @return density of aqueous phase (kg/m3) or NaN if the phase is not present
    double64 AqueousPhaseDensity() const;
  
    /// @return density of carbonic phase (kg/m3) or NaN if the phase is not present
    double64 CarbonicPhaseDensity() const;
  
    /// @return viscosity of aqueous phase (Pa.s) or NaN if the phase is not present
    double64 AqueousPhaseViscosity() const;
  
    /// @return viscosity of carbonic phase (Pa.s) or NaN if the phase is not present
    double64 CarbonicPhaseViscosity() const;
  
    /// @return compressibility of aqueous phase (1/Pa) or NaN if the phase is not present
    double64 AqueousPhaseCompressibility() const;

    /// @return compressibility of carbonic phase (1/Pa) or NaN if the phase is not present
    double64 CarbonicPhaseCompressibility() const;
  
    /// computes water saturation from the mass balance; @attention use only if all 3 components are there
    double64 WaterSaturation() const;
  
    /// saturation of halite in the presence of water and carbonic phase else NaN is returned
    double64 HaliteVolumeFraction() const;

  private:
    // copy construction is forbidden
    PVTX_Calculator_H2O_CO2_NaCl( const PVTX_Calculator_H2O_CO2_NaCl& ) = delete;
    PVTX_Calculator_H2O_CO2_NaCl() = delete;
  
    /// as in Node version but with mass interpolated to barycentre of the element and properties calculated there
    void InterpolateInputVariablesToBaryCenter( const Element<dim>* );

    /// finds the highest capillary entry pressure among the elements sharing the node
    double64 MaxEntryPressureOfParentElements( const Node<dim>* ) const;

  private:
    // constants
    const variables::VariableSet_CO2GeoSequestration& props_;
    static constexpr double64 rhoNaCl_ = 2170.;             ///< density = 2170 kg/m3
    // variables that PhaseStateFinder_H2O_CO2_NaCl has constant references to
    double64       Pf_, ToC_;
    enum           XBULK { TMASS=0, XAQ=1, XCARB=2, XSALT=3 };
    ArrayVariable  Xbulk_;                                  ///< total_mass and Xaq, Xcarb, Xsalt mass fractions of the unit-volume system
    // class state variables and variables to avoid copying
    std::vector<double64>  IPOL_;                           ///< for the interpolation of properties
    ArrayVariable aq_ph_composition_, carb_ph_composition_; ///< mass fractions: aqueous phase with 3 vs. carbonic phase with 2 components
    // computational engine for phase-stability calculations; contains EOS
    PhaseStateFinder_H2O_CO2_NaCl flash_;                   ///< including EOS and HaliteLiquidus
};


// Member functions

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::AqueousPhaseDensity() const { return flash_.rho_aq(); }

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::CarbonicPhaseDensity() const { return flash_.rho_carb(); }

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::AqueousPhaseViscosity() const { return flash_.mu_aq(); }

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::CarbonicPhaseViscosity() const { return flash_.mu_carb(); }

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::AqueousPhaseCompressibility() const { return flash_.beta_aq(); }

template<size_t dim>
inline double64 PVTX_Calculator_H2O_CO2_NaCl<dim>::CarbonicPhaseCompressibility() const { return flash_.beta_carb(); }

} // end csmp


#endif /* CSMP_PVTX_CACULATOR_H2O_CO2_NACL_H */
