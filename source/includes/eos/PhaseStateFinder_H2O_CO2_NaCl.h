#ifndef PHASESTATEFINDER_H2O_CO2_NACL_H
#define PHASESTATEFINDER_H2O_CO2_NACL_H

#include "CSMP_definitions.h"
#include "HaliteLiquidus.h"
#include "EOS_CO2H2ONaCl_Spycher2004.h"

namespace csmp {

/**
    Diagnostics on the co-existance of phases captured by an enumeration class.
    'aq' is the default, 'full' allows to calculate a state in the full compositional triangle for testing purposes
    and will not occur in an actual application.
*/
enum class SYSTEM_STATE {aq=0,carb=1,salt=2,aq_salt=3,aq_carb=4,carb_salt=5,aq_carb_salt=6,full=7,undefined=8};

std::string parseState( SYSTEM_STATE );

/**

    @attention the mass fractions of non existing phases are set to zero.
*/
class PhaseStateFinder_H2O_CO2_NaCl
{
    // enum corner uses a convention similar to csmp for the corners of a triangle
    enum corner{a=0,b=1,c=2};
  
    // sorry for all the "publicity"
    // this is a result of the way I developed it, will be cleaned
public:
    PhaseStateFinder_H2O_CO2_NaCl( const double64& t, const double64& p, const double64& total_mass,
                                   const double64& massfraction_h2o, const double64& massfraction_co2, const double64& massfraction_nacl );

    //  Equilibrate is almost all you need - will find the correct  phase state and compute all fluid properties
    //  that come directly from eos
    SYSTEM_STATE Equilibrate();

    // Access to results, always run Equilibrate before
    double64 massCarbonicPhase() const;
    double64 rho_carb() const;
    double64 mu_carb() const;
    double64 beta_carb() const;
    double64 Y_h2o() const;
    double64 Y_co2() const;

    double64 massAqueousPhase() const;
    double64 rho_aq() const;
    double64 mu_aq() const;
    double64 beta_aq() const;
    double64 X_h2o() const;
    double64 X_co2() const;
    double64 X_nacl() const;
    double64 D_Co2() const;

    double64 massHalite() const;

    // i/o for debugging, can later be abandoned
    void DisplayFractions( SYSTEM_STATE );
    void InputComposition();
    void InputXcompYcomp();
  
    SYSTEM_STATE State() const;
  
    /// output the current state in as much as it is of interest to user
    void Out( bool print_private_state_as_well=false ) const;

private:
    PhaseStateFinder_H2O_CO2_NaCl();
    // tracking the outside world
    const double64& t;                                      ///< temperature [Celsius]
    const double64& p;                                      ///< pressure [Pascal]
    const double64& total_mass;                             ///< total mass of H2O+CO2+NaCl in fluid phase and NaCl crystalline salt in the control volume [kg]
    const double64& bulk_massfraction_h2o;                  ///< fraction of total_mass that is H2O [-]
    const double64& bulk_massfraction_co2;                  ///< fraction of total_mass that is CO2 [-]
    const double64& bulk_massfraction_nacl;                 ///< fraction of total mass that is NaCl (dissolved and solid) [-]
  
    // constants used
    static constexpr double64 PI = 3.14159265358979323846;  ///< the number PI
    static constexpr double64 h2o_x = 0.;                   ///< cartesian x-coordinate of the H2O corner of the triangular phase diagram
    static constexpr double64 h2o_y = 0.;                   ///< cartesian y-coordinate of the H2O corner of the triangular phase diagram
    static constexpr double64 co2_x = 1.;                   ///< cartesian x-coordinate of the CO2 corner of the triangular phase diagram
    static constexpr double64 co2_y = 0.;                   ///< cartesian y-coordinate of the CO2 corner of the triangular phase diagram
    static constexpr double64 nacl_x = 0.5;                 ///< cartesian x-coordinate of the NaCl corner of the triangular phase diagram
    const double64 nacl_y = sin(60./180.*PI);               ///< cartesian y-coordinate of the NaCl corner of the triangular phase diagram
    static constexpr double64 molar_mass_h2o  = 18.015e-3;  ///< [kg mol-1]
    static constexpr double64 molar_mass_co2  = 44.010e-3;  ///< [kg mol-1]
    static constexpr double64 molar_mass_nacl = 58.443e-3;  ///< [kg mol-1]

    // private variables
    double64 weight_phase_a,weight_phase_b,weight_phase_c;  ///< weights (mole-based) of phases on corner points in triangular regions
    double64 bulk_xh2o,bulk_xco2,bulk_xnacl;                ///< mole fractions of bulk composition
    double64 x_bulkcomp,y_bulkcomp;                         ///< cartesian coordinates of bulk composition
    double64 x[3],y[3];                                     ///< corner points of triangular regions within phase diagram
    double64 xsat;                                          ///< mole fraction of NaCl in salt-saturated but CO2-free water
    double64 msat;                                          ///< molality of NaCl in salt-saturated aqueous phase (this does not depend on CO2) [mole NaCl / kg H2O]

    double64 massCarbonicPhase_;                            ///< [kg]
    double64 rho_carb_;                                     ///< Density of carbonic phase [kg m-3]
    double64 mu_carb_;                                      ///< dynamic viscosity of carbonic phase [Pa s]
    double64 beta_carb_;                                    ///< isothermal compressibility of carbonic phase [Pa-1]
    double64 Y_h2o_;                                        ///< mass fraction of H2O in carbonic phase [-]
    double64 Y_co2_;                                        ///< mass fraction of CO2 in carbonic phase [-]

    double64 massAqueousPhase_;                             ///< [kg]
    double64 rho_aq_;                                       ///< Density of aqueous phase [kg m-3]
    double64 mu_aq_;                                        ///< dynamic viscosity of aqueous phase [Pa s]
    double64 beta_aq_;                                      ///< isothermal compressibility of aqueous phase [Pa-1]
    double64 X_h2o_;                                        ///< mass fraction of H2O in aqueous phase [-]
    double64 X_co2_;                                        ///< mass fraction of CO2 in aqueous phase [-]
    double64 X_nacl_;                                       ///< mass fraction of NaCl in aqueous phase [-]
    double64 D_Co2_;                                        ///< Diffusivity of CO2 in the aqueous phase (just copy/pasted from SKM's code, no further info)

    double64 massHalite_;                                   ///< [kg]

    SYSTEM_STATE   mystate;
    HaliteLiquidus liquidus;
    EOS_CO2H2ONaCl_Spycher04 eos;

    // mole fraction composition of phases in triangular or lower dimensional order phase regions, order as "corner" enum
    double64 corner_xh2o[3];
    double64 corner_xco2[3];
    double64 corner_xnacl[3];
    // same thing in mass fraction
    double64 corner_massfraction_h2o[3];
    double64 corner_massfraction_co2[3];
    double64 corner_massfraction_nacl[3];

  private:
    // helper functions
    void MoleToMassfractionAtTriangularCorners();
    void ConvertMolefractionTriangularCornersToCartesian(corner/*,double64 xh2o, double64 xco2, double64 xnacl*/);
    void ConvertMolefractionToCartesian(double64& x, double64& y, double64 xh2o, double64 xco2, double64 xnacl);
    bool EvaluateFractionsInTriangularRegion();
    void EvaluateSaltyAqCarb();
    void SetCornerPointsTriangularRegions( SYSTEM_STATE );
};

inline double64 PhaseStateFinder_H2O_CO2_NaCl::massCarbonicPhase() const { return massCarbonicPhase_; }
inline double64 PhaseStateFinder_H2O_CO2_NaCl::rho_carb() const { return rho_carb_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::mu_carb() const { return mu_carb_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::beta_carb() const { return beta_carb_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::Y_h2o() const { return Y_h2o_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::Y_co2() const { return Y_co2_;}

inline double64 PhaseStateFinder_H2O_CO2_NaCl::massAqueousPhase() const { return massAqueousPhase_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::rho_aq() const { return rho_aq_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::mu_aq() const { return mu_aq_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::beta_aq() const { return beta_aq_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::X_h2o() const { return X_h2o_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::X_co2() const { return X_co2_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::X_nacl() const { return X_nacl_;}
inline double64 PhaseStateFinder_H2O_CO2_NaCl::D_Co2() const { return D_Co2_; }

inline SYSTEM_STATE  PhaseStateFinder_H2O_CO2_NaCl::State() const {return mystate;}

inline double64 PhaseStateFinder_H2O_CO2_NaCl::massHalite() const { return massHalite_;}

} // end csmp

#endif // PHASESTATEFINDER_H2O_CO2_NACL_H
