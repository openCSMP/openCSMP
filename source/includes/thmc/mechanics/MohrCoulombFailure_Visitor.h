/*
 *  MohrCoulombFailure_Visitor.h
 *  Revised and corrected.
 */

#ifndef MOHR_COULOMB_FAILURE_VISITOR_H
#define MOHR_COULOMB_FAILURE_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"
#include "ScalarVariable.h"
#include "StressInvariants.h"

namespace csmp {

class ScalarVariable;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;
template<uint32_t> class Region;
template<uint32_t> class Boundary;

/**
    @brief Evaluates Mohr-Coulomb shear and tensile failure criteria at
    element integration points or element barycentres.

    Reads the Cartesian stress tensor, friction angle, cohesion, tensile
    strength, and Biot coefficient from the model database. Fluid pressure
    is included via the Biot effective stress if defined.

    The smooth Mohr-Coulomb yield envelope (Zienkiewicz Vol. 2, Ch. 4.5.11)
    is used to avoid the corners of the hexagonal surface.

    @par Failure criterion
    Yielding occurs when:

        F_mc = p·sin φ - c·cos φ + q·K(θ) ≥ 0

    where p is mean effective stress (compression positive), q is the von
    Mises deviatoric stress, K(θ) is the smooth Mohr-Coulomb angular factor,
    φ is the friction angle, and c is cohesion.

    Tensile failure occurs when the least principal effective stress σ₃ᵉᶠᶠ
    exceeds the tensile strength T_s (i.e. σ₃ᵉᶠᶠ < -T_s).

    @par Sign convention
    Compressive stresses are positive (geotechnical convention) throughout.
    The constructor parameter @p positive_compressive_stress_convention
    controls whether the input stress tensor follows this convention.
    If false, all stress components are negated before processing.

    @par Output variables
    - "failure"        — F_mc value at each integration point / element
    - "tensile failure" — 1 if tensile failure, 0 otherwise
    
    @attention this visitor does consider the effect of pore pressure on shear and tensile failure

    @author S.K. Matthai
*/
template<uint32_t dim>
class MohrCoulombFailure_Visitor final : public Visitor<dim> {
  public:
    /**
        @brief Constructs the visitor and resolves all database keys.

        @param model      The model to visit.
        @param analysis_input_variables    where stresses, strengths etc. are stored; must all be the same.

        @note For 3D models use this constructor. PlaneAssumption is not
              applicable to 3D stress states.

        @throws csmp::Exception if any required variable is missing or has
                the wrong type or placement.
    */
    explicit MohrCoulombFailure_Visitor( Model<dim>& model,
                                         PLACEMENT analysis_input_variables=ELEMENT_INTEGRATION_POINT );

    /**
        @brief Constructs the visitor for 2D models with explicit plane assumption.
        (named constructor for 2D models)

        @param model      The model to visit.
        @param assumption PLANE_STRESS (σ_zz = 0) or
                          PLANE_STRAIN (σ_zz = ν(σ_xx + σ_yy)).
        @param analysis_input_variables    where stresses, strengths etc. are stored; must all be the same.

        @note This overload is only meaningful for dim == 2. For dim == 3
              use the single-argument constructor.

        @throws csmp::Exception if any required variable is missing or has
                the wrong type or placement.
    */
    static MohrCoulombFailure_Visitor<dim>
    For2D( Model<dim>& model, PlaneAssumption assumption, PLACEMENT analysis_input_variables=ELEMENT_INTEGRATION_POINT );
                                         
    ~MohrCoulombFailure_Visitor() = default;

    void Visit( Element<dim>* ) override final;

//    void Visit( Region<dim>* ) override final;
//    void Visit( Boundary<dim>* ) override final;
//    void Visit( Model<dim>* ) override final;

  private:
    /// private constructor used by For2D
    MohrCoulombFailure_Visitor( Model<dim>&      model,
                                PlaneAssumption  assumption,
                                PLACEMENT        target_placement );

    static void ComputeFailure( const StressInvariants& si,
                                 double                  cohesion,
                                 double                  friction_deg,
                                 double                  tensile_str,
                                 double                  biot_alpha,
                                 double                  pf,
                                 double&                 Fmc,
                                 double&                 F01 ) noexcept;

    /**
        @brief Evaluates failure criteria from a stress tensor and material
        properties at a single point.

        Constructs a StressInvariants object, applies Biot effective stress,
        and evaluates both shear and tensile failure criteria.

        @param stress       Cartesian stress tensor (compression positive).
        @param cohesion     Cohesion c [Pa].
        @param friction_deg Friction angle φ [degrees].
        @param tensile_str  Tensile strength T_s [Pa].
        @param biot_alpha   Biot coefficient α [-].
        @param nu           Poisson's ratio
        @param condition       PLANE_STRAIN/PLANE_STRESS (2D only)
        @param pf           Fluid pressure p_f [Pa].
        @param[out] Fmc     Mohr-Coulomb criterion value (≥ 0 means failure).
        @param[out] F01     greater 0 if tensile failure
    */
    void EvaluateFailure( const TensorVariable<dim>& stress,
                          double                     cohesion,
                          double                     friction_deg,
                          double                     tensile_str,
                          double                     biot_alpha,
                          double                     nu,
                          PlaneAssumption            condition,
                          double                     pf,
                          double&                    Fmc,
                          double&                    F01 ) const noexcept;

    static constexpr bool   verbose_ = true;  ///< set to true if you want extra diagnostics
    PlaneAssumption         plane_assumption_;
    PLACEMENT               analysis_var_placement_;

    csmp::Index  stress_key_;        ///< Cartesian stress tensor
    csmp::Index  pf_key_;            ///< fluid pressure (optional)
    csmp::Index  cohesion_key_;      ///< cohesion [Pa]
    csmp::Index  friction_key_;      ///< friction angle [degrees]
    csmp::Index  tensile_str_key_;   ///< tensile strength [Pa]
    csmp::Index  biot_alpha_key_;    ///< Biot coefficient [-]
    csmp::Index  nu_key_;            ///< Poisson's ratio
    csmp::Index  failure_key_;       ///< F_mc output
    csmp::Index  failure01_key_;     ///< tensile failure flag output
};

// specialisation declarations

/**
    @brief Explicit specialisation of the constructor for 2D models.

    Accepts a PlaneAssumption argument which is not meaningful for
    3D or 1D models and is therefore absent from the primary template
    constructor.

    This specialisation is declared here to prevent implicit instantiation
    of the primary template constructor for dim == 2U before this
    specialisation is seen by the compiler.
*/
template<>
MohrCoulombFailure_Visitor<2U>::MohrCoulombFailure_Visitor(
    Model<2U>&      model,
    PlaneAssumption assumption,
    PLACEMENT       target_placement );

} // end csmp

#endif

