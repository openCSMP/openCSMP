// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef IMPES_Setup_H
#define IMPES_Setup_H

#include "CSMP_definitions.h"

namespace csmp {

enum RELPERM_MODEL {LINEARMODEL, BROOKSCOREY, VANGENUCHTEN, EXPERIMENTAL, HETEROGENEITY_AWARE };

enum FORMULATION {GLOBAL, PHASE};

enum SOLUTION {IMPES, QEIMPES, QFIMPES};

template<uint32_t> class PropertyDatabase;

template<uint32_t dim>
class IMPES_Setup {

public:
    IMPES_Setup();
    IMPES_Setup(bool with_gravity, bool with_capillary);
    ~IMPES_Setup();

    void WriteVariableFile( const char* file_name );

    void WithCapillaryForces(bool);
    bool WithCapillaryForces() const;

    void WithGravitationalForces(bool);
    bool WithGravitationalForces() const;

    void WithDivergenceCorrection(bool);
    bool WithDivergenceCorrection() const;

    void TwoPhaseFlowModel(RELPERM_MODEL relperm);
    RELPERM_MODEL TwoPhaseFlowModel() const;
    
    void Formulation(FORMULATION formulation);
    FORMULATION Formulation() const;

    void Solution(SOLUTION solution);
    SOLUTION Solution() const;

    void NormalSaturationChanges(double);
    double NormalSaturationChanges() const;
    
    void MaximumSaturationChanges(double);
    double MaximumSaturationChanges() const;

    void CFL_Multiplier(double);
    double CFL_Multiplier() const;

    void MaximumTransportTimesteps(size_t);
    size_t MaximumTransportTimesteps() const;

    void MaximumTransportTimestepSize(double);
    double MaximumTransportTimestepSize() const;

    void ActivationCriteria(double);
    double ActivationCriteria() const;

    void Limiter(double);
    double Limiter() const;

    void Deactivate(size_t);
    size_t Deactivate() const;

    void DefaultSetting();
    void ReadSettingsFromFile(const PropertyDatabase<dim>& database);

#ifdef CSMP_WITH_SAMG_SOLVER
    void SAMG_Reconstruction(size_t);
    size_t SAMG_Reconstruction() const;
#endif

    void ReportSetting();
    
    std::set<std::string>::const_iterator   OutputOnceBegin();
    std::set<std::string>::const_iterator   OutputOnceEnd();

    std::set<std::string>::const_iterator   OutputAlwaysBegin();
    std::set<std::string>::const_iterator   OutputAlwaysEnd();


private:
    bool                 with_capillary_forces_,
                         with_gravitational_forces_,
                         with_divergence_correction_;
    RELPERM_MODEL        relperm_model_;
    FORMULATION          formulation_;
    SOLUTION             solution_;
    double             normal_saturation_changes_;
    double             maximum_saturation_changes_;
    double             CFL_multiplier_;
    size_t               maximum_number_transport_timesteps_;
    double             maximum_transport_timestep_size_;
    double             activation_criteria_;
    double             limiter_;
    size_t               deactivate_;
    size_t               samg_reconstruction_;
    std::set<std::string>      output_once_;
    std::set<std::string>      output_always_;

};

template<uint32_t dim>
inline
void IMPES_Setup<dim>::WithCapillaryForces(bool y_or_n)
 { with_capillary_forces_ = y_or_n; }


template<uint32_t dim>
inline
bool IMPES_Setup<dim>::WithCapillaryForces() const
 { return with_capillary_forces_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::WithGravitationalForces(bool y_or_n)
 { with_gravitational_forces_ = y_or_n; }


template<uint32_t dim>
inline
bool IMPES_Setup<dim>::WithGravitationalForces() const
 { return with_gravitational_forces_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::WithDivergenceCorrection(bool y_or_n)
 { with_divergence_correction_ = y_or_n; }
 
 
template<uint32_t dim>
inline
bool IMPES_Setup<dim>::WithDivergenceCorrection() const
 { return with_divergence_correction_; }
 

template<uint32_t dim>
inline
void IMPES_Setup<dim>::TwoPhaseFlowModel( RELPERM_MODEL relperm )
 { relperm_model_ = relperm; }


template<uint32_t dim>
inline
RELPERM_MODEL IMPES_Setup<dim>::TwoPhaseFlowModel() const
 { return relperm_model_;}


template<uint32_t dim>
inline
void IMPES_Setup<dim>::Formulation(FORMULATION formulation)
 { formulation_ = formulation; }
 

template<uint32_t dim>
inline 
FORMULATION IMPES_Setup<dim>::Formulation() const
 { return formulation_; }
 

template<uint32_t dim>
inline
void IMPES_Setup<dim>::Solution(SOLUTION solution)
{ solution_ = solution; }


template<uint32_t dim>
inline 
SOLUTION IMPES_Setup<dim>::Solution() const
{ return solution_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::NormalSaturationChanges(double normal_saturation_changes)
 { normal_saturation_changes_ = normal_saturation_changes; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::NormalSaturationChanges() const
 { return normal_saturation_changes_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::MaximumSaturationChanges(double maximum_saturation_changes)
 { maximum_saturation_changes_ = maximum_saturation_changes; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::MaximumSaturationChanges() const
{ return maximum_saturation_changes_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::CFL_Multiplier(double CFL_multiplier)
 { CFL_multiplier_ = CFL_multiplier; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::CFL_Multiplier() const
 { return CFL_multiplier_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::MaximumTransportTimesteps(size_t maximum_number_transport_timesteps)
 { maximum_number_transport_timesteps_ = maximum_number_transport_timesteps; }


template<uint32_t dim>
inline
size_t IMPES_Setup<dim>::MaximumTransportTimesteps() const
 { return maximum_number_transport_timesteps_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::MaximumTransportTimestepSize(double maximum_transport_timestep_size)
 { maximum_transport_timestep_size_ = maximum_transport_timestep_size; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::MaximumTransportTimestepSize() const
 { return maximum_transport_timestep_size_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::ActivationCriteria(double criteria)
{ activation_criteria_ = criteria; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::ActivationCriteria() const
{ return activation_criteria_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::Limiter(double limiter)
{ limiter_ = limiter; }


template<uint32_t dim>
inline
double IMPES_Setup<dim>::Limiter() const
{ return limiter_; }


template<uint32_t dim>
inline
void IMPES_Setup<dim>::Deactivate(size_t number_of_transport_step)
{ deactivate_ = number_of_transport_step; }


template<uint32_t dim>
inline
size_t IMPES_Setup<dim>::Deactivate() const
{ return deactivate_; }


#ifdef CSMP_WITH_SAMG_SOLVER
template<uint32_t dim>
inline
void IMPES_Setup<dim>::SAMG_Reconstruction(size_t number_of_pressure_step)
{ samg_reconstruction_ = number_of_pressure_step; }

template<uint32_t dim>
inline
size_t IMPES_Setup<dim>::SAMG_Reconstruction() const
{ return samg_reconstruction_; }
#endif

template<uint32_t dim>
inline
std::set<std::string>::const_iterator IMPES_Setup<dim>::OutputOnceBegin()
 { return output_once_.begin(); }


template<uint32_t dim>
inline
std::set<std::string>::const_iterator IMPES_Setup<dim>::OutputOnceEnd()
 { return output_once_.end(); }


template<uint32_t dim>
inline
std::set<std::string>::const_iterator IMPES_Setup<dim>::OutputAlwaysBegin()
 { return output_always_.begin(); }


template<uint32_t dim>
inline
std::set<std::string>::const_iterator IMPES_Setup<dim>::OutputAlwaysEnd()
 { return output_always_.end(); }


 
} // end csmp

#endif
