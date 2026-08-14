//
//  FaultFlowPropertyCalculator.h
//
//  Porosity - permeability relations for fractures and faults under stress.
//
//  Created by Stephan Matthai on 12/18/12.
//  Copyright (c) 2012 Stephan Matthai. All rights reserved.
//
#ifndef FAULT_FLOW_PROPERTY_CALCULATOR_H
#define FAULT_FLOW_PROPERTY_CALCULATOR_H

#include "CSMP_definitions.h"
#include "InSituStress.h"

namespace csmp {

double  background_k_ManningIngebritsen99( double subsurface_depth );

class Index;
struct MechanicalProperties;
template<uint32_t> class Element;
template<uint32_t> class Model;
 
/** Developed for 3D this method will also perform on 2D models, this calculator also aims to also handle 2D fracture models, and
    it includes procedures for fracture aperture calculation
*/
template<uint32_t dim>
class FaultFlowPropertyCalculator {
  public:
      FaultFlowPropertyCalculator() = default;
      ~FaultFlowPropertyCalculator() = default;
    
    /// quick and dirty variation of fluid pressure with depth for the mechanical calculations
    void CalculateLinearFluidPressure( Model<dim>&, double pf_model_top, double fluid_density );
    
    /// Evaluate shear-, normal- and effective stress acting on fractures, then evaluate failure criteria
    void EvaluateStressesAndFailurePotential( Model<dim>&,
                                              const std::set<std::string>& fractures,
                                              const InSituStress& insitu_stres );
    
    /// Evaluates failure criteria and phi-max-limited fault dilation from far-field stress, pf, and burial depth
    void EvaluateFailureAndDilatation( Model<dim>&,
                                       const std::set<std::string>& faults,
                                       const MechanicalProperties& faulted_rock_props,
                                       const InSituStress& insitu_stress,
                                       bool recompute_hydrostatic_fluid_pressure=false );
    
    /// semi-analytic aperture model using maximum fracture dimension to compute constant value
    void ConstantFractureApertureFromCruikshankModel( Model<dim>&, const std::set<std::string>& fractures );

    /// semi-analytic aperture model computing values from distance to fracture tip and centerline
    void VariableFractureApertureFromCruikshankModel( Model<dim>&, const std::set<std::string>& fractures );
    
    /// semi-analytic aperture model in which only the most remote tip of the fracture is closed
    void VariableApertureForElongatedFractures( Model<dim>&, const std::set<std::string>& fractures );
    
    /// compacts fractures where s_eff >= UCS, to mimimum aperture value
    void CompressiveFailureAperture( Model<dim>&,const std::set<std::string>& fractures, double minimum_aperture );
    
    /// calculate fault permeability from distance to tip, fault size, effective stress, porosity and dilatation
    void PermeabilityFromStressAndDilatation( Model<dim>&,
                                              const std::set<std::string>& faults,
                                              double max_aperture,
                                              double flow_tortuosity=1.5,
                                              bool compute_halo_permeability=false );
   
    /// contouring faults with distance from tipline and centerline attributes
    void CalculateTipAndCenterLineDistances( Model<dim>&, const std::set<std::string>& faults );
    
    /// failure criteria: evaluates tensile and shear failure as well as dilation vs. sliding criteria; returns cumulative area of slip patches
    double IsCriticallyStressed( Model<dim>&, const MechanicalProperties&, const char* fault, const InSituStress& );
    
    private:
     /// calculates shortest distance of any element barycenter to fault tip and returns it
     double  DistanceToPerimeter( Model<dim>& model,
                                  const char* fault,
                                  const char* distance_to_tipline_variable );

     /// computes shortest distance of any given point to the barycenter of a given region
     double  DistanceFromCenter( Model<dim>& model,
                                 const char* fault_region,
                                 const char* distance_to_tipline_variable,
                                 const char* distance_from_center_variable );
    
     /// Nodal values of the VectorVariable 'displacement' are interpolated across fault assuming that u=v=0 at the tipline
     bool InterpolateDisplacement() { throw std::logic_error("FaultFlowPropertyCalculator::InterpolateDisplacement: method not implemented yet!"); }
    
     /// extrapolates 'fracture permeability' and 'breccia zone permeability' along normal into damage zone for given fault width, taking into account cumulative displacement
     void CalculateLateralPermeabilityVariation( Element<dim>* surface_e, double breccia_width, double damage_zone_width, const csmp::Index& );
 };

 } // end csmp

#endif /* FAULT_FLOW_PROPERTY_CALCULATOR_H */
