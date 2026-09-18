// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Tractions_Example.h
//  CSMP_API_library2025
//
//  Created by Stephan Matthai on 31/8/2025.
//

#ifndef CSMP_TRACTIONS_EXAMPLE_H
#define CSMP_TRACTIONS_EXAMPLE_H

#include "CSMP_definitions.h"
#include "Example.h"
#include <Eigen/Dense>

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Point;
template<uint32_t> class Face;

/// return type for data from traction calculation
template<uint32_t dim>
struct ForceAndTorque {
    Eigen::Matrix<double,dim,1> force  = Eigen::Matrix<double,dim,1>::Zero();
    Eigen::Matrix<double,dim,1> torque = Eigen::Matrix<double,dim,1>::Zero();
};

struct TractionResult {
    Eigen::Vector3d velocity;     // sample velocity
    Eigen::Matrix3d gradU;        // estimated velocity gradient
    Eigen::Matrix3d viscousStress; 
    Eigen::Matrix3d cauchyStress;
    Eigen::Vector3d traction;     // total traction
    Eigen::Vector3d tractionNormal;
    Eigen::Vector3d tractionTangential;
};


class  Tractions_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();

  private:
  
    // Verification steps in developing this formulation using CSMP++ and Eigen
    bool Test_NormalStressCalculation();
    bool Test_ViscousShearStressCalculation();
    bool Test_TractionsCalculation();
    bool Test_ComputationOfDragForces();
    bool Test_NoSlipTangentialTractionComputation();

    /// calculation example illustrating the no-slip BC and tensor algebra involved
    void DragForceOnInclinedPlane();
    
    /// demo: first rank (linear U-gradient) based shear tractions acting on inclined surface for no-slip boundary conditions
    TractionResult FlowInducedTangentialTraction( const Eigen::Vector3d& u_sample, ///< flow velocity at point above surface
                                                  const Eigen::Vector3d& n,        ///< surface normal
                                                  double d,                        ///< distance from surface along n
                                                  double mu,                       ///< dynamic viscosity
                                                  double p );
};


/// Face wise integration of normal and shear stresses as well as contributions to torque from boundary
template<uint32_t dim>
ForceAndTorque<dim> integrateFluidTractionOnFaceUsingCellBarycentre(
                                            const Face<dim>* face,
                                            const csmp::Index& vel_key,    ///< from volumetric element
                                            const csmp::Index& press_key,  ///< from Face  on the Boundary
                                            double viscosity,
                                            const Eigen::Matrix<double,dim,1>& obj_centroid );


/// Only those tractions that arise from non-face-parallel velocity gradients 
template<uint32_t dim>
ForceAndTorque<dim> computeNormalAndShearTractionsOnBoundary( const Model<dim>&, const std::string& region_boundary, double viscosity );



/**
    Computes normal and shear tractions acting on the user defined Boundary in the model as the consequence of the flow of a Newtonian fluid and pressure gradients acting on the boundary surface.
 */
template<uint32_t dim>
void computeBoundaryTractions( Model<dim>&, std::string boundary );

/// returns pressure force due to pressure gradient acting on 'object'
std::array<double,3u> pressureForceActingOnObject( const Model<3u>&, const std::string& object );




} // csmp

#endif // CSMP_TRACTIONS_EXAMPLE_H
