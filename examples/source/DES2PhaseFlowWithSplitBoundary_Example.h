// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H
#define DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H

#include "Example.h"
#include "Model.h"
#include "FlowFunctionsModule.h"

namespace csmp {

class  DES2PhaseFlowWithSplitBoundary_Example : public Example
{

public:
    virtual void Run();
    virtual void Specifications();

private:

  template<uint32_t dim>
  void RunSimulation(Model<dim>& model);

  template<uint32_t dim, template<uint32_t> class FLOW_FUNCTIONS>
  void Compute2PhaseFlowProperties( Model<dim>& mdl, FLOW_FUNCTIONS<dim>& flowfunctions, bool with_gravity, bool with_tensor_k );

  template<uint32_t dim>
  void ComputeSteadyStatePressure(Model<dim>& mdl, bool with_gravity, bool with_tensor_k);

};

/// calculates dip vectors for lower-dimensional elements in the model
template<uint32_t dim>
void computeGravityDipVectors( Model<dim>& model );

/**
      loops over the finite volumes of the supplied lower-dimensional supplied region, computing:
      
      -  diameter:        average of half of the distances to the nodes connected in the plane of the lower-dimensional FV
      - normal:            average of the normals of the surface elements making up the sectors of the the FV
      - vertical extent: of FV as determined from the midpoints of the segments connecting the lower-dim FV to its neighbors
      
      @todo the "dip vector" variables and thickness attributes need to be computed elsewhere
*/
template<uint32_t dim>
void computeFV_Diameter_Normal_VerticalExtent( Model<dim>&, const std::string& region_name );


/**
    Computes the CO2 saturation value above which there a wedge of CO2 reaches from the lowest point of the FV to the highest point
*/
template<uint32_t dim>
void computeSpillPointSaturation( Model<dim>&, const std::string& region_name );

/**
      Checks whether thickness is large enough to permit spillage.
      Checks whether thickness is << element diameter, else lower-dimensional representation makes little sense.
*/
template<uint32_t dim>
void lowerDimensionalLayerDiagnostics( Model<dim>&, const std::string& region_name );


} // csmp

#endif // DES_2PHASE_FLOW_WITH_SPLIT_BOUNDARY_EXAMPLE_H
