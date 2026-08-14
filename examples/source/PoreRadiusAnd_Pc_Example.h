#ifndef CSMP_PORE_RADIUS_AND_PC_EXAMPLE_H
#define CSMP_PORE_RADIUS_AND_PC_EXAMPLE_H

#include "Example.h"
#include "Box.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Region;

class  PoreRadiusAnd_Pc_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};



void scaleModel( Model<3>&, double scale_factor );

template<uint32_t dim>
double equivalentPermeability( Model<dim>&, const char* region,
                               BOX_BOUNDARY inflow_boundary, BOX_BOUNDARY outflow_boundary,
                               double& flux_through_model );

double boundaryFluxFEM( const Model<3U>& model, const char* boundary );

template<uint32_t dim>
double boxBoundaryFluxFEM( const Region<dim>& flow_domain, const Index& velo_key, BOX_BOUNDARY target_boundary );

void scaleModel( Model<3>& mdl, double scale_factor );

} // csmp

#endif // CSMP_PORE_RADIUS_AND_PC_EXAMPLE_H
