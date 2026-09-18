// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_STAGGERED_GRID_STOKES_SOLVER_EXAMPLE_H
#define CSMP_STAGGERED_GRID_STOKES_SOLVER_EXAMPLE_H

#include "Example.h"
#include "Box.h"
#include "ModelSubDomain.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Region;

/** @brief StaggeredGridStokesSolver
     Since the Stokes equation is linear, it can be solved directly in a two step process:
     1) compute pressure on a fine grid using LFEM (linear finite element method)
     2) compute velocity on a coarser but node-matched grid using QFEM (quadratic finite element method)
     This example illustrates this consistent (ability to resolve parabolic flow profile within a single quadratic FEM) and therefore highly accurate approach which does not require stabilisation by comparison with Tutorial4_Example
*/
class  StaggeredGridStokesSolver_Example : public Example {
  public:
    virtual void Run();
    virtual void Specifications();
  
  private:
    void DisplayTitle();

    void AssignNoSlipBoundaryConditions( Model<3>&,
                                         Region<3>& stokes_flow_domain,
                                         Region<3>& grain_edges,
                                         size_t flow_direction );

    void ConstructVelocityVector( Model<3>& );

    double VolumeAveragedVelocity( const Model<3U>&, const std::string& flow_domain, double& velovity_integral );

    double EquivalentPermeabilityQuadraticFEM( Model<3U>&, const std::string& flow_domain,
                                               BOX_BOUNDARY inflow_boundary,
                                               BOX_BOUNDARY outflow_boundary,
                                               double& flux_through_model );
      
    /// compute seepage forces and fluid induced tractions acting on pore walls / solid skeleton
    void SeepageForces( Model<3U>& linear_FEM_model );
    void SeepageTractions( Model<3U>& quadratic_FEM_model );

    void PrintAnalysisResults( const Model<3U>&, double porosity, const std::vector<double>&  k_equiv );

    // TESTING
    /// quadratic fluid pressure computation
    void ComputePressureLEFT_RIGHT( Model<3>& );
    void NodesToFile( const Model<3>&, const char* region, SUBDOMAIN_PART );
    void FlagReadWriteTest( const Model<3>&, const char* variable );
    void ChangeBoundaryFlags( Model<3>&, const char* boundary, SUBDOMAIN_PART, const char* node_variable, VARIABLE_FLAG );
    void ChangeBoxBoundaryFlags( Model<3>&, BOX_BOUNDARY, const char* node_variable, VARIABLE_FLAG );
    double  LinearModelElementVolume( const Model<3U>& model, const std::string& flow_domain );
};


// non member functions
template<size_t dim>
double boxBoundaryFluxFEM_Quadratic( const Boundary<dim>& boundary, const csmp::Index& velo_key, double& boundary_area );

double  volumeAveragedVelocity( const Model<3>& model, const std::string& flow_domain, double& velo_integral );

void nodesToFile( const Model<3>& model, const char* region, SUBDOMAIN_PART part );

void nodeCoordinatesToText( const char* file, const std::set<Node<3U>*>& points );

void flagReadWriteTest( const Model<3>& model, const char* variable );


} // csmp

#endif // CSMP_STAGGERED_GRID_STOKES_SOLVER_EXAMPLE_H

