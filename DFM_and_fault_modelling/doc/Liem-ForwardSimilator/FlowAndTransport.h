// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
// Created by Michael Liem on 08.06.22.
//

#ifndef CSMP_API_FLOWANDTRANSPORT_H
#define CSMP_API_FLOWANDTRANSPORT_H

#include "CSMP_definitions.h"

// the CSMP model
#include "Model.h"
#include "Region.h"
#include "Boundary.h"

// input interfaces
#include "TextFileIO.h"
#include "TextInterface.h"

// utility functions
#include "ConstantFactor.h"

// SAMG
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#include "EigenSolver.h"

// the FE algorithm
#include "PDE_Integrator.h"
// PDE operators building the FE algorithm
#include "NumIntegral_NT_op_N_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_lhsop_N_dV.h"
#include "VelocityAndVolumeFlux.h"

// Transport scheme
//#include "NodeCenteredFiniteVolumeTransport.h"
//#include "DESAdvectionDiffusion.h"

#include <filesystem>

namespace csmp{

  class FlowAndTransport {

  public:
    static void LoadInputData( Model<2U>& model, const std::string& input_path_and_name );
    static void SolveForPressure( Model<2U>& model, const PropertyDatabase<2>& p_ref );

    static void SolveTransientPressure( Model<2U>& model, PDE_Integrator<2U,Element>& transient_pressure);

  };

}// namespace csmp


#endif //CSMP_API_FLOWANDTRANSPORT_H
