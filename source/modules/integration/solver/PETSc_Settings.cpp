// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  PETSc_Settings.cpp
//  open-csmp-2024
//

#include "PETSc_Settings.h"

namespace csmp {

void PETSc_Settings::SetKSPType(const string& value)
{
    ksp_type_ = value;
}

void PETSc_Settings::SetPCType(const std::string& value)
{
    pc_type_ = value;
}

void PETSc_Settings::SetRelativeTolerance(double value)
{
    relative_tolerance_ = value;
}

void PETSc_Settings::SetAbsoluteTolerance(double value)
{
    absolute_tolerance_ = value;
}

void PETSc_Settings::SetMaximumIterations(int value)
{
    maximum_iterations_ = value;
}

void PETSc_Settings::SetMonitorTrueResidual(bool value)
{
    monitor_true_residual_ = value;
}

void PETSc_Settings::SetPrintConvergedReason(bool value)
{
    print_converged_reason_ = value;
}

void PETSc_Settings::SetUseInitialGuess(bool value)
{
    use_initial_guess_ = value;
}


const string& PETSc_Settings::GetKSPType() const
{
    return ksp_type_;
}

const string& PETSc_Settings::GetPCType() const
{
    return pc_type_;
}

double PETSc_Settings::GetRelativeTolerance() const
{
    return relative_tolerance_;
}

double PETSc_Settings::GetAbsoluteTolerance() const
{
    return absolute_tolerance_;
}

int PETSc_Settings::GetMaximumIterations() const
{
    return maximum_iterations_;
}

bool PETSc_Settings::GetMonitorTrueResidual() const
{
    return monitor_true_residual_;
}

bool PETSc_Settings::GetPrintConvergedReason() const
{
    return print_converged_reason_;
}

bool PETSc_Settings::GetUseInitialGuess() const
{
    return use_initial_guess_;
}

} // end namespace csmp
