//
//  PETSc_Settings.h
//  open-csmp-2024
//

#ifndef CSMP_PETSC_SETTINGS_H
#define CSMP_PETSC_SETTINGS_H

#include "CSMP_definitions.h"
#include "SolverSettings.h"

using namespace std;

namespace csmp {

class PETSc_Settings : public SolverSettings {
public:
    PETSc_Settings() = default;
    ~PETSc_Settings() override = default;

    
    void SetKSPType( const string& value );
    void SetPCType(const string& value);
    void SetRelativeTolerance(double value);
    void SetAbsoluteTolerance(double value);
    void SetMaximumIterations(int value);
    void SetMonitorTrueResidual(bool value);
    void SetPrintConvergedReason(bool value);

    const string& GetKSPType() const;
    const string& GetPCType() const;
    double GetRelativeTolerance() const;
    double GetAbsoluteTolerance() const;
    int GetMaximumIterations() const;
    bool GetMonitorTrueResidual() const;
    bool GetPrintConvergedReason() const;

private:
    
    std::string ksp_type_ = "gmres";
    std::string pc_type_ = "ilu";

    double relative_tolerance_ = 1.0e-10;
    double absolute_tolerance_ = 1.0e-50;

    int maximum_iterations_ = 1000;

    bool monitor_true_residual_ = false;
    bool print_converged_reason_ = false;
    
};

} // namespace csmp

#endif




