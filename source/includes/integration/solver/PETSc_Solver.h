//
//  PETSc_Solver.h
//  open-csmp-2024
//
//  Created by Mariia Dubnytska on 04.05.2026.
//

#ifndef CSMP_PETSC_SOLVER_H
#define CSMP_PETSC_SOLVER_H

#include "Solver.h"
#include "PETSc_Settings.h"

#include <petsc.h>

#include <array>
#include <charconv>
#include <optional>
#include <regex>
#include <string_view>


/**
  @brief a csmp-petsc interface solver class
 */
namespace csmp {

class SparseMatrix;
class CompressedRowMatrix;

class PETSc_Solver : public Solver
  {
  public:
    PETSc_Solver();
    explicit PETSc_Solver(PETSc_Settings* settings);
    ~PETSc_Solver() override;

    void InputSolverSettings(SolverSettings& settings) override final;
    PETSc_Settings* GetSolverSettings() override final;

    std::string Name() const override { return "PETSc_Solver"; }

  protected:

   /**
     @brief solves matrix equation A x = b using PETSc Matrix solver.
     
     @param A csmp sparse matrix
     @param b RHS vector
     @param x vector of unknowns
     @param no_unknowns degrees of freedom of the problem
   */
#if defined(CSMP_WITH_PETSC_SOLVER)
    virtual void SolveMatrixEquation(SparseMatrix& A,
                                     std::vector<double>& b,
                                     std::vector<double>& x,
                                     size_t no_unknowns) override;

    virtual void SolveMatrixEquation(CompressedRowMatrix& A,
                                     std::vector<double>& b,
                                     std::vector<double>& x,
                                     size_t no_unknowns) override;
#endif

  private:

#if defined(CSMP_WITH_PETSC_SOLVER)
    void ConfigureKSP(KSP ksp);
    void CheckPETSc_Initialized() const;
#endif
    PETSc_Settings* settings_ = nullptr;
    bool owns_settings_ = false;
};



// ERROR CODES

struct ErrorInfo {
    int code;
    std::string_view symbol;
    std::string_view description;
};

inline constexpr std::array kKnownErrors{
    ErrorInfo{0,   "PETSC_SUCCESS",              "Success"},
    ErrorInfo{55,  "PETSC_ERR_MEM",              "Out of memory"},
    ErrorInfo{56,  "PETSC_ERR_SUP",              "Operation unsupported for this object type"},
    ErrorInfo{57,  "PETSC_ERR_SUP_SYS",          "Operation unsupported on this system"},
    ErrorInfo{58,  "PETSC_ERR_ORDER",            "Operation performed in the wrong order"},
    ErrorInfo{59,  "PETSC_ERR_SIG",              "Signal received"},
    ErrorInfo{60,  "PETSC_ERR_ARG_SIZ",          "Nonconforming object sizes"},
    ErrorInfo{61,  "PETSC_ERR_ARG_IDN",          "Argument aliasing is not permitted"},
    ErrorInfo{62,  "PETSC_ERR_ARG_WRONG",        "Invalid argument"},
    ErrorInfo{63,  "PETSC_ERR_ARG_OUTOFRANGE",   "Argument out of range"},
    ErrorInfo{64,  "PETSC_ERR_ARG_CORRUPT",      "Null or corrupted PETSc object argument"},
    ErrorInfo{65,  "PETSC_ERR_FILE_OPEN",        "Unable to open file"},
    ErrorInfo{66,  "PETSC_ERR_FILE_READ",        "Unable to read file"},
    ErrorInfo{67,  "PETSC_ERR_FILE_WRITE",       "Unable to write file"},
    ErrorInfo{68,  "PETSC_ERR_ARG_BADPTR",       "Invalid pointer argument"},
    ErrorInfo{69,  "PETSC_ERR_ARG_NOTSAMETYPE",  "Arguments must have the same object type"},
    ErrorInfo{70,  "PETSC_ERR_POINTER",          "Invalid memory address"},
    ErrorInfo{71,  "PETSC_ERR_MAT_LU_ZRPVT",     "Zero pivot during LU factorization"},
    ErrorInfo{72,  "PETSC_ERR_FP",               "Floating-point exception"},
    ErrorInfo{73,  "PETSC_ERR_ARG_WRONGSTATE",   "Object is in the wrong state"},
    ErrorInfo{74,  "PETSC_ERR_COR",              "Corrupted PETSc object"},
    ErrorInfo{75,  "PETSC_ERR_ARG_INCOMP",       "Incompatible arguments"},
    ErrorInfo{76,  "PETSC_ERR_LIB",              "Error in a library called by PETSc"},
    ErrorInfo{77,  "PETSC_ERR_PLIB",             "PETSc internal inconsistency"},
    ErrorInfo{78,  "PETSC_ERR_MEMC",             "Memory corruption"},
    ErrorInfo{79,  "PETSC_ERR_FILE_UNEXPECTED",  "Unexpected data in file"},
    ErrorInfo{80,  "PETSC_ERR_ARG_NOTSAMECOMM",  "Arguments use different communicators"},
    ErrorInfo{81,  "PETSC_ERR_MAT_CH_ZRPVT",     "Zero pivot during Cholesky factorization"},
    ErrorInfo{82,  "PETSC_ERR_CONV_FAILED",      "Iterative method failed"},
    ErrorInfo{83,  "PETSC_ERR_USER",             "Required user function was not provided"},
    ErrorInfo{84,  "PETSC_ERR_INT_OVERFLOW",     "Integer overflow"},
    ErrorInfo{85,  "PETSC_ERR_ARG_NULL",         "Required argument is null"},
    ErrorInfo{86,  "PETSC_ERR_ARG_UNKNOWN_TYPE", "Unknown registered type"},
    ErrorInfo{87,  "PETSC_ERR_MPI_LIB_INCOMP",   "Incompatible MPI library at runtime"},
    ErrorInfo{88,  "PETSC_ERR_SYS",              "System-call error"},
    ErrorInfo{89,  "PETSC_ERR_ARG_TYPENOTSET",   "Object type has not been set"},
    ErrorInfo{91,  "PETSC_ERR_NOT_CONVERGED",    "Solver did not converge"},
    ErrorInfo{95,  "PETSC_ERR_USER_INPUT",       "Missing or incorrect user input"},
    ErrorInfo{96,  "PETSC_ERR_GPU_RESOURCE",     "Unable to load a GPU resource"},
    ErrorInfo{97,  "PETSC_ERR_GPU",              "GPU API error"},
    ErrorInfo{98,  "PETSC_ERR_MPI",              "MPI error"},
    ErrorInfo{100, "PETSC_ERR_MEM_LEAK",         "Memory allocation/free imbalance"},
    ErrorInfo{101, "PETSC_ERR_PYTHON",           "Python exception"},
};

inline constexpr ErrorInfo kUnknown{
    0, "PETSC_ERR_UNKNOWN", "Unknown, non-standard, or version-specific PETSc error code"
};

[[nodiscard]] constexpr ErrorInfo describe(int code) noexcept {
    for (const auto& error : kKnownErrors)
        if (error.code == code) return error;

    auto result = kUnknown;
    result.code = code;
    return result;
}


/** Usage
 
    petsc::printError(63);
    // PETSC_ERR_ARG_OUTOFRANGE (63) — Argument out of range

    petsc::printError(-3);
    // PETSC_ERR_UNKNOWN (-3) — Unknown, non-standard, or version-specific PETSc error code
 */
inline void printError(int code)
{
    const auto error = describe(code);

    std::println(
        "{} ({}) — {}",
        error.symbol,
        error.code,
        error.description
    );
}

// Extracts codes from lines such as:
//   "PETSc error code 63"
//   "error code: 82"
//   "PETSC ERROR: error code = 95"
[[nodiscard]] inline std::optional<ErrorInfo>
parse_error(std::string_view text)
{
    static const std::regex pattern{
        R"(\b(?:PETSc\s+)?error\s+code\s*[:=]?\s*(-?\d+)\b)",
        std::regex::icase
    };

    std::match_results<std::string_view::const_iterator> match;
    if (!std::regex_search(text.begin(), text.end(), match, pattern))
        return std::nullopt;

    int code{};
    const auto number = std::string_view{match[1].first,
                                         static_cast<std::size_t>(match[1].length())};
    const auto [end, ec] =
        std::from_chars(number.data(), number.data() + number.size(), code);

    if (ec != std::errc{} || end != number.data() + number.size())
        return std::nullopt;

    return describe(code);
}

} // end namespace csmp

#endif
