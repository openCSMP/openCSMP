#ifndef MESCHACH_SOLVER_H
#define MESCHACH_SOLVER_H
#include "Solver.h"
#include "SparseMatrix.h"
//#include "meschach.h"

#ifdef CSMP_WITH_MESCHACH

namespace csmp {

class Meschach_Solver : public Solver {
public:
    Meschach_Solver();
    virtual ~Meschach_Solver();

protected:
    virtual void SolveWithMeschach(SparseMatrix& A,
                                   std::vector<double64>& b,
                                   std::vector<double64>& x,
                                   double64 solver_tolerance) = 0;

private:
    virtual void SolveMatrixEquation(SparseMatrix& A,
                                     std::vector<double64>& b,
                                     std::vector<double64>& x,
                                     size_t no_unknowns);

    double64 GuessResidual(const SparseMatrix& A,
                            const std::vector<double64>& b,
                            std::vector<double64>& x);
    
    double64 CalculateResidual(const SparseMatrix& A,
                                const std::vector<double64>& b,
                                const std::vector<double64>& x);

    double64 residual_factor_;
};
} // end namespace csmp

#endif  //CSMP_WITH_MESCHACH

#endif
