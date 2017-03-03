#ifndef GAUSSIAN_ELIMINATION_SOLVER_H
#define GAUSSIAN_ELIMINATION_SOLVER_H
#include "Meschach_Solver.h"
#include "SparseMatrix.h"
#ifdef CSMP_WITH_MESCHACH

namespace csmp {

class Gauss_Solver : public Meschach_Solver {

public:
    Gauss_Solver();
    virtual ~Gauss_Solver();
protected:
    virtual void SolveWithMeschach(csmp::SparseMatrix& A,
                                   std::vector<double64>& b,
                                   std::vector<double64>& x,
                                   double64 solver_tolerance);

private:
    double64 pivot_factor;
};

} // end namespace csmp

#endif //CSMP_WITH_MESCHACH

#endif
