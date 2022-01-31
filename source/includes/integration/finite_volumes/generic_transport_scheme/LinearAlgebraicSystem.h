//
//  LinearAlgebraicSystem.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 6/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_LINEAR_ALGEBRAIC_SYSTEM_H
#define CSMP_LINEAR_ALGEBRAIC_SYSTEM_H

#include "SparseMatrix.h"

namespace csmp {

/**
   linear system LHS X = RHS
*/
struct LinearAlgebraicSystem {
    explicit LinearAlgebraicSystem( size_t m_x_n );
    void     Resize( size_t m_x_n );
    
    SparseMatrix           LHS;
    std::vector<double>  RHS;
    std::vector<double>  X;
};

} // end csmp

#endif /* CSMP_LINEAR_ALGEBRAIC_SYSTEM_H */
