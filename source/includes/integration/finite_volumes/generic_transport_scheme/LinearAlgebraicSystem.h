// Copyright © 2020 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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
