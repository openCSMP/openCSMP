//
//  PDE_Integrator_new_functions.cpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 18/5/2026.
//  Copyright © 2026 Stephan Matthai. All rights reserved.
//

#include "PDE_Integrator.h"

using namespace std;

namespace csmp {

// Inside your PDE_Integrator element assembly loop
// Assuming:
// - 'element' is the current finite element.
// - 'prop_key' is the CSMP::Index for the variable being assembled.
// - 'var_size' is the number of components (e.g., 2 or 3 for VECTOR).
// - 'offset' is the starting index in the global block for this variable.
// - 'local_K' is your tensor element stiffness matrix.

for (uint32_t local_n1 = 0; local_n1 < element->Nodes(); ++local_n1) {
    auto* node1 = element->N(local_n1);
    const size_t global_idx1 = node1->Idx();
    const size_t reduced_row_pos = DOF_indexes[global_idx1];

    // LOOP 1: The Equation Row (Component 'i' of Node 1)
    for (uint32_t dim_i = 0; dim_i < var_size; ++dim_i) {
        
        // Determine the correct flag index for the ROW component
        const int flag_i = (prop_key.type == SCALAR || prop_key.type == ARRAY) ? 0 : static_cast<int>(dim_i);
        
        // If THIS specific component is Dirichlet-flagged, we do not assemble an equation row for it.
        if (reduced_row_pos == NULL_IDX || node1->Status(prop_key, flag_i) == DIRICH) continue;

        // Calculate the global reduced row index for G.AddIf()
        const size_t target_row = reduced_row_pos * var_size + dim_i + offset;
        
        // Calculate the global unreduced row index for pivotVector_
        const size_t unreduced_row = global_idx1 * var_size + dim_i + offset;

        for (uint32_t local_n2 = 0; local_n2 < element->Nodes(); ++local_n2) {
            auto* node2 = element->N(local_n2);
            const size_t global_idx2 = node2->Idx();
            const size_t reduced_col_pos = DOF_indexes[global_idx2];

            // LOOP 2: The Variable Column (Component 'j' of Node 2)
            for (uint32_t dim_j = 0; dim_j < var_size; ++dim_j) {
                
                // Get the tensor stiffness coupling between (Node1, dim_i) and (Node2, dim_j)
                double K_val = local_K(local_n1, dim_i, local_n2, dim_j);
                
                // Skip zero-entries early to save branch/memory overhead
                if (K_val == 0.0) continue; 

                // Determine the correct flag index for the COLUMN component
                const int flag_j = (prop_key.type == SCALAR || prop_key.type == ARRAY) ? 0 : static_cast<int>(dim_j);

                // Check if the COLUMN component is Dirichlet-flagged
                if (reduced_col_pos == NULL_IDX || node2->Status(prop_key, flag_j) == DIRICH) {
                    
                    // COLUMN IS DIRICHLET: Perform Exact Algebraic Elimination
                    // Fetch the known boundary value for this specific component
                    // (Replace GetDirichletValue with your actual CSMP boundary value fetcher)
                    double u_b = node2->GetDirichletValue(prop_key, flag_j); 
                    
                    // Subtract the cross-coupling force from the unreduced pivot vector
                    pivotVector_[unreduced_row] -= K_val * u_b;
                } 
                else {
                    // BOTH COMPONENTS ACTIVE: Assemble safely into the reduced global matrix
                    const size_t target_col = reduced_col_pos * var_size + dim_j + offset;
                    
                    // AddIf guarantees we only write to slots allocated by the SparsityPattern
                    G.AddIf(target_row, target_col, K_val);
                }
            }
        }
    }
}

} // end csmp
