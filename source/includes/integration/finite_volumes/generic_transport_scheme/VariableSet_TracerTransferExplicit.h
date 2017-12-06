//
//  VariableSet_TracerTransferExplicit
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#ifndef VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H
#define VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H

#include "VariableSet_TracerTransfer.h"

namespace csmp {

/**
    Variables used by the tracer-transfer scheme.
*/
struct VariableSet_TracerTransferExplicit : public VariableSet_TracerTransfer {
     template<size_t dim>
     explicit VariableSet_TracerTransferExplicit( const PropertyDatabase<dim>& );

   private:
     VariableSet_TracerTransferExplicit() = delete;
     VariableSet_TracerTransferExplicit( const VariableSet_TracerTransferExplicit& ) = delete;
     VariableSet_TracerTransferExplicit& operator=( const VariableSet_TracerTransferExplicit& ) = delete;

     /// verifies types and placements
     void CheckVariables() const;
};
 
} // end csmp

#endif /* defined(VARIABLE_SET_TRACER_TRANSFER_EXPLICIT_H) */
