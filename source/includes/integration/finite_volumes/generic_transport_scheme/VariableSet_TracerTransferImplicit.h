//
//  VariableSet_TracerTransferImplicit
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#ifndef VARIABLE_SET_TRACER_TRANSFER_IMPLICIT_H
#define VARIABLE_SET_TRACER_TRANSFER_IMPLICIT_H

#include "VariableSet_TracerTransfer.h"

namespace csmp {

/**
    Variables used by the tracer-transfer scheme.
*/
struct VariableSet_TracerTransferImplicit : public VariableSet_TracerTransfer {
     template<size_t dim>
     explicit VariableSet_TracerTransferImplicit( const PropertyDatabase<dim>& );

   private:
     VariableSet_TracerTransferImplicit() = delete;
     VariableSet_TracerTransferImplicit( const VariableSet_TracerTransferImplicit& ) = delete;
     VariableSet_TracerTransferImplicit& operator=( const VariableSet_TracerTransferImplicit& ) = delete;

     /// verifies types and placements
     void CheckVariables() const;
};
 
} // end csmp

#endif /* defined(VARIABLE_SET_TRACER_TRANSFER_IMPLICIT_H) */
