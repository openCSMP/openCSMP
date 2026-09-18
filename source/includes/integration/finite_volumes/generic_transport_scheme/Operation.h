// Copyright © 2020 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_OPERATION_H
#define CSMP_OPERATION_H

#include "CSMP_definitions.h"

namespace csmp {

/// time multiplication is handled inside of the MathOperator
enum OPERATION_TYPE { ADD, SUBTRACT, MULTIPLY, DIVIDE };
std::string parseOperationType( OPERATION_TYPE );

/**
      Specifies how to accumulate the Matrix- or Vector operators.
      TODO: add factor for the case where multiplication or division is the operation type?
*/
class Operation {
  public:
    /// the operation type, execution level and name, are strung together into a key
    Operation( const std::string& name, OPERATION_TYPE, int execution_level, bool multiply_with_dt ); 
    
    /// to determine place of operation in the sequence of operations
    bool operator<( const Operation& op ) const { return key_ < op.key_; }
    
    /// addition, subtraction, multiplication, division
    OPERATION_TYPE OperationType() const { return operation_type_; }
    
    /// when (in a sequence) the operation shall be performed
    int ExecutionLevel() const { return execution_level_; }
    
    void MultiplyWithTimeIncrement( bool multiply ) { multiply_with_time_increment_ = multiply; }
    bool MultiplyWithTimeIncrement() const { return multiply_with_time_increment_; }
    
    void Out() const;
    
  private:
    std::string    name_, key_;
    OPERATION_TYPE operation_type_  = ADD;
    int            execution_level_ = 1;
    bool           multiply_with_time_increment_ = false;
};

} // end csmp

#endif /* CSMP_OPERATION_H */
