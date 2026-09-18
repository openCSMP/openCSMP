// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Operation.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 22/7/20.
//

#include "Operation.h"

using namespace std;

namespace csmp {

/**
 enum OPERATION_TYPE { ADD, SUBTRACT, MULTIPLY, DIVIDE, MULTIPLY_WITH_TIME_INCREMENT };
*/
string  parseOperationType( OPERATION_TYPE type )
 { 
    switch( type ) {
         case ADD:
           return "ADD";
         case SUBTRACT:
           return "SUBTRACT";
         case MULTIPLY:
           return "MULTIPLY";
         case DIVIDE:
           return "DIVIDE";
      }
    return "ADD";
 }



    /// the operation type, execution level and name, are strung together into a key
Operation::Operation( const string& name, OPERATION_TYPE type, int execution_level, bool multiply_with_dt )
 :  name_(name), 
    key_(parseOperationType(type) + "_" + to_string(execution_level) + "_" + name),
    operation_type_(type),
    execution_level_(execution_level),
    multiply_with_time_increment_(multiply_with_dt)
 {
 } 
    

    
void Operation::Out() const
 {
     cout <<"\nOperation: "<< name_;
     cout <<"\n\tkey:             "<< key_;
     cout <<"\n\ttype:            "<< parseOperationType(operation_type_);
     cout <<"\n\texecution level: "<< execution_level_ << endl;
     if ( multiply_with_time_increment_ )
       cout <<"\n\toperator will be multiplied with time increment.\n";
 }


} // end csmp
