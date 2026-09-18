// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SAMG_Exception.h"

namespace csmp {
  SAMG_Exception::SAMG_Exception(int ierr) throw() :
    ierr_(ierr)
    {
      switch(ierr/100) {
        case 0 : err_string_ = "Generic error"; break;
        case 1 : err_string_ = "Error during AMG setup phase"; break;
        case 2 : err_string_ = "Error during AMG solution phase"; break;
        case 8 : err_string_ = "Error in auxilary component"; break;
        case 9 : err_string_ = "Error during treatment of coarsest level"; break;
        default: err_string_ = "Error code that is not captured by this exception"; break;
      }
    }
    
  SAMG_Exception::SAMG_Exception(int ierr, std::string err_string) throw() :
    ierr_(ierr),
    err_string_(err_string)
    {}
    
  SAMG_Exception::~SAMG_Exception() throw() {}



  const char* SAMG_Exception::what() const throw() {
    return err_string_.c_str();
  }
  
  int SAMG_Exception::ierr() const throw() {
    return ierr_;
  }
    
  
} // end namespace csmp
