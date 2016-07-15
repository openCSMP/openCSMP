#ifndef SAMG_EXCEPTION_H
#define SAMG_EXCEPTION_H

#include <exception>
#include <stdexcept>
#include <string>

namespace csmp {

/// translation of SAMG error codes into human readable strings
class SAMG_Exception : public std::exception {
  public:
    explicit SAMG_Exception(int ierr) throw();
    SAMG_Exception(int ierr, std::string err_string) throw();
    virtual ~SAMG_Exception() throw();
    
    virtual const char* what() const throw();
    int ierr() const throw();
  private:
    int ierr_;
    std::string err_string_;
};

} // end namespace csmp

#endif
