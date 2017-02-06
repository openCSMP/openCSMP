#include <iostream>
#include "Exception.h"

using namespace std;

namespace csmp {

Exception::Exception()
 : csmp_exception_(ERROR),
   originator_("probably ErrorHandler"),
   message_("unspecified")
 {
 }


Exception::Exception( CSMP_MESSAGE err, 
                      const string& orig,
                      const string& msg )
 : csmp_exception_(err),
   originator_(orig),
   message_(msg)
 {
#ifndef NDEBUG 
    Out(cout);
    cout <<"\nHit return to continue."<< endl;
    getchar();
#endif
 }


Exception::Exception( CSMP_MESSAGE err, 
                      const string& orig,
                      const string& param,
                      const string& msg )
 : csmp_exception_(err),
   originator_(orig),
   message_(param)
 {
    message_ += "  ";
    message_ += msg;
#ifndef NDEBUG 
    Out(cout);
    cout <<"\nHit return to continue."<< endl;
    getchar();
#endif
 }


Exception::Exception( const Exception& excp )
 {
    *this = excp;
 }
 
 
Exception& Exception::operator=( const Exception& excp )
 {
    if ( &excp != this ) {
         csmp_exception_ = excp.csmp_exception_;
         originator_     = excp.originator_;
         message_        = excp.message_;
      }
    return *this;
 }        


string   Exception::What()  const
 {
    return message_;
 }


string   Exception::Originator()  const
 {
    return originator_;
 }
 
 
CSMP_MESSAGE  Exception::Message() const
 {
    return csmp_exception_;
 }


const char* Exception::what() const throw()
  {
    return message_.c_str();
  }


void Exception::Out(std::ostream& os) const
 {
    os <<"\n"<< string(parseMessage(csmp_exception_)) <<": "<< originator_ << endl;
    os << message_ << endl;
 }

Exception::~Exception() throw()
  {

  }


std::string  parseMessage( CSMP_MESSAGE msg )
 {
     if ( msg == INFO ) return string("INFO");
     if ( msg == EXCEPTION ) return string("EXCEPTION");
     if ( msg == WARNING ) return string("WARNING");
     if ( msg == ERROR ) return string("ERROR");
     if ( msg == FATAL_ERROR ) return string("FATAL_ERROR");
     
     return string("parseMessage(CSMP_MESSAGE): cannot parse message.");
 }



} // end namespace csmp 
