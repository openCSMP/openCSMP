#ifndef CSMP_STANDARD_IO_HANDLER_H
#define CSMP_STANDARD_IO_HANDLER_H

#include "CSMP_definitions.h"

namespace csmp {

/// Handles the command line and file IO for changeable input parameters during runs(may also be read from file)
class Standard_IO_Handler {
  public:
    explicit Standard_IO_Handler( const char* doc_file );
    Standard_IO_Handler();
    ~Standard_IO_Handler();
    
    bool        YesNo( const char question[150] );
    void        Protocol( const char* info );
    void        Protocol( const char* info, double64 parameter );
    double64    RecordChoice( const char* question );
    long        RecordIntChoice( const char* question );
    bool        RecordLogicalChoice( const char* question );
    std::string RecordLiteralChoice( const char* question );
    void        RecordInformation( const std::string& info );    
    void        Erase();
    void        Out() const;
    
  private:
    std::string             protocol_file_;
    std::list<std::string>  input_output;
};


} // end namespace csp


#endif
