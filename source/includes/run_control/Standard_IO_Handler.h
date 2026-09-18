// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_STANDARD_IO_HANDLER_H
#define CSMP_STANDARD_IO_HANDLER_H

#include "CSMP_definitions.h"

namespace csmp {

/** @brief Primitive handler for command line and file input/output 
    of user-defined input parameters during runs (input may also be read from file)
    
    @author Stephan K. Matthai
    @date 1996
*/
class Standard_IO_Handler {
  public:
    /// constructs handler, specifying which file the collected information shall be written to
    explicit Standard_IO_Handler( const char* doc_file );
    Standard_IO_Handler();
    ~Standard_IO_Handler();
  
    /// prints question onto screen and collects input from user (stdin)
    bool        YesNo( const std::string& question, const char* help_message = nullptr );
  
    /// prompts user to enter a floating point value which is collected from stdin returned
    double    RecordChoice( const char* question );
  
    /// prompts user to enter an integer value which is collected from stdin returned
    long        RecordIntChoice( const char* question );
  
    /// prompts user to make a binary choice which is collected from stdin returned
    bool        RecordLogicalChoice( const char* question );
  
    /// prompts user to input a word of text in response to the posed question
    std::string RecordLiteralChoice( const char* question );
  
    /// prompts user to input a string of text in response to the posed question
    std::string RecordMultipleLiteralChoice( const char* question );
  
    /// collects text from stdin and stores it in internal container of strings
    void        RecordInformation( const std::string& info );
  
    /// removes all logged information
    void        Erase();

    /// appends info string to doc_file (specified during handler construction)
    void        Protocol( const char* info );
  
    /// appends info string + numeric value to doc_file
    void        Protocol( const char* info, double parameter );
  
    /// prints all the recorded information to the output text file 'doc_file'
    void        Out() const;
    
    /// sets default help message fot this handler
    void        SetDefaultHelp( const std::string& default_help );

  private:
    std::string             protocol_file_;  ///< the ASCII text file to which the log will be written
    std::list<std::string>  input_output;    ///< the collected information
    std::string             default_help_{ "No available help..." }; /// default help message if none is provided
};


} // end namespace csp


#endif
