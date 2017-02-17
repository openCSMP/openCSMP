#ifndef CSMP_PARAMETER_H
#define CSMP_PARAMETER_H

#include "CSMP_number_types.h"
#include "Index.h"

namespace csmp {

/**
    @brief Record that describes the physical variables that are stored on a model:
    
    Stores variable name and notation (potential in correspondance with publications),
    SI unit, type, range (simulaton-specific physically meaningful minimum and maximum value, 
    and background information on the variable and its usage.
    
    Usage can be used to tell a GUI which
    parameters have to be input by the user, distinguishing them from 
    dependent (computed) values and/or the result of constitutive relation
    calculations. However, there is no fixed terminology communicating this
    yet.
    
    Parameter records are managed by the PropertyDatabase class.
    
    @author Stephan K. Matthai
    @date 1994
*/
struct Parameter {
    Parameter();
    Parameter( const Parameter& );
    ~Parameter();
    Parameter& operator=( const Parameter& );
  
    /// comparitor (to verify uniqueness of a new parameter definition
    bool    operator==( const csmp::Parameter& ) const;
    bool    operator!=( const csmp::Parameter& ) const;
  
    /// less-than operator so that Parameters can be stored in STL associative containers
    bool    operator<( const csmp::Parameter& ) const;

    /// define a physical variable at runtime, using console input (stdin)
    void           DefineFromStdin();
  
    /// checks whether the supplied value lies within min/max defined for this parameter
    bool           IsWithinRange( double64 ) const;
  
    /// expected range of this parameter in the specific simulation
    void           Range( double64& vmin, double64& vmax ) const;
  
    /// prints parameter record to screen
    void           Out() const;
  
    /// writes parameter record to a binary file
    bool           Out( FILE* fp ) const;
  
    /// reads parameter record from a binary file
    bool           In( FILE* fp );
    
    std::string    name;        ///< typically a human-readable name like 'fluid pressure' that can contain blanks (no tabs or line breaks)
    std::string    notation;    ///< e.g., k for permeability, v for velocity etc.
    std::string    unit;        ///< normally SI unit like kg/m2
    double64       min, max;    ///< physically meaningful value range, specific to simulation problem
    std::string    usage;       ///< with regard to computation: INPUT, COMPUTED etc.
    std::string    explanation; ///< how property is used, e.g., stress calculation etc.
    std::string    reference;   ///< to a paper that describes a related calculation
    mutable csmp::Index  key;   ///< variable accessor for this parameter as calculated and assigned by PropertyDatabase
 };

/// to print parameter description in compact form to an output stream
std::ostream&  operator<<( std::ostream&, const Parameter& );

} // csmp

#endif




