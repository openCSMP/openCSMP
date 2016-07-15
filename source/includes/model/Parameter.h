#ifndef CSMP_PARAMETER_H
#define CSMP_PARAMETER_H

#include "CSMP_number_types.h"
#include "Index.h"

namespace csmp {

/**
    Public record with information about physical variables:
    
    Stores notation (in correspondance with publications),
    SI unit, type, range, name, and background information on the variable
    and its usage. 
    
    Usage has already been applied to tell a GUI which
    Parameters are input, dependent and/or the result of constitutive relation
    calculations. However, there is no fixed terminology communicating this
    yet.
    
    Is managed by the PropertyDatabase class.
    
    @author SKM (1994)
*/
struct Parameter {
    Parameter();
    Parameter( const Parameter& );
    ~Parameter();
    Parameter& operator=( const Parameter& );

    bool    operator==( const csmp::Parameter& ) const;
    bool    operator!=( const csmp::Parameter& ) const;
    bool    operator<( const csmp::Parameter& ) const;

    /// define a physical variable at runtime, using console input
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
    
    std::string    name;
    std::string    notation;    ///< e.g., k for permeability, v for velocity etc.
    std::string    unit;        ///< normally SI unit like kg/m2
    double64       min, max;    ///< physically meaningful value range, specific to simulation problem
    std::string    usage;       ///< with regard to computation: INPUT, COMPUTED etc.
    std::string    explanation; ///< how property is used, e.g., stress calculation etc.
    std::string    reference;   ///< to a paper that describes a related calculation
    mutable csmp::Index  key;   ///< variable accessor for this parameter as calculated and assigned by PropertyDatabase
 };

std::ostream&  operator<<( std::ostream&, const Parameter& );

inline bool Parameter::IsWithinRange( double64 value ) const
 {
    if ( value > max || value < min ) 
      return false;
    return true;
 }

} // csmp

#endif




