#include <iomanip>
#include <algorithm>
#include "Parameter.h"
#include "binaryReadWrite.h"

using namespace std;

namespace csmp {

/**
    default constructor
    
    leaves Parameter 'unspecified'
    without notation or unit,
    range -1e25 to 1e25,
    no usage, explanation, reference etc.
*/
Parameter::Parameter()
 :  name("unspecified"),
    notation("-"),
    unit("-"),
//    min(numeric_limits<double64>::quiet_NaN()), - all hell breaks loose!
//    max(numeric_limits<double64>::quiet_NaN()),
    min(-1.0e+25), // must be numbers that can be represented
    max(1.0e+25),
    usage("???"),
    explanation("???"),
    reference("???")
    // key
 {
 }

  
Parameter::~Parameter()
  {
  }
  

Parameter::Parameter( const Parameter& p )
 :  name(p.name),
    notation(p.notation),
    unit(p.unit),
    min(p.min),
    max(p.max),
    usage(p.usage),
    explanation(p.explanation),
    reference(p.reference),
    key(p.key)
 {
 }
 

Parameter& Parameter::operator=( const Parameter& p )
 {
    if ( &p != this ) {
        name        = p.name;
        notation    = p.notation;
        unit        = p.unit;
        min         = p.min;
        max         = p.max;
        usage       = p.usage;
        explanation = p.explanation;
        reference   = p.reference;
        key         = p.key;
      }
    return *this;
 }


bool Parameter::operator==( const csmp::Parameter& p ) const
  {
    return ( name==p.name && min==p.min && max==p.max &&
             notation==p.notation && unit==p.notation && usage==p.usage &&
             explanation==p.explanation && reference==p.reference && key==p.key );
  }


bool Parameter::operator!=( const csmp::Parameter& p ) const
  {
     if ( name!=p.name ) return false;
     if ( key!=p.key ) return false;
     if ( min!=p.min ) return false;
     if ( max!=p.max ) return false;
     if ( name!=p.name ) return false;
     if ( notation!=p.notation )  return false;
     if ( unit!=p.unit ) return false;
     if ( usage!=p.usage ) return false;
     if ( explanation!=p.explanation ) return false;
     if ( reference!=p.reference ) return false;
     return true;
  }


/// comparison operator for associative STL containers, relying on a key comparison
bool Parameter::operator<( const csmp::Parameter& p ) const
 {
    return ( key < p.key );
 }


void  Parameter::Range( double64& vmin, double64& vmax ) const
 {
     vmin = min;
     vmax = max;
 }


ostream&  operator<<( ostream& stream, const Parameter& p )
 {
    stream << endl << setiosflags( ios::left ); 
    stream << setw(40) << p.name;
    stream << setw(6) << p.notation;
    string  str("[");
    str += p.unit;
    str += "]"; 
    stream <<" range: ";
    stream.setf(ios::scientific);
    stream << setiosflags( ios::right );
    const long prec = cout.precision(2);
    stream << setw(9) << p.min <<" - "<< setw(9) << p.max <<", "<< p.key; 
    stream.unsetf( ios::scientific );
    stream.precision(prec);
    stream << resetiosflags( ios::adjustfield );
    return stream;
 }

/// checks whether the supplied value is in the range stored in the parameter data
bool Parameter::IsWithinRange( double64 value ) const
 {
    if ( value > max || value < min ) return false;
    return true;
 }


void Parameter::Out(std::ostream& os) const 
 {
    os <<"\nParameter::Out: "<< name <<"\n\n" << name; 
    os <<"  " << notation;
    os <<"  " << unit;
    os <<"  " << min <<".."<< max;
    os <<"\n\n" << key << endl;
    os <<"  usage:       " << usage << endl;
    os <<"  explanation: " << explanation << endl;
    os <<"  reference:   " << reference << endl;
 }


bool Parameter::Out( FILE* fp ) const
  {
  skm_C_fwrite( fp,  name.c_str() );
  skm_C_fwrite( fp,  notation.c_str() );
  skm_C_fwrite( fp,  unit.c_str() );
  key.Out(fp);
  fwrite( (void*) &min, sizeof(double64), 1, fp );
  fwrite( (void*) &max, sizeof(double64), 1, fp );
  skm_C_fwrite( fp,  usage.c_str() );
  skm_C_fwrite( fp,  explanation.c_str() );
  skm_C_fwrite( fp,  reference.c_str() );

  return true; /// @todo (1-C) Meaningless return statement
  }


bool Parameter::In( FILE* fp )
  {
  char buf[399];
  skm_C_fread( fp,  buf );
  name = buf;
  skm_C_fread( fp,  buf );
  notation = buf;
  skm_C_fread( fp,  buf );
  unit = buf;
  key.In(fp);
  fread( (void*) &min, sizeof(double64), 1, fp );
  fread( (void*) &max, sizeof(double64), 1, fp );
  skm_C_fread( fp,  buf );
  usage = buf;
  skm_C_fread( fp,  buf );
  explanation = buf;
  skm_C_fread( fp,  buf );
  reference = buf;

  return true; /// @todo (1-C) Meaningless return statement
  }
 
 
void Parameter::DefineFromStdin()
 {
    string in_text;
    string vtype;
    cout <<"\n\nParameter::DefineFromStdin: PLEASE ENTER NEW PROPERTY DATA\n";
    cout <<"\nvariable name: ";
    cin >> name;
    cout <<"\nnotation: ";
    cin >> notation;
    cout <<"\nunit: ";
    cin >> unit;
    cout <<"\nmininum plausible value: ";
    cin >> min;
    cout <<"\nmaximum plausible value: ";
    cin >> max;
    cout <<"\nplacement: ";
    cin >> in_text; 
    key.place = parsePlacement(in_text.c_str());
    /// Roman, 2013: Added explicit type identifications
    cout <<"\nType of VARIABLE: SCALAR, VECTOR, TENSOR, ARRAY, FLAGGEDARRAY";
    cin  >>vtype;
    std::transform(vtype.begin(),vtype.end(),vtype.begin(),::toupper);
    key.type = parseType( vtype.c_str());
    cout <<"\nData depth (1-Scalar, dim-Vector, dim*dim-Tensor, size- Array or FlaggedArray";
    cin >> key.dataDepth;
    cout <<"\nusage: ";
    cin >> usage;
    cout <<"\nexplanation: ";
    cin >> explanation;
    cout <<"\nreference: ";
    cin >> reference;
 }    
 
} // end namespace csmp
 
 
 
