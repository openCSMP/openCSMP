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
//    min(numeric_limits<double>::quiet_NaN()), - all hell breaks loose!
//    max(numeric_limits<double>::quiet_NaN()),
    min(-1.0e+25), // must be numbers that can be represented
    max(1.0e+25),
    usage("???"),
    explanation("???"),
    reference("???")
    // key
 {
 }

  
  
/*
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
*/




void  Parameter::Range( double& vmin, double& vmax ) const noexcept
 {
     vmin = min;
     vmax = max;
 }


ostream&  operator<<( ostream& stream, const Parameter& p )
 {
    stream << endl << setiosflags( ios::left );
    if ( p.name.size() > 20 )
      stream << setw(static_cast<int>(p.name.size()+2)) << p.name;
    else
      stream << setw(20) << p.name;
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


/// comparison operator for associative STL containers, relying on name comparison (alphabetical order)
bool Parameter::operator<( const csmp::Parameter& p ) const noexcept
 {
    // until 30/6/2026
    // return ( key < p.key );
    return name < p.name;
 }


// for testing whether a new parameter is unique
bool Parameter::operator==( const csmp::Parameter& p ) const noexcept
{
    return ( name==p.name && min==p.min && max==p.max &&
             notation==p.notation && unit==p.unit &&
             usage==p.usage &&
             explanation==p.explanation && reference==p.reference && key==p.key );
}


bool Parameter::operator!=( const csmp::Parameter& p ) const noexcept
{
    if ( name!=p.name )        return true;
    if ( key!=p.key )          return true;
    if ( min!=p.min )          return true;
    if ( max!=p.max )          return true;
    if ( notation!=p.notation ) return true;
    if ( unit!=p.unit )        return true;
    if ( usage!=p.usage )      return true;
    if ( explanation!=p.explanation ) return true;
    if ( reference!=p.reference )     return true;
    return false;
}



/// checks whether the supplied value is in the range stored in the parameter data
bool Parameter::IsWithinRange( double value ) const noexcept
 {
    if ( value > max || value < min ) return false;
    return true;
 }


void Parameter::Out() const 
 {
    cout <<"\nParameter::Out: "<< name <<"\n\n" << name; 
    cout <<"  " << notation;
    cout <<"  " << unit;
    cout <<"  " << min <<".."<< max;
    cout <<"\n\n" << key << endl;
    cout <<"  usage:       " << usage << endl;
    cout <<"  explanation: " << explanation << endl;
    cout <<"  reference:   " << reference << endl;
 }


bool Parameter::Out( std::fstream& fp ) const
  {
  binaryFileWrite( fp,  name.c_str() );
  binaryFileWrite( fp,  notation.c_str() );
  binaryFileWrite( fp,  unit.c_str() );
  key.Out(fp);
  fp.write( (char*) &min, sizeof(double) );
  fp.write( (char*) &max, sizeof(double) );
  binaryFileWrite( fp,  usage.c_str() );
  binaryFileWrite( fp,  explanation.c_str() );
  binaryFileWrite( fp,  reference.c_str() );

  return true; /// @todo (1-C) Meaningless return statement
  }


bool Parameter::In( fstream& fp )
  {
  char buf[399];
  binaryFileRead( fp,  buf );
  name = buf;
  binaryFileRead( fp,  buf );
  notation = buf;
  binaryFileRead( fp,  buf );
  unit = buf;
  key.In(fp);
  fp.read( (char*) &min, sizeof(double) );
  fp.read( (char*) &max, sizeof(double) );
  binaryFileRead( fp,  buf );
  usage = buf;
  binaryFileRead( fp,  buf );
  explanation = buf;
  binaryFileRead( fp,  buf );
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
 
 
 
