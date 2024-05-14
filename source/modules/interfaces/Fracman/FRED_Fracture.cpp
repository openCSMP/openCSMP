#include "FRED_Fracture.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

// an ID which is icremented as a fracture is constructed or
// decremented when it is destructed to count all fractures
uint32_t FRED_Fracture::global_id = 0;


FRED_Fracture::FRED_Fracture()
 :  id(++global_id), 
    fracture_set_id(1),
    aperture(0.001),
    compressibility(1.0e-7),
    permeability(1.0e-9)
 {
    // 'unit_normal' done with default constructor
 }
 

FRED_Fracture::FRED_Fracture( double        ap,
                              double        compr,
                              double        perm,
                              const list<double>& properties,
                              const list<mjl::Point3D>& bry )
 :  id(++global_id),
    fracture_set_id(1),
    aperture(ap),
    compressibility(compr),
    permeability(perm),
    props(properties),
    boundary(bry)
 {
    // 'unit_normal' done with default constructor
 }

 
FRED_Fracture::~FRED_Fracture()
 {
    global_id--; 
 }
 
 
FRED_Fracture::FRED_Fracture( const FRED_Fracture& ffr )
 :  id(++global_id) 
 {
    *this = ffr;
 }
 

uint32_t  FRED_Fracture::ID() const
 {
    return id;
 }

string  FRED_Fracture::TextID() const
 {
    char num[30];
    snprintf( num, sizeof(num), "%ul", id );
    
    return string(num);
 }
 
 
int FRED_Fracture::SetID() const
 {
    return fracture_set_id;
 }
 
 
string  FRED_Fracture::TextSetID() const
 {
    char num[30];
    snprintf( num, sizeof(num), "%ul", fracture_set_id );
    
    return string(num);
 }

 
FRED_Fracture& FRED_Fracture::operator=( const FRED_Fracture& ffr )
 {
    if ( &ffr != this ) {
         id              = ffr.id;
         fracture_set_id = ffr.fracture_set_id;
         aperture        = ffr.aperture;
         compressibility = ffr.compressibility;
         permeability    = ffr.permeability;
         props           = ffr.props;
         boundary        = ffr.boundary;
         unit_normal     = ffr.unit_normal;
      }
    return *this;
 }
 
 
    
/**


Explanation and file entries:


@code
 ID#  #points #fset permeability      compressibility   aperture
 
    1    4    1     1.000000e+000     3.000000e-006     1.000000e-001
    1     3.333330e+001    -1.000000e+002     1.000000e+002
    2     3.333330e+001    -1.000000e+002    -1.000000e+002
    3     3.333330e+001     1.000000e+002    -1.000000e+002
    4     3.333330e+001     1.000000e+002     1.000000e+002
    0     1.000000e+000     0.000000e+000     0.000000e+000
@endcode
*/    
bool  FRED_Fracture::InitializeFrom( int nprops, ifstream& ifs )
 {
    char         text_line[256];
    char        *token;
    const char* const  delims = " ,=#%<>"; 
    uint32_t      n_pts, p_id;
    mjl::Point3D  pt;
    
    // 0. zapping all previous storage
    Erase();

    // ---------------------------------------------------
    // 1. reading headline
    // ---------------------------------------------------
    ifs.getline( text_line, 256 );
    // 1.1 ID number of fracture
    token           = strtok( text_line, delims ); 
    if ( token == NULL ) return false;
    id              = static_cast<uint32_t>(atoi( token )); 
    // 1.2 number of points defining it
    token           = strtok( NULL, delims );
    if ( token == NULL ) return false;
    n_pts           = atoi( token ); 
    // 1.3 ID number of fracture set to which fracture beuint32s 
    token           = strtok( NULL, delims );
    if ( token == NULL ) return false;
    fracture_set_id = atoi( token ); 
    // 1.4 permeability 
    token           = strtok( NULL, delims );
    if ( token == NULL ) return false;
    permeability    = atof( token );
    //     converting permeability from D to m2 (SI)
    permeability   -= 13.0;
    permeability    = pow( 10.0, permeability ); 
    props.push_back(permeability);
    // 1.5 compressibility 
    token           = strtok( NULL, delims );
    if ( token == NULL ) return false;
    compressibility = atof( token ); 
    //     converting compressibility from 1/kPa to 1/Pa (SI)
    compressibility *= 1000.0;
    props.push_back(compressibility);
    // 1.6 aperture 
    token           = strtok( NULL, delims );
    if ( token == NULL ) return false;
    aperture        = atof( token ); 
    props.push_back(aperture);
    
    // getting additional properties if such are specified
    if ( nprops > 3 ) {
         for ( int i=0; i<(nprops-3); i++ ) {
              token           = strtok( NULL, delims );
              if ( token == NULL ) return false;
              props.push_back( atof( token ) );
           }
      }

    if ( !boundary.empty() ) 
      boundary.erase( boundary.begin(), boundary.end() );

    // ---------------------------------------------------
    // 2. reading point lines
    // ---------------------------------------------------
    for ( uint32_t i=0u; i<n_pts; i++ ) {
         ifs.getline( text_line, 256 );
         // ID of point
         token = strtok( text_line, delims ); 
         if ( token == NULL ) return false;
         p_id  = atoi( token ); 
         // coordinates of point
         token = strtok( NULL, delims ); // X
         if ( token == NULL ) return false;
         pt(0) = atof( token ); 
         token = strtok( NULL, delims ); // Y
         if ( token == NULL ) return false;
         pt(1) = atof( token );         
         token = strtok( NULL, delims ); // Z
         if ( token == NULL ) return false;
         pt(2) = atof( token );
         // perimeter point list
         boundary.push_back( pt );   
      }
      
    // ---------------------------------------------------
    // 3. reading unit normal to fracture
    // ---------------------------------------------------
    ifs.getline( text_line, 256 );
    // ID of unit normal = 0
    token = strtok( text_line, delims ); 
    if ( token == NULL ) return false;
    // coordinate lengths of unit normal
    token = strtok( NULL, delims ); // X
    if ( token == NULL ) return false;
    pt(0) = atof( token ); 
    token = strtok( NULL, delims ); // Y
    if ( token == NULL ) return false;
    pt(1) = atof( token );         
    token = strtok( NULL, delims ); // Z
    if ( token == NULL ) return false;
    pt(2) = atof( token );

    mjl::Point3D  orig(0.,0.,0.);
    
    unit_normal.Set(orig,pt);

    return true;
    
 } // end InitializeFrom

    
    
    
    
// compare perimeter length
bool FRED_Fracture::operator<( const FRED_Fracture& ffr ) const
 {
    return ( Perimeter() < ffr.Perimeter() );
 }
    
    
void FRED_Fracture::BaryCenter( mjl::Point3D& ctr ) const
 {
    list<mjl::Point3D>::const_iterator  ita;
    double                             points;
  
    for ( ctr.Set(0.,0.,0.), points=0.0, 
          ita=boundary.begin(); ita!=boundary.end(); ita++, 
          points+=1.0 ) 
      ctr += (*ita);
      
    ctr /= points;
 }


void FRED_Fracture::BaryCenter( double& x, double& y, double& z ) const
 {
    mjl::Point3D  ctr(0.,0.,0.);
    
    BaryCenter( ctr );
    x = ctr.X();
    y = ctr.Y();
    z = ctr.Z();
 }
 
 
double  FRED_Fracture::Perimeter() const
 {
    list<mjl::Point3D>::const_iterator  ita, itb = boundary.begin();
    mjl::Edge3D                         edge;
    double                             perim(0.0);
  
    // getting segments from point to point
    for ( itb++, ita=boundary.begin(); itb!=boundary.end(); ita++, itb++ ) 
      {
         edge.Set( (*ita), (*itb) );
         perim += edge.Length();
      }
      
    // adding segment from last point to first point
    edge.Set( (*boundary.rbegin()), (*boundary.begin()) );
    perim += edge.Length();
    
    return perim; 
 }


// returns diameter of assuming that the fracture has a circular shape
double  FRED_Fracture::Diameter() const
 {
    // u = 2 pi r
    return Perimeter() / PI; 
 }
 
 

void FRED_Fracture::BoundingBox( mjl::Point3D& cnr1, mjl::Point3D& cnr8 ) const
 {
    cnr1 = (*min_element( boundary.begin(), boundary.end() ) );
    cnr8 = (*max_element( boundary.begin(), boundary.end() ) );
 }
 


void FRED_Fracture::Erase()
 {
    aperture        = 0.0;        
    compressibility = 0.0;
    permeability    = 0.0;
    id              = 0;
    fracture_set_id = 0;
    
    props.erase( props.begin(), props.end() );
    boundary.erase( boundary.begin(), boundary.end() );
    
    // unit_normal;
 }

// move entire fracture by specified amount 
void FRED_Fracture::Move( double dx, double dy, double dz )
 {
    list<mjl::Point3D>::iterator  ita;
    mjl::Point3D                  displacement(dx,dy,dz);
    
    for ( ita=boundary.begin(); ita!=boundary.end(); ita++ ) 
      (*ita) += displacement;
 }

 
// move entire fracture by specified amount 
void FRED_Fracture::Scale( double xfac, double yfac, double zfac )
 {
    list<mjl::Point3D>::iterator  ita;
    mjl::Point3D                  scale(xfac,yfac,zfac);
    
    for ( ita=boundary.begin(); ita!=boundary.end(); ita++ ) 
      (*ita) *= scale;
      
    // scaling the unit normal as well
    unit_normal.dest_ *= scale;
    unit_normal.Normalize();
 }


void FRED_Fracture::Out() const
 {
    cout <<"\nFRED_Fracture::Out: ID: "<< id;
    cout <<"\nAperture:               "<< aperture;
    cout <<"\nCompressibility:        "<< compressibility;
    cout <<"\nPermeability:           "<< permeability;
    
    if ( props.size() > 3 ) {
         list<double>::const_iterator  pit = props.begin();
         cout <<"\nOther property data:    "<< endl;
         pit = next(pit,3);
         
         while ( pit != props.end() ) {
              cout << (*pit) <<"  ";
              pit++;
           }
         cout << endl;
      }
    
    cout <<"\n\nPoints defining the perimeter of fracture: "<< endl;
    list<mjl::Point3D>::const_iterator  ita;
    int                                a;
  
    for ( a=1, ita=boundary.begin(); ita!=boundary.end(); ita++, a++ ) 
      cout <<"\t"<< a <<": "<< (*ita).X() <<"\t"<< (*ita).Y() <<"\t"<< (*ita).Z() << endl;
    
    cout <<"\nUnit normal to fracture plane:"<< endl;
    cout <<"\t"<< unit_normal.dest_.x_ <<"\t"<< unit_normal.dest_.y_;
    cout <<"\t"<< unit_normal.dest_.z_ << endl;
    
    cout <<"\nPerimeter of fracture:  "<< Perimeter() << endl;
      
    cout << endl; 
 }

} // end namespace csp 
 
 
 
