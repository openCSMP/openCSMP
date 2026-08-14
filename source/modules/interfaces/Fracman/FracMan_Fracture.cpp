#include "FracMan_Fracture.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {

// an ID which is icremented as a fracture is constructed or
// decremented when it is destructed to count all fractures
uint32_t FRACMAN_Fracture::global_id = 0;


FRACMAN_Fracture::FRACMAN_Fracture()
 :  id(++global_id), 
    fracture_set_id(1),
    aperture(0.001),
    compressibility(1.0e-7),
    permeability(1.0e-9)
 {
    // 'unit_normal' done with default constructor
 }
 

FRACMAN_Fracture::FRACMAN_Fracture( double        ap,
                              double        compr,
                              double        perm,
                              const list<double>& properties,
                              const list<Point<3>>& bry )
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

 
FRACMAN_Fracture::~FRACMAN_Fracture()
 {
    global_id--; 
 }
 
 
FRACMAN_Fracture::FRACMAN_Fracture( const FRACMAN_Fracture& ffr )
 :  id(++global_id) 
 {
    *this = ffr;
 }
 

uint32_t  FRACMAN_Fracture::ID() const
 {
    return id;
 }

string  FRACMAN_Fracture::TextID() const
 {
    char num[30];
    snprintf( num, sizeof(num), "%ul", id );
    
    return string(num);
 }
 
 
int FRACMAN_Fracture::SetID() const
 {
    return fracture_set_id;
 }
 
 
string  FRACMAN_Fracture::TextSetID() const
 {
    char num[30];
    snprintf( num, sizeof(num), "%ul", fracture_set_id );
    
    return string(num);
 }

 
FRACMAN_Fracture& FRACMAN_Fracture::operator=( const FRACMAN_Fracture& ffr )
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
bool  FRACMAN_Fracture::InitializeFrom( long nprops, ifstream& ifs )
 {
    char         text_line[256];
    char        *token;
    const char* const  delims = " ,=#%<>"; 
    long          n_pts, p_id;
    Point<3>  pt;
    
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
    for ( int i=0; i<n_pts; i++ ) {
         ifs.getline( text_line, 256 );
         // ID of point
         token = strtok( text_line, delims ); 
         if ( token == NULL ) return false;
         p_id  = atol( token ); 
         // coordinates of point
         token = strtok( NULL, delims ); // X
         if ( token == NULL ) return false;
         pt[0] = atof( token );
         token = strtok( NULL, delims ); // Y
         if ( token == NULL ) return false;
         pt[1] = atof( token );
         token = strtok( NULL, delims ); // Z
         if ( token == NULL ) return false;
         pt[2] = atof( token );
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
    pt[0] = atof( token );
    token = strtok( NULL, delims ); // Y
    if ( token == NULL ) return false;
    pt[1] = atof( token );
    token = strtok( NULL, delims ); // Z
    if ( token == NULL ) return false;
    pt[2] = atof( token );

    Point<3>  orig(0.,0.,0.);
    
    unit_normal = pt - orig;

    return true;
    
 } // end InitializeFrom

    
    
    
    
// compare perimeter length
bool FRACMAN_Fracture::operator<( const FRACMAN_Fracture& ffr ) const
 {
    return ( Perimeter() < ffr.Perimeter() );
 }
    
    
void FRACMAN_Fracture::BaryCenter( Point<3>& ctr ) const
 {
    list<Point<3>>::const_iterator  ita;
    double                             points;
  
    for ( ctr.Set(0.,0.,0.), points=0.0, 
          ita=boundary.begin(); ita!=boundary.end(); ita++, 
          points+=1.0 ) 
      ctr += (*ita);
      
    ctr /= points;
 }


void FRACMAN_Fracture::BaryCenter( double& x, double& y, double& z ) const
 {
    Point<3>  ctr(0.,0.,0.);
    
    BaryCenter( ctr );
    x = ctr[0];
    y = ctr[1];
    z = ctr[2];
 }
 
 
double  FRACMAN_Fracture::Perimeter() const
 {
    list<Point<3>>::const_iterator  ita, itb = boundary.begin();
    double                          perim(0.0);
  
    // getting segments from point to point
    for ( itb++, ita=boundary.begin(); itb!=boundary.end(); ita++, itb++ ) 
      {
         perim += distance( (*ita), (*itb) );
      }
      
    // adding segment from last point to first point
    perim += distance( (*boundary.rbegin()), (*boundary.begin()) );
    
    return perim; 
 }


// returns diameter of assuming that the fracture has a circular shape
double  FRACMAN_Fracture::Diameter() const
 {
    // u = 2 pi r
    return Perimeter() / PI; 
 }
 
 

void FRACMAN_Fracture::BoundingBox( Point<3>& cnr1, Point<3>& cnr8 ) const
 {
    cnr1 = (*min_element( boundary.begin(), boundary.end() ) );
    cnr8 = (*max_element( boundary.begin(), boundary.end() ) );
 }
 


void FRACMAN_Fracture::Erase()
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
void FRACMAN_Fracture::Move( double dx, double dy, double dz )
 {
    list<Point<3>>::iterator  ita;
    Point<3>                  displacement(dx,dy,dz);
    
    for ( ita=boundary.begin(); ita!=boundary.end(); ita++ ) 
      (*ita) += displacement;
 }

 
// move entire fracture by specified amount 
void FRACMAN_Fracture::Scale( double xfac, double yfac, double zfac )
 {
    list<Point<3>>::iterator  ita;
    Point<3>                  scale(xfac,yfac,zfac);
    
    for ( ita=boundary.begin(); ita!=boundary.end(); ita++ ) 
      (*ita) *= scale;
      
    // scaling the unit normal as well
    unit_normal *= scale;
    unit_normal.NormalizeLengthTo(1.);
 }


void FRACMAN_Fracture::Out() const
 {
    cout <<"\nFRACMAN_Fracture::Out: ID: "<< id;
    cout <<"\nAperture:               "<< aperture;
    cout <<"\nCompressibility:        "<< compressibility;
    cout <<"\nPermeability:           "<< permeability;
    
    if ( props.size() > 3 ) {
         list<double>::const_iterator  pit = props.begin();
         cout <<"\nOther property data:    "<< endl;
         pit++; pit++; pit++;
         
         while ( pit != props.end() ) {
              cout << (*pit) <<"  ";
              pit++;
           }
         cout << endl;
      }
    
    cout <<"\n\nPoints defining the perimeter of fracture: "<< endl;
    list<Point<3>>::const_iterator  ita;
    int                                a;
  
    for ( a=1, ita=boundary.begin(); ita!=boundary.end(); ita++, a++ ) 
      cout <<"\t"<< a <<": "<< (*ita)[0] <<"\t"<< (*ita)[1] <<"\t"<< (*ita)[2] << endl;
    
    cout <<"\nUnit normal to fracture plane:"<< endl;
    cout <<"\t"<< unit_normal[0] <<"\t"<< unit_normal[1];
    cout <<"\t"<< unit_normal[2] << endl;
    
    cout <<"\nPerimeter of fracture:  "<< Perimeter() << endl;
      
    cout << endl; 
 }

} // end namespace csmp
 
 
 
