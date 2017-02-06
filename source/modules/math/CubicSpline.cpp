#include "CubicSpline.h"

#include <cmath>

using namespace std;

namespace csmp {

// function prototypes
void spline( const std::vector<double64>& x, // x values    (0..n-1) 
             const std::vector<double64>& y, // f(x) values (0..n-1)
             double64 yp1, double64 ypn,    // slope at beginning and end
	           std::vector<double64>& y2 );    // f''(x) at above points


CubicSpline::CubicSpline()
 : x_range_(std::numeric_limits<double64>::quiet_NaN()),
   y_range_(std::numeric_limits<double64>::quiet_NaN())
 {
 }


CubicSpline::CubicSpline( const char* datafile )
 {
    Initialize( datafile );
 }


CubicSpline::CubicSpline( const CubicSpline& cs )
 : xa_(cs.xa_), ya_(cs.ya_), y2a_(cs.y2a_), 
   x_range_(cs.x_range_), y_range_(cs.y_range_)
 {
 }



CubicSpline::~CubicSpline()
 {
 }



CubicSpline& CubicSpline::operator=( const CubicSpline& cs )
 {
    if ( &cs != this ) {
         xa_      = cs.xa_;
         ya_      = cs.ya_;
         y2a_     = cs.y2a_;
         x_range_ = cs.x_range_; 
         y_range_ = cs.y_range_;
      }
    return *this;
 }
    

/** Reads two-column datafile with the format:

header line - just filename and explanation what's tabulated
blank line
number of x,y value pairs (the tabulated function)
derivate1 and derivative_n at the function origin and endpoint
x-y value pairs (one per line) 1..n 
*/
void CubicSpline::Initialize( const char* datafile )
 {
    string    str(datafile); str +=".txt";
    ifstream  ifs( str.c_str() );
    assert( ifs.is_open() );
    
    if ( !xa_.empty() ) xa_.erase( xa_.begin(), xa_.end() );
    if ( !ya_.empty() ) ya_.erase( ya_.begin(), ya_.end() );
    if ( !y2a_.empty() ) y2a_.erase( y2a_.begin(), y2a_.end() );
    
    // reading the title echoing it to screen
    char cstr[256U];
    ifs.getline( cstr, 256U );
    cout <<"\nCubicSpline::Initialize: reading: "<< cstr << endl;
    // blankline
    ifs.getline( cstr, 256U );
    // number of value pairs
    uint32 value_pairs, n(0U);
    ifs >> value_pairs;
    xa_.reserve( value_pairs );
    ya_.reserve( value_pairs );
    y2a_.resize( value_pairs );
    double64 x, y, y1, yn;
    // reading the function derivatives at the origin and the endpoint
    ifs >> y1 >> yn;    
    // reading the value pairs
    while ( n < value_pairs ) {
         ifs >> x >> y;
         xa_.push_back( x );
         ya_.push_back( y );
         n++;
      }
    ifs.close();
    

    // finding the ranges of x and f(x) values

    xa_min = *( min_element( xa_.begin(), xa_.end() ) );
    xa_max = *( max_element( xa_.begin(), xa_.end() ) );

    ya_min = *( min_element( ya_.begin(), ya_.end() ) );
    ya_max = *( max_element( ya_.begin(), ya_.end() ) );

    x_range_ = std::fabs(xa_max - xa_min);
    y_range_ = std::fabs(ya_max - ya_min);

    // creating the cubic spline from input points
    spline( xa_, ya_, y1, yn, y2a_ );    
 
 } // end Initialize

void  CubicSpline::Initialize( const std::vector<double64>& xa, const std::vector<double64>& ya,
                               const double64 y1, const double64 yn )
{

  if ( !xa_.empty() ) xa_.erase( xa_.begin(), xa_.end() );
  if ( !ya_.empty() ) ya_.erase( ya_.begin(), ya_.end() );
  if ( !y2a_.empty() ) y2a_.erase( y2a_.begin(), y2a_.end() );

   xa_ = xa;
   ya_ = ya;
   y2a_.resize( xa_.size() );

   // finding the ranges of x and f(x) values

   xa_min = *( min_element( xa_.begin(), xa_.end() ) );
   xa_max = *( max_element( xa_.begin(), xa_.end() ) );

   ya_min = *( min_element( ya_.begin(), ya_.end() ) );
   ya_max = *( max_element( ya_.begin(), ya_.end() ) );

   x_range_ = fabs(xa_max - xa_min);
   y_range_ = fabs(ya_max - ya_min);

   spline( xa_, ya_, y1, yn, y2a_ );

}


/// write internal data to screen
void CubicSpline::Out(std::ostream& os) const
 {
    os <<"\nCubicSpline: internal data: ";
    os <<"\n x, f(x) and f'(x) at "<< xa_.size() <<" user defined points.";
    vector<double64>::const_iterator  it1(ya_.begin());
    vector<double64>::const_iterator  it2(y2a_.begin());
    
    for ( vector<double64>::const_iterator
          it=xa_.begin(); it!=xa_.end(); it++, it1++, it2++ )
      os <<"\n"<< *it <<" "<< *it1 <<" "<< *it2;
    os << endl;
 }




/// creates cubic spline from input points
void spline( const std::vector<double64>& x, // x values    (0..n-1) 
             const std::vector<double64>& y, // f(x) values (0..n-1)
             double64 yp1, double64 ypn,    // slope at beginning and end
	           std::vector<double64>& y2 )     // f''(x) at above points
{
  assert( x.size() == y.size() );
  assert( x.size() >= 3U ); 
  y2.resize( x.size() );

	double64 p, qn, sig, un;

	const size_t n(y2.size());
	std::vector<double64>  u(n-1);
	
	if (yp1 > 0.99e30) y2[0]=u[0]=0.;
	else {
  		y2[0] = -0.5;
  		u[0]  = (3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  	}

	for ( uint32 i=1;i<n-1;i++) {
  		 sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  		 p=sig*y2[i-1]+2.;
  		 y2[i]=(sig-1.)/p;
  		 u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
  		 u[i]=(6.*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  	}

	if (ypn > 0.99e30) qn=un=0.;
	else {
  		 qn=0.5;
  		 un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  	}
  	
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.);
	
	// is this legitimate with an unsigned int?
	for ( int32 k=static_cast<int32>(n-2U); k>=0; k-- )
		y2[static_cast<uint32>(k)] = 
		  y2[static_cast<uint32>(k)] * y2[static_cast<uint32>(k+1)] + u[static_cast<uint32>(k)];
		
} // end spline





/// spline based interpolation, NumRecipes Chapter 3, p.118 modified
double64 splint( const vector<double64>& xa, 
                  const vector<double64>& ya, 
                  const vector<double64>& y2a, 
                  double64 x )
{
  assert( xa.size() >= 3U );

	long  klo(0U);
	long  khi(xa.size()-1U);
	
	while ( khi-klo > 1U ) {
  		long k=(khi+klo) >> 1;
  		if ( xa[k] > x ) khi=k;
  		else klo=k;
  	}
  	
	const double64 h(xa[khi] - xa[klo]);
	
	if ( fabs(h) < numeric_limits<double64>::epsilon() ) {
	     string err("\nsplint: Bad xa input to routine splint: ");
	     cout << err << x <<" ("<< h <<"), xa:"<< std::endl;
	     for ( uint32 i=0U; i<xa.size(); i++ ) cout <<" "<< xa[i];
	     cout << endl;
	     throw range_error( err.c_str() );
	  }
	  
	double64 a((xa[khi]-x) / h);
	double64 b((x-xa[klo]) / h);
	
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
	     	 +(b*b*b-b)*y2a[khi])*(h*h)/6.;
		
} // end splint

} // end csmp
