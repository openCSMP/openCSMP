#include "CSMP_mathUtilities.h"
#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

using namespace std;

namespace csmp {

/** rsqrt -- computes reciprocal square root */
double64
rsqrt( double64 val )
{
    return 1.0 / sqrt(val);
}


/** _v_norm1 -- computes (scaled) 1-norms of vectors */
double64 vector_norm1( vector<double64>& x, vector<double64>& scale )
{
	size_t	 i, dim = x.size();
	double64 s, sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += fabs(x[i]);
      }
      
    for ( i=0; i<dim; i++ ) {	
         s = scale[i];
		 sum += ( s == 0.0 ) ? fabs(x[i]) : fabs(x[i]/s);
	  }

	return sum;
}


double64 vector_norm2( vector<double64>& x, vector<double64>& scale )
{
	size_t	 i, dim = x.size();
	double64 s, sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += x[i] * x[i];
      }

    for ( i=0; i<dim; i++ ) {	
         s = scale[i];
		 sum += ( s == 0.0 ) ? x[i] * x[i] : (x[i]/s) * (x[i]/s);
	  }
	  
	return sqrt(sum);
}





/**
 
Converts (global) interpolation function derivative matrix into form where it
can be multiplied, for instance, with a 2DOF material property matrix or
for computation of 2D stiffness matrix.  

Original B: 
 
    | dN0dx dN1dx ... dN(nodes)dx |
B = | dN0dy dN1dy ... dN(nodes)dy |


Transformed B (see Zienkiewicz Vol.1 p. 48, eq. 3.10b): 
 
    | dN0dx   0   dN1dx  ... dN(nodes)dx   0   |
B = |   0   dN0dy   0    ...   0   dN(nodes)dy |
    | dN0dy dN0dx dN1dy  ... dN5dy dN(nodes)dx |
  
 
@section arguments Input Arguments 

The shape function derivative matrix DN for 1 degree of freedom per node. 

@return The expanded matrix is returned into the input matrix.

@section application Application

To pre-process DN matrices in multi-DOF computations.  

    @test OK SKM refactored and retested 2/12/2014
*/
void dN_To2DOF( size_t nodes, DenseMatrix<DM_MIN>& B )
 {
    assert( B.Cols() == nodes );
    const size_t dof(B.Cols() * 2U);

    // farming old matrix out into entries of new one
    // starting with last row
    B.Resize(3,B.Cols() * 2U);
    for ( size_t i=0; i<nodes; i++ )
      {
         B(2,i*2)   = B(1,i);
         B(2,i*2+1) = B(0,i);
      }
    
    // spreading out first row
    for ( size_t i=nodes; i>0; i-- ) B(0,(i-1)*2) = B(0,(i-1));
    // zeroing intermediate positions
    for ( size_t i=1; i<dof; i+=2 ) B(0,i) = 0.;
      
    // spreading out second row
    for ( size_t i=nodes; i>0; i-- ) B(1,(i-1)*2+1) = B(1,(i-1));
    // zeroing intermediate positions
    for ( size_t i=0; i<dof; i+=2 ) B(1,i) = 0.;
      
 } // end DNto2DOF




/**
    Transformation of interpolation function derivative matrix to 
    3 vectorial degrees of freedom per node as shown in Zienkewicz,
    vol. 1, p.133.
    
    @test OK SKM retested 2/12/2014
*/
void  dN_To3DOF( size_t nodes, DenseMatrix<DM_MIN>& B )
 {
    assert( B.Cols() == nodes );
    const size_t dof = B.Cols() * 3U;
    const size_t old_nodes = B.Cols();

    B.Resize(6,dof);
    
    // inserting values into the extra three bottom rows
    for ( size_t i=0; i<old_nodes; i++ )
      {
         // fourth row
         B(3,i*3)   = B(1,i); // d/dy
         B(3,i*3+1) = B(0,i); // d/dx
         B(3,i*3+2) = 0.0;    // 0
         
         // fifth row
         B(4,i*3)   = 0.0;    // 0
         B(4,i*3+1) = B(2,i); // d/dz
         B(4,i*3+2) = B(1,i); // d/dy

         // sixth row
         B(5,i*3)   = B(2,i); // d/dz
         B(5,i*3+1) = 0.0;    // 0
         B(5,i*3+2) = B(0,i); // d/dx
      }
      
    // spreading out the values in the first three rows
    // (going backward in order not to overwrite values in old
    //  position, while these are still needed)
    for ( size_t i=(old_nodes-1U); (i+1U)>0U; i-- )
      {
         // first row
         B(0,i*3)   = B(0,i);
         B(0,i*3+1) = 0.0;
         B(0,i*3+2) = 0.0; // zeros in the first row
         
         // second row
         B(1,i*3+1) = B(1,i);
         B(1,i*3)   = 0.0;
         B(1,i*3+2) = 0.0; // zeros in the first row
         
         // third row
         B(2,i*3+2) = B(2,i);
         B(2,i*3)   = 0.0;
         B(2,i*3+1) = 0.0; // zeros in the first row
      }
    
 } // end dN_To3DOF




/** _v_norm_inf -- computes (scaled) infinity-norm (supremum norm) of vectors */
double64 vector_norm_inf(  vector<double64>& x, vector<double64>& scale )
{
	size_t	 i, dim = x.size();
	double64 sum(0.0);

    if ( scale.empty() ) {
         scale.resize(dim);
         for ( i=0; i<dim; i++ ) sum += fabs(x[i]);
      }

cout <<"\nvector_norm_inf (in CSMP_math): function not implemented yet."<< endl;
return 0.;
/*
	double64 s, tmp, maxval = 0.0;
	if ( scale == (VEC *)NULL )
		for ( i = 0; i < dim; i++ )
		{	tmp = fabs(x->ve[i]);
			maxval = ms_max(maxval,tmp);
		}
	else if ( scale->dim < dim )
		error(E_SIZES,"_v_norm_inf");
	else
		for ( i = 0; i < dim; i++ )
		{	s = scale->ve[i];
			tmp = ( s== 0.0 ) ? fabs(x->ve[i]) : fabs(x->ve[i]/s);
			maxval = ms_max(maxval,tmp);
		}

	return maxval;
*/
}


void  vector_randomize( vector<double64>& x, double64 scale_fac )
 {
   // get a different seed every time 
    std::random_device rd;
    std::mt19937::result_type seed = rd() ^ (
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::system_clock::now().time_since_epoch()
                ).count() +
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> rndist(0,scale_fac);
 
    for ( vector<double64>::iterator it=x.begin();
          it!=x.end(); it++ )
      *it = rndist(gen);

 } // end vector_randomize




template<typename Var>
void out( const Var& obj )
 {
   cout <<"\ncontainer: "<< obj.size() << endl;
   
   int row_break(1), split_after(10); 
   for ( size_t i=0; i<obj.size(); i++, row_break++ )
     {
        if ( obj[i] > 0 ) cout <<" ";
        cout << obj[i] <<" ";
        if ( row_break == split_after )
          {
             cout << endl;
             row_break = 0;
          }
      }
   cout << endl;
 }




/// approximation of the complementary error function using a Chebyshev polynomial
double64  erfc_Chebyshev( double64 x )
{
	double64 z=fabs(x);
	double64 t=1.0/(1.0+0.5*z);
	double64 ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	
	return std::max( x >= 0.0 ? ans : 2.0-ans, 0. );
}


double64  erf_Chebyshev( double64 x )
{
	double64 z=fabs(x);
	double64 t=1.0/(1.0+0.5*z);
	double64 ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	
	// ascertain that the result is positive
	return std::max( 1. - (x >= 0. ? ans : 2. - ans), 0. );
}




// template explicit instantiations

template void out( const vector<ONE_BYTE_NUMBER>& );
template void out( const vector<size_t>& );
template void out( const vector<int32>& );
template void out( const vector<double64>& );



/// calculates arithmetic average of variables
template<typename Var> void average( const vector<Var>& var_vec, Var& var )
 {
    typename vector<Var>::const_iterator  it(var_vec.begin());
    var = *it++;
    while ( it != var_vec.end() ) var += *it++;
    var /= static_cast<double64>(var_vec.size());
 }
   
template void average( const vector<ScalarVariable >&, ScalarVariable& );

template void average( const vector<VectorVariable<1U> >&, VectorVariable<1U>& );
template void average( const vector<VectorVariable<2U> >&, VectorVariable<2U>& );
template void average( const vector<VectorVariable<3U> >&, VectorVariable<3U>& );

template void average( const vector<TensorVariable<1U> >&, TensorVariable<1U>& );
template void average( const vector<TensorVariable<2U> >&, TensorVariable<2U>& );
template void average( const vector<TensorVariable<3U> >&, TensorVariable<3U>& );





/// spline interpolation of values
double64 splineValue( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (1. - t)*y1 + t*y2 + t*(1.-t)*( a*(1.-t) + b*t);
}



double64 splineDerivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return (y2-y1)/( x2-x1 ) + (1.-2.*t)*( a*(1.-t)+b*t)/(x2-x1) + t*(1.-t)*(b-a)/(x2-x1);
}



double64 splineSecondDerivative( double64 x, double64 x1, double64 x2, double64 y1, double64 y2, double64 k1, double64 k2)
{
    const double64 a =  k1*( x2-x1 ) - ( y2 - y1 );
    const double64 b = -k2*( x2-x1 ) + ( y2 - y1 );
    const double64 t = ( x - x1) / ( x2 - x1 );

    return 2.*( b-2.*a +(a-b)*3.*t)/(x2-x1)/(x2-x1);
}




/** 
    Finding the root of the function by the Secant method

*/
double64 secant_method( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance)
{
  double64 xm, x0;
  double64 c;
  
  if (function(xmin) * function(xmax) < 0) {  // check the range for finding root
    do {
      
      x0 = (xmin * function(xmax) - xmax * function(xmin)) / (function(xmax) - function(xmin));
      c = function(xmin) * function(x0);

      xmin = xmax;
      xmax = x0;

      if (c == 0) break;
      
      xm = (xmin * function(xmax) - xmax * function(xmin)) / (function(xmax) - function(xmin));
      
    } while (fabs(xm - x0) >= tolerance); // repeat the loop
    
  } else
  {
    cout << "Error at CSMP_mathUtilities:the_secant_method Can not find a root in the given interval";
    xm =  std::numeric_limits<double>::quiet_NaN() ;
  }
  
  return xm;
}




/**
 
 finding maximum of a function using golden-section search, see:
 
 https://en.wikipedia.org/wiki/Golden-section_search
 
 */
double64 maximum_of_function( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance )
{
  
  double64 gr((sqrt(5.)+1.0)/2.0);  // golden ratio
  
  double64 xup(xmax - (xmax - xmin) / gr);
  double64 xlow(xmin + (xmax - xmin) / gr);
  
  do {
    if (function(xup) > function(xlow))   // check the function is maximized
    {
      xmax = xlow;
    } else
    {
      xmin = xup;
    }
    
    xup  = xmax - (xmax - xmin) / gr;
    xlow = xmin + (xmax - xmin) / gr;
    
  } while(fabs(xup - xlow) >= tolerance) ;
  
  return (xup + xlow)/2. ;
  
}
 
  
  
  
  
/**
 
 finding minimum of a function using golden-section search, see:
 
 https://en.wikipedia.org/wiki/Golden-section_search
 
 */
double64 minimum_of_function( double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 tolerance )
{
  
  double64 gr((sqrt(5.)+1.0)/2.0); // golden ratio
  
  double64 xup(xmax - (xmax - xmin) / gr);
  double64 xlow(xmin + (xmax - xmin) / gr);
  
  do {
    if (function(xup) < function(xlow))  // check the function is minimized
    {
      xmax = xlow;
    } else
    {
      xmin = xup;
    }
    
    xup  = xmax - (xmax - xmin) / gr;
    xlow = xmin + (xmax - xmin) / gr;
    
  } while(fabs(xup - xlow) >= tolerance) ;
  
  return (xup + xlow)/2. ;
  
}

  
  
  
/**
 
 finding the x value, in the range [xmin,xmax],  where function the df(x)/dx-(f(x)-f(x1))/(x-x1).. this line is the secant and tangent to f(x) at the shock point
 Mahyar: this is a basic function for two phase shock velocity analysis.

*/
double64 g( double64 x, double64 x1, double64 (*function)( double64 ), double64 (*dfunction)( double64 ))   // "g(x) = df(x)/dx-(f(x)-f(x1))/(x-x1)"
  {
    return dfunction(x)-(function(x)-function(x1))/(x-x1) ;
  }
  
  
  
double64 secant_line( double64 x1, double64 xmin, double64 xmax, double64 (*function)( double64 ), double64 (*dfunction)( double64 ), double64 tolerance)
  {
    
    double64 xm, x0;
    double64 c;
    
     if (g(xmin,x1,function ,dfunction) * g(xmax,x1,function ,dfunction) < 0) {  // check the range for finding root
       do {
        
        x0 = (xmin * g(xmax,x1,function ,dfunction) - xmax * g(xmin,x1,function ,dfunction)) / (g(xmax,x1,function ,dfunction) - g(xmin,x1,function ,dfunction));
        c  = g(xmin,x1,function ,dfunction) * g(x0,x1,function ,dfunction);
        
        xmin = xmax;
        xmax = x0;
        
        if (c == 0.0) break;
        
        xm = (xmin * g(xmax,x1,function ,dfunction) - xmax * g(xmin,x1,function ,dfunction)) / (g(xmax,x1,function ,dfunction) - g(xmin,x1,function ,dfunction));
        
       } while (fabs(xm - x0) >= tolerance); // repeat the loop
      
     }
    else
    {
      cout << "Error at CSMP_mathUtilities:the_secant_method Can not find a root in the given interval";
      xm =  std::numeric_limits<double>::quiet_NaN() ;
    }
    
    return xm;
  }





} // csmp









