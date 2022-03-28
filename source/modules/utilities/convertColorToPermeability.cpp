#include "convertColorToPermeability.h"
#include "Matrix.h"

using namespace std;

namespace csmp {

const double AZURE		   =1.0e-12;		/*  26, range 1-44 RGB 256 */
const double  BLUE		   =1.0e-13;		/*  62, range 44-73 */
const double  AQUAMARINE =1.0e-14;		/*  85, range 74-124 */
const double  GREEN		   =1.0e-15;		/* 164, range 125-181 */
const double  YELLOW		 =1.0e-16;		/* 199, range 182-204 */
const double  PEACH		   =1.0e-17;		/* 210, range 205-214 */
const double  ORANGE		 =1.0e-18;		/* 220, range 215-235 */
const double  RED	       =1.0e-19;		/* 252, range 236-254 */
const double  WHITE		   =1.0e-20;		/*   0 */
const double  BLACK		   =1.0e-21;		/* 255 */

/**

Function converts 256 colors (rainbow scheme)
to permeability values ranging from 1nD m2 to 1D.
'sc' is scaling factor for permeability, if
sc = 1.0 above scale applies.
*/
void convertColorToPermeability( double sc, Matrix& perm )
 {
   const size_t m = perm.Rows(), n = perm.Cols();

   for ( size_t i{0U}; i <m; i++)
     for ( size_t j{0U}; j <n; j++)
       {
 		    /* COLOR - PERMEABILITY CONVERSION */ 
 		    if      ( perm(i,j) == 0   ) perm( i,j) = (sc * WHITE);
  		  else if ( perm(i,j) == 255 ) perm( i,j) = (sc * BLACK);
  		  else if ( perm(i,j) <   44 ) perm( i,j) = (sc * AZURE);    
  		  else if ( perm(i,j) <   74 ) perm( i,j) = (sc * BLUE);
  		  else if ( perm(i,j) <  125 ) perm( i,j) = (sc * AQUAMARINE);
  		  else if ( perm(i,j) <  182 ) perm( i,j) = (sc * GREEN);
  		  else if ( perm(i,j) <  205 ) perm( i,j) = (sc * YELLOW);
  		  else if ( perm(i,j) <  215 ) perm( i,j) = (sc * PEACH);
  		  else if ( perm(i,j) <  236 ) perm( i,j) = (sc * ORANGE);
  		  else if ( perm(i,j) >= 236 ) perm( i,j) = (sc * RED);
  		  else {
  		        std::cerr <<"\nWARNING, convertColorToPermeability: ";
              std::cerr <<" conversion value not identified." << std::endl;
              perm( i,j) = std::numeric_limits<double>::quiet_NaN();
           }
       }
                                   
 } // end convertColorToPermeability


} // end namespace csmp
