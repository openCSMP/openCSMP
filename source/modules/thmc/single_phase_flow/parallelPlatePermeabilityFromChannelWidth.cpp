#include "parallelPlatePermeabilityFromChannelWidth.h"
#include "Region.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "PropertyHandle.h"
#include "NumIntegral_dNT_dN_dV.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "LinearSolver.h"
#include "ExtractVectorVariableLength.h"
#include "CSMP_definitions.h"

using namespace std;

namespace csmp {


/*! \file parallelPlatePermeabilityFromChannelWidth.h */

/**
     Solves the Laplace equation with a source term of 1 to find find the maximum
     transect across the target region by means of differentiation.
     
     @author SKM 8/2005
*/
template<uint32_t dim>
void parallelPlatePermeabilityFromChannelWidth( Model<dim>& sg, 
                                                const char* channel_region, 
                                                const char* channel_width,
                                                double minimum_channel_width )
 {
    if ( sg.Database().Type(channel_width) != SCALAR or 
         sg.Database().Placement(channel_width) != NODE )
         throw csmp::Exception( ERROR, "parallelPlatePermeabilityFromChannelWidth: property '", 
                                channel_width, "' must be a scalar node property." );

    PropertyHandle<dim>  src( sg, "source term", SCALAR, ELEMENT ); src = 1.; // gives integral_x1_x2 d^3 / 12 mu
    PropertyHandle<dim>  lap( sg, "parabolic function", SCALAR, NODE );
    lap = 1.0e-30;

    #ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Solver solver;
    PDE_Integrator<dim,Element>  parabolic_profile( solver );
    #else
    /// add extra functionality for alternative solver if needed
    CSMP_DEFAULT_LINEAR_SOLVER solver;
    PDE_Integrator<dim,Element>  parabolic_profile( solver );
    #endif
  
    // the maximum computed 'parabolic function' value is equivalent to the pore radius of the corresponding pore space segment
    // laplacian matrix [L] on the left-hand side
    NumIntegral_dNT_dN_dV<dim>  laplacian( sg.Database(), "parabolic function", "parabolic function" );
    // source vector {q} on the right-hand side = 1
    NumIntegral_NT_op_N_dV<dim> rhs( sg.Database(),"source term", "parabolic function" );
    parabolic_profile.Add( &laplacian );
    parabolic_profile.Add( &rhs );

    // boundary conditions: we set value of parabolic function at walls of pore space to zero=1e-50 
    // and fix this value by treating it as a Dirichlet boundary condition
    Region<dim>&  gref(sg.Region(channel_region));
    gref.InputPropertyValue( "parabolic function", makeScalar(DIRICH,1.0e-30), PERIMETER );
  
    // calculation
    parabolic_profile.IntegrateOver( sg, gref );
    parabolic_profile.Reset();
    printRangeOfVariable( sg, channel_region, "parabolic function" );
  
    // the variable 'parabolic function' is now mapped to the variable permeability. This is possible since we 
    // have verified that the ensuing flow matches the parallel plate approximation k = d^2 / 12
    csmp::Index  k_key = sg.Database().StorageKey("permeability");
    ScalarVariable  sc;
    for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ ) {
         (*eit)->PropertyValueAtBaryCenter( lap.Key(), sc );
         (*eit)->Store( k_key, sc );
      }

    // imposing a lower bound on channel permeability according to user-specified minimum channel width
    const bool use_upper_limit(false);
    const double min_channel_k = (minimum_channel_width * minimum_channel_width) / 12.;
    imposeLimitOn( sg, channel_region, "permeability", use_upper_limit, min_channel_k );

    printRangeOfVariable( sg, channel_region, "permeability" );
    printRangeOfVariable( sg, "permeability" );

    PropertyHandle<dim>  grad( sg, "grad parabolic function", VECTOR, ELEMENT );
    PropertyHandle<dim>  magn( sg, "magnitude grad parabolic function", SCALAR, ELEMENT ); 

    // (fracture) cavity radius = 1/2-aperture (=derivative of "parabolic function")
    sg.CopyGradientOfProperty_A_To_B( "parabolic function", "grad parabolic function" );
    ExtractVectorVariableLength<dim>  grad_magnitude( sg.Database(), "grad parabolic function", "magnitude grad parabolic function" );
    gref.Apply( grad_magnitude );
//    sg.InputPropertyValue( channel_width, makeScalar(PLAIN,0.) );
    gref.ExtrapolateCellToNodeProperty( "magnitude grad parabolic function", channel_width );
    PropertyHandle<dim>  cwidth( sg, channel_width ); cwidth *= 2.;
//    printRangeOfVariable( sg, channel_region, channel_width );
   
    // limiting the channel width to user specified value
    imposeLimitOn( sg, channel_region, channel_width, use_upper_limit, minimum_channel_width );
    printRangeOfVariable( sg, channel_region, channel_width );

 } // end parallelPlatePermeabilityFromChannelWidth
 
template void parallelPlatePermeabilityFromChannelWidth( Model<2U>&, const char*, const char*, double );
template void parallelPlatePermeabilityFromChannelWidth( Model<3U>&, const char*, const char*, double );



} // end namespace csmp

