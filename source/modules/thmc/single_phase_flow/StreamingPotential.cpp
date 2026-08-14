#include "StreamingPotential.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "PDE_Integrator.h"
#include "NumIntegral_dNT_lhsop_dN_dV.h"
#include "NumIntegral_dNT_rhsop_dN_dV.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
PotentialSource<dim>::PotentialSource( const PropertyDatabase<dim>& p ) 
  : Interrelation<dim>(p),
    mu( Interrelation<dim>::GlobalProperty("fluid viscosity") ),
    k( Interrelation<dim>::GlobalProperty("permeability") ),
    rhof( Interrelation<dim>::GlobalProperty("fluid density") ),
    R( Interrelation<dim>::GlobalProperty("R factor") )
 {
    Interrelation<dim>::Name("PotentialSource");
    Interrelation<dim>::OutputCondition( R, PLAIN );
    Interrelation<dim>::ResultProperty("R factor");
 } 



template<uint32_t dim>
void PotentialSource<dim>::Calculate()
 {
    // numerator of R factor, R = mu / (k rhof)
    mu.AssignTo( temp1 );
    
    R = -temp1(); // - since the C-factor is negative
    
    // denominator of r factor
    k.AssignTo( temp1 );
    rhof.AssignTo( temp2 );
    
    R /= temp1() * temp2();
 } 


template class PotentialSource<2U>;
template class PotentialSource<3U>;




/// the streaming potential method
template<uint32_t dim>
StreamingPotential<dim>::StreamingPotential( Model<dim>& sg )
 : ones( sg, "one coefficient", SCALAR, ELEMENT ),
   ccoeff( sg, "coupling coefficient", SCALAR, ELEMENT )
 {
 }
 
 
template<uint32_t dim>
StreamingPotential<dim>::~StreamingPotential()
 {
 }
    
 
/**

0 = div^2 phi  -  c div^2 p

an arbitrary point somewhere in the model is used as a reference "electrode"
it is assigned a Dirichlet potential value of zero.
*/
template<uint32_t dim>
void StreamingPotential<dim>::EvaluatePotential( Model<dim>& sg,
                                                 const char* potential )
 {
    // 1. settin the coupling coefficient to a uniform value of 
    ones   =  1.0;
    ccoeff = -1.0e-6; // V Pa-1
  
    // 2. computing R coefficients for equation
//    PotentialSource<dim>  potential_source( sg.Database() );
//    sg.Pass( potential_source );
    
    // 3. computation of the streamling potential
    PDE_Integrator<dim,Element>  streaming_potential;
    //streaming_potential.IncreaseMultiGridVectorStorage( 10 );
  
    NumIntegral_dNT_lhsop_dN_dV<dim>     lap2( sg.Database(),
                                            "one coefficient",
                                            potential,
                                            potential );
                                            
    NumIntegral_dNT_rhsop_dN_dV<dim>  c_lap2p( sg.Database(),
                                              "coupling coefficient",
                                              "fluid pressure",
                                               potential );
    
    streaming_potential.Add( &lap2 );
    streaming_potential.Add( &c_lap2p );
    
    sg.Apply( streaming_potential );

 } // end EvaluatePotential
    


template class StreamingPotential<1U>;
template class StreamingPotential<2U>;
template class StreamingPotential<3U>;

} // csmp


