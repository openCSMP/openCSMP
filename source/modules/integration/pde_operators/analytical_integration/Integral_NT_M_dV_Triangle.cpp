// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_NT_M_dV_Triangle.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {


/** For the mapping between a finite volume and linear finite element
discretization.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
Integral_NT_M_dV_Triangle<dim,CELL>::Integral_NT_M_dV_Triangle( const PropertyDatabase<dim>& pref,
                                                                const char* test )
  : MathOperatorRHS<dim,CELL>(pref,"permeability",test),
    NPROP(3)
 {
    MathOperatorRHS<dim,CELL>::Name("Integral_NT_M_dV_Triangle", "__", test );
    
        // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || 
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "Integral_NT_M_dV_Triangle<dim>::(constructor)", 
                             test, "Mapped variable must be a scalar property placed on the nodes." );
 }





/** Reads the Operand values from the elements.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_M_dV_Triangle<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    // this integral is only for analytically integrated finite elements
    assert( e.FE()->UsesLocalCoordinates() == false );

   // reading Young's modulus (must be an element variables)
   e.NodePropertyVector( MathOperatorRHS<dim,CELL>::TestOperandKey(), NPROP );
}




/** Computes the volume (area) integral over the testfunction products
multiplied with the Operand.  


@section implementation Implementation

N Phi_CV = sum_j Ni(at area barycenter) * area_j * Phi_i
*/
template<uint32_t dim, template<uint32_t> class CELL>
void Integral_NT_M_dV_Triangle<dim,CELL>::ComputeContribution( const CELL<dim>& e )
  {
     MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());

     Point<dim>     ctr(e.BaryCenter());
     vector<double> IPOL(e.Nodes());

     // sector surrounding node 0
     // -------------------------
     Point<dim>  xyz(ctr);
     // point average of node 0, sgm_mid_points 01, 20 and barycenter
     xyz += e.N(0)->Coordinate();
     // segment midpoint 01
     xyz += (e.N(0)->Coordinate() + e.N(1)->Coordinate()) / 2.;
     // segment midpoint 20
     xyz += (e.N(0)->Coordinate() + e.N(2)->Coordinate()) / 2.;
     // averaging coordinates to find barycenter
     xyz /= 4.; 
     // evaluating test functions at barycenter of quadrilateral sector
     e.N_AtGlobalPoint( IPOL, xyz.Coordinates() );
     // assigning nodal contribution of M-sector integral 1
     //                                Ni        Mj        
     MathOperatorRHS<dim,CELL>::RHS[0] = IPOL[0] * NPROP[0]();
     MathOperatorRHS<dim,CELL>::RHS[1] = IPOL[1] * NPROP[0]();
     MathOperatorRHS<dim,CELL>::RHS[2] = IPOL[2] * NPROP[0]();

     // sector surrounding node 1
     // -------------------------
     xyz  = ctr;
     xyz += e.N(1)->Coordinate();
     // segment 01
     xyz += (e.N(0)->Coordinate() + e.N(1)->Coordinate()) / 2.;
     // segment 12
     xyz += (e.N(1)->Coordinate() + e.N(2)->Coordinate()) / 2.;
     xyz /= 4.; 
     e.N_AtGlobalPoint( IPOL, xyz.Coordinates() );
     MathOperatorRHS<dim,CELL>::RHS[0] += IPOL[0] * NPROP[1]();
     MathOperatorRHS<dim,CELL>::RHS[1] += IPOL[1] * NPROP[1]();
     MathOperatorRHS<dim,CELL>::RHS[2] += IPOL[2] * NPROP[1]();

     // sector surrounding node 2
     // -------------------------
     xyz  = ctr;
     xyz += e.N(2)->Coordinate();
     // segment 12
     xyz += (e.N(1)->Coordinate() + e.N(2)->Coordinate()) / 2.;
     // segment 20
     xyz += (e.N(0)->Coordinate() + e.N(2)->Coordinate()) / 2.;
     xyz /= 4.; 
     e.N_AtGlobalPoint( IPOL, xyz.Coordinates() );
     MathOperatorRHS<dim,CELL>::RHS[0] += IPOL[0] * NPROP[2]();
     MathOperatorRHS<dim,CELL>::RHS[1] += IPOL[1] * NPROP[2]();
     MathOperatorRHS<dim,CELL>::RHS[2] += IPOL[2] * NPROP[2]();
     
     // Evaluating the M-Volume integrals for the FV solution
     vol_div3 = e.Volume() / 3.;
     MathOperatorRHS<dim,CELL>::RHS[0] *= vol_div3;
     MathOperatorRHS<dim,CELL>::RHS[1] *= vol_div3;
     MathOperatorRHS<dim,CELL>::RHS[2] *= vol_div3;

} // end ComputeContribution



template class Integral_NT_M_dV_Triangle<2U,Element>;
template class Integral_NT_M_dV_Triangle<3U,Element>;

template class Integral_NT_M_dV_Triangle<2U,Face>;
template class Integral_NT_M_dV_Triangle<3U,Face>;

} // csmp
