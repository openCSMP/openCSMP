// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PermeabilityTensor_Example.h"
#include "ScalarVariable.h"
#include "TensorVariable.h"
#include "VectorVariable.h"

#define DIM 3U

using namespace std;

namespace csmp{

void PermeabilityTensor_Example::Specifications()
{
  SetTitle( "Permeability tensor calculations" );
  SetDifficulty( 1 );
  SetCategory( "Software Functionality" );
  AddAuthor( "P. Lang & S. Azizmohammadi" );
  AddDescription( "source in: PermeabilityTensor_Example.cpp" );
  AddDescription( "General operations with csmp variable types scalar, vector, and tensor" );
  AddDescription( "Darcy flow calculations using a full tensor permeability" );
  AddDescription( "and pressure gradient" );
}


/** DARCY FLOW CALCULATIONS USING A PERMEABILITY TENSOR AND PRESSURE GRADIENT

   | qx |    A    | Kxx Kxy Kxz |  | dPHI / dX |
   | qy | =  --   | Kyx Kyy Kyz |  | dPHI / dY |
   | qz |    mu   | Kzy Kzx Kzz |  | dPHI / dZ |

*/
void PermeabilityTensor_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  // creating an empty flow rate vector( 3D - template parameter)
  VectorVariable<DIM> q;
  // which gives us a flow vector with qx, qy, qz

  // creating ScalarVariables for cross sectional area and viscosity respectively
  ScalarVariable a( PLAIN, 1. ), mu( PLAIN, 1.0e-3 );
  //                                        ^^^^^^ CSMP uses SI units, e.g. Pa.s

  // now we prepare for establishing the permeability tensor by creating a
  // cache variable, which will represent our horizontal permeability
  double kHorizontal( 1.0e-12 );
  //                    ^^^^^^ [m2] ~ 1 Darcy

  // now we establish the permeability tensor, filled with (PLAIN, 0.) values
  TensorVariable<DIM> kTensor( PLAIN, 0. );

  // we set the Kxx and Kzz permeability
  kTensor( 0,0 ) = kHorizontal;
  kTensor( 2,2 ) = kHorizontal;
  // the vertical permeability is set to a third of the horizontal
  kTensor( 1,1 ) = kHorizontal / 3.;

  // next we create the potential gradient vector of our conceptual flow
  VectorVariable<DIM> grad( PLAIN, 0. );
  //                         ^^^^^^^^^ initialize all values
  // dPHI / dX
  grad( 0 ) = 10.;
  // dPHI / dY
  grad( 1 ) = 0.;
  // dPHI / dZ
  grad( 2 ) = 0.1;

  // calculating flow vector
  for( uint32_t i = 0; i < 3; ++i )
    q( i ) = a() / mu() * (   kTensor( i,0 )*grad( 0 )
                                        + kTensor( i,1 )*grad( 1 )
                                        + kTensor( i,2 )*grad( 2 ) );
  //           ^^^^^        ^^^^^                 ^^^^^             ^^^^^
  // the Value() function of a variable returns its value, but does not allow
  // to manipilate it. ( in contrast to variable() which allows manipulation,
  // see above

  // Outputting results
  cout << "\nExample: Cross sectional area: " << a() << endl;
  cout << "\nExample: Fluid viscosity: " << mu() << endl;
  cout << "\nExample: Permeability tensor: \n";
  kTensor.Out();
  cout << "\nExample: Potetntial gradient: \n";
  grad.Out();
  cout << "\nExample: Resulting flow rates: \n";
  q.Out();

} // end Run

} // csmp

