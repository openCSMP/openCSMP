// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NUM_INTEGRAL_DNT_OP_DN_DN_INTERFACE_H
#define NUM_INTEGRAL_DNT_OP_DN_DN_INTERFACE_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "InterFace.h"

namespace csmp {

class FiniteElement;
template<uint32_t> class Model;

/**
       Couples the nodes on the opposite sides of the InterFace,
       using an operand (thermal or hydraulic conductivity) and the thickness attribute "thickness" for the interface.
       (the latter needs to be defined as a separate variable with the placement INTER_FACE).
       
       The coupling is created by emulating a thin-shell element (prism or hexahedron).
       For this purpose the opposing collocated nodes are separated from one-another an amount equivalent to the thickness
       before the coordinate matrix is initialised.
       
       The separation of the nodes creates prism or hexahedral elements, from triangular and quadrilateral InterFace objects.
       The InterFace normal is used in the construction of the edges of these elements.
       
       @author SKM
       @date 18/8/22
*/
template<uint32_t dim>
class NumIntegral_dNT_op_dN_dV_InterFace : public MathOperatorLHS<dim,InterFace> {
  public:
    NumIntegral_dNT_op_dN_dV_InterFace( const Model<dim>&,
                                        const char* oper,
                                        const char* basic,
                                        const char* test );
    
    void GetOperands( const InterFace<dim>& ) override final;
    void ComputeContribution( const InterFace<dim>& ) override final;
    
  protected:
      /// separates opposing nodes and writes new coordinates to XY matrix of corresponding finite element type
      void InitialiseCoordinateMatrix( const InterFace<dim>&, DenseMatrix<DM_MIN>& ) const;
  
  private:
      csmp::INDEX<SCALAR,INTER_FACE> thi_key_;            ///< 'thickness' parameter for the interface
      double                         thickness_;          ///< zone thickness modelled by interface
      FiniteElement*                 hexa_ptr_ = nullptr; ///< element type used for quadrilateral Interfaces
      FiniteElement*                 pris_ptr_ = nullptr; ///< element type for triangles

      DenseMatrix<DM_MIN>  DN_, DNT_;  ///< matrices needed for the accumulation
};

} // csmp

#endif
