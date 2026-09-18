// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h"
//#include "UpwindControlVisitor.h"
//#include "ExplicitFiniteVolumeTransportPHX.h"
#include <iostream>
#include <chrono>

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::~CVFEM_Upwind_NumIntegral_dNT_op_dN_dV() {}


template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::CVFEM_Upwind_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                                        UpwindControlVisitor<dim>& upwind_visitor,
                                                                                        ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                                                                        const char* oper,
                                                                                        const char* basic,
                                                                                        const char* test,
                                                                                        const char* upwind,
                                                                                        const char* grav_trigger,
                                                                                        const char* thickness)//Benoit 2025 add

    : CVFEM_MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
      UpwindVisitor(upwind_visitor),
      finite_volume(fv_transport),
      B(dim,3),
      uvar_(pref.StorageKey(upwind)),
      gtvar_(pref.StorageKey(grav_trigger)),
      thickness_(pref.StorageKey(thickness))//Benoit 2025 add

{

    string name = "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV ";
    name += upwind;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim,CELL>::Name(cname, oper, basic, test);

    // testing the Operands

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT )
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               oper, "Operand must be placed on the element.");
    if ( MathOperatorLHS<dim,CELL>::BasicOperandType() == SCALAR )
        operand_is_tensor = false;
    else if ( MathOperatorLHS<dim,CELL>::BasicOperandType() == TENSOR )
        operand_is_tensor = true;
    else
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               oper, "Operand must be a scalar or tensor property.");

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
        throw csmp::Exception( ERROR,  "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
    
    if (uvar_.place != NODE || uvar_.type != SCALAR)
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               upwind, "Upwind variable must be a scalar property placed on the nodes.");

    if (gtvar_.place != NODE || gtvar_.type != SCALAR)
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)",
                               grav_trigger, "Upwind variable must be a scalar property placed on the nodes.");

}

template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    MathOperatorLHS<dim,CELL>::GetOperands(e);

    if (uvar_.place == NODE && uvar_.type == SCALAR)
        e.NodePropertyVector( uvar_, el_uvar);
    else
    {
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperands",
                               "Only nodal properties allowed" );
    }
} // end GetOperands


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::GetOperandsCVFEM( const CELL<dim>& e,
                                                                        csmp::Index upwind_var_key )
{
    MathOperatorLHS<dim,CELL>::GetOperands(e);

    if (uvar_.place == NODE && uvar_.type == SCALAR &&
            upwind_var_key.place == NODE && upwind_var_key.type == SCALAR)
    {
        e.NodePropertyVector( uvar_, el_uvar);
        e.NodePropertyVector( upwind_var_key, upwind_var_multiplier);
    }
    else {
        throw csmp::Exception( ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperands",
                               "Only nodal properties allowed" );
    }

    for (size_t i = 0; i < el_uvar.size(); i++)
        el_uvar[i]() *= upwind_var_multiplier[i]();

} // end GetOperandsCVFEM


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    GetUpwindMatrix( e );

    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();

    if ( operand_is_tensor )
    {
        operand_tensor     = 0;

        if constexpr(dim==1U) // JK added: constexpr Loop will jump to the correct if condition at runtime.
        {
            //permeability tensor diagonal components
            operand_tensor(0,0)= MathOperatorLHS<dim, CELL>::MTRL[0](0,0);
        }
        else if constexpr(dim==2U)
        {
            //permeability tensor diagonal components
            operand_tensor(1,1)= MathOperatorLHS<dim, CELL>::MTRL[0](1,1);
            operand_tensor(0,0)= MathOperatorLHS<dim, CELL>::MTRL[0](0,0);
        }
        else if constexpr(dim==3U)
        {
            //permeability tensor diagonal components
            operand_tensor(1,1)= MathOperatorLHS<dim, CELL>::MTRL[0](1,1);
            operand_tensor(0,0)= MathOperatorLHS<dim, CELL>::MTRL[0](0,0);
            operand_tensor(2,2)= MathOperatorLHS<dim, CELL>::MTRL[0](2,2);
        }
    }

    else
    {
        operand = MathOperatorLHS<dim,CELL>::MTRL[0](0,0);
    }

    // loop over facets of the element
    for ( auto i{0U}; i<e.FV()->Facets(); i++ )
    {
        if ( operand_is_tensor )
        {
            e.dN(BB);
            B = operand_tensor;
            B *= BB;
        }

        else
        {
            e.dN(B);
            B *= operand;
        }

        e.FV()->FacetEdgeNodes( i, inside_node_, outside_node_ );

        area = finite_volume.GetFacetArea( e.Idx(), i );
        area*=e.Read(thickness_);//Benoit 2025 add


        for ( uint32_t d = 0; d < dim; d++)
        {
            normal_component = finite_volume.GetFacetNormalComponent( e.Idx(), i, d );

            for ( uint32_t n = 0; n < e.Nodes(); n++)
            {
                contribution = B(d,n)*area*normal_component;
                if      (upwind(inside_node_, outside_node_) == 1) contribution *= el_uvar[inside_node_]();
                else if (upwind(inside_node_, outside_node_) == 2) contribution *= el_uvar[outside_node_]();
                else contribution = 0.0;
                MathOperatorLHS<dim,CELL>::LHS(inside_node_,n)  -= contribution;
                MathOperatorLHS<dim,CELL>::LHS(outside_node_,n) += contribution;
            }
        }
    }
} // end ComputeContribution


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim,CELL>::GetUpwindMatrix( const CELL<dim>& e )
{
    upwind = UpwindVisitor.UpwindMatrix( gtvar_, e.Idx() );
} // GetUpwindMatrix


template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<1>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<2>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<3>;

template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<1,Face>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<2,Face>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<3,Face>;

} // csmp
