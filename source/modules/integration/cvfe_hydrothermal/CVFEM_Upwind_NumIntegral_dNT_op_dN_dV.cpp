#include "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV.h"

using namespace std;

namespace csmp {

template<size_t dim>
CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::~CVFEM_Upwind_NumIntegral_dNT_op_dN_dV() {}

template<size_t dim>
CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::CVFEM_Upwind_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref, 
                                                                            UpwindControlVisitor<dim>& upwind_visitor,
                                                                            ExplicitFiniteVolumeTransportPHX<dim>& fv_transport,
                                                            				const char* oper, 
                                                            				const char* basic, 
                                                            				const char* test,
                                                            				const char* upwind,
                                                            				const char* grav_trigger)
 
  : CVFEM_MathOperatorLHS<dim>(pref,oper,basic,test),
    UpwindVisitor(upwind_visitor),
    finite_volume(fv_transport),
    B(dim,3),
    uvar_(pref.StorageKey(upwind)),
    gtvar_(pref.StorageKey(grav_trigger))
{

    string name = "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV ";
    name += upwind;
    char * cname = new char[name.length()+1];
    strcpy(cname, name.c_str());
    MathOperatorLHS<dim>::Name(cname, oper, basic, test);

    // testing the Operands 

    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT )
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)", 
                    oper, "Operand must be placed on the element.");

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR,  "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)", 
                    test, "Basic (dependent) variable must be a scalar property placed on the nodes.");

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || MathOperatorLHS<dim>::TestOperandType() != SCALAR )
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
    
    if (uvar_.place != NODE || uvar_.type != SCALAR)
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)", 
                    upwind, "Upwind variable must be a scalar property placed on the nodes.");
                    
    if (gtvar_.place != NODE || gtvar_.type != SCALAR)
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV::(constructor)", 
                    grav_trigger, "Upwind variable must be a scalar property placed on the nodes.");
                    
}

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperands( Element<dim>& e )
 {


 	MathOperatorLHS<dim>::GetOperands(e);
      
  if (uvar_.place == NODE && uvar_.type == SCALAR)
    e.NodePropertyVector( uvar_, el_uvar);
  else {
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperands", 
                                      "Only nodal properties allowed" );
  }

    
 } // end GetOperands

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperandsCVFEM( Element<dim>& e,
                                          				      csmp::Index upwind_var_key )
 {


  MathOperatorLHS<dim>::GetOperands(e);
      
  if (uvar_.place == NODE && uvar_.type == SCALAR && 
      upwind_var_key.place == NODE && upwind_var_key.type == SCALAR)
    {
     e.NodePropertyVector( uvar_, el_uvar);
     e.NodePropertyVector( upwind_var_key, upwind_var_multiplier);
    }
  else {
    throw csmp::Exception( CSMP_ERROR, "CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetOperands", 
                                      "Only nodal properties allowed" );
  }

 for (size_t i = 0; i < el_uvar.size(); i++)
    el_uvar[i]() *= upwind_var_multiplier[i]();
    
 } // end GetOperandsCVFEM

template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::ComputeContribution( Element<dim>& e )
 {
 
    GetUpwindMatrix( e );
    
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

       e.dN( B );
       operand = MathOperatorLHS<dim>::MTRL[0](0,0);
//       B *= operand;
//    Point<dim> facet_point;
//    std::vector<csmp_float> facet_point_vector;
//   csmp_float detJ;

    // loop over facets of the element
    for ( size_t i=0U; i<e.FV_Stencil()->Facets(); i++ )
       {

//       detJ = e.dN( B, e.ConnectedFiniteVolumeStencil()->FacetEdgeMidPoint(i).Coordinates() );

       e.dN_AtIntegrationPoint( B, i, 1 );
       //B.Out();

       B *= operand;
    

          e.FV_Stencil()->FacetEdgeNodes( i, inside_node_, outside_node_ );
          area = finite_volume.GetFacetArea( e.Idx(), i );
          
          for (size_t d = 0; d < dim; d++)
           {
            normal_component = finite_volume.GetFacetNormalComponent( e.Idx(), i, d );
            for (size_t n = 0; n < e.Nodes(); n++)
               {
                 contribution = B(d,n)*area*normal_component;
    			 if (upwind(inside_node_, outside_node_) == 1) contribution *= el_uvar[inside_node_]();
    			 else if (upwind(inside_node_, outside_node_) == 2) contribution *= el_uvar[outside_node_]();
    			 else contribution = 0.0;
                 MathOperatorLHS<dim>::LHS(inside_node_,n)  -= contribution;
                 MathOperatorLHS<dim>::LHS(outside_node_,n) += contribution;
               }
            }
       }


} // end ComputeContribution


template<size_t dim>
void CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<dim>::GetUpwindMatrix( Element<dim>& e )
{

  upwind = UpwindVisitor.UpwindMatrix( gtvar_, e.Idx() );    

} // GetUpwindMatrix


template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<1>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<2>;
template class CVFEM_Upwind_NumIntegral_dNT_op_dN_dV<3>;

} // csmp
