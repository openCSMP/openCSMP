#include "FiniteElementVariablePlacement_Example.h"

#include "ANSYS_Model2D.h"
#include "InputDataManager.h"
#include "PropertyHandle.h"

#include "FiniteVolumeStencil.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

#include "Variables_TracerTransfer.h"
#include "VTK_Interface.h"

#define DIM 2U

using namespace std;

namespace csmp {

void FiniteElementVariablePlacement_Example::Specifications()
    {
        SetTitle( "Finite Element Variable Placement Test" );
        SetDifficulty( 2 );
        SetCategory( "Software Functionality" );
        AddAuthor( "Mahyar Madadi" );
        AddDescription( "source in: FiniteElementVariablePlacement_Example.cpp" );
    }


void FiniteElementVariablePlacement_Example::Run()
{
    // reading the model name
    string modelName = "xsection";
  
    // reading 2D ANSYS model and its variable
  
    ANSYS_Model2D  model( modelName.c_str(), "CSMP-Operand_Test-variables.txt" );
  
    model.InstantiateFiniteVolumes(); // SKM comment: If you do not have any finite volume variables in your property database file, the stencils will not be created automatically.
                                      // Therefore  this = NULL, in other words, there is no finite volume stencil connected to your FE and therefore you cannot access that.
  
    // reference to the Master Region in the model; this Region contains  all the elements and nodes.
    Region<2U>& rref = model.Region( "Model" );
  
  
  /* Summary: by Mahyar 6/July/2018
     In this part .. we check and compare two scheme of interpolations from Nodes to Bary Center, IP, Sector IPs and Facet IPs...
     The new scheme works well for the Bary Center, Sector IP's and FAcet IP's but it is wrong for the IP...
     The bug reported to SKM.
  */
  
  
  // 1) Initialization of values: visosity oil (Scalar) and nodal velocity (Vector)
  
    double64 viscosity_oil = 4.0e-4 ;
    double64 vx = 0.01 ;
    double64 vy = 0.01 ;
    double64 pi = 4.0*atan(1.0) ;
  
    double64 sxx  = 0.01 ;
    double64 syy  = 0.01 ;
    double64 sxy  = 0.01 ;

  
    // define the viscosity of oil as a scalar in the model with constant value
    model.InputPropertyValue("nodal scalar variable", makeScalar(PLAIN, viscosity_oil)) ;  // oil vicosity at 20 degree
    model.InputPropertyValue("nodal vector variable",makeVector(PLAIN, PLAIN, vx, vy)) ;
    model.InputPropertyValue("nodal tensor variable",makeTensor(PLAIN, PLAIN, sxx, sxy, syy, sxy));

  

    // define the key parameter for the "viscosity of oil"
    const csmp::INDEX< SCALAR, NODE> key_node_Parameter (model.Database().StorageKey("nodal scalar variable"));
    const csmp::INDEX< VECTOR, NODE> key_node_VectorParameter (model.Database().StorageKey("nodal vector variable"));
    const csmp::INDEX< TENSOR, NODE> key_node_TensorParameter (model.Database().StorageKey("nodal tensor variable"));

     double64 i = 0 ;
    // printing the whole values for the viscosity
    for (  auto it = rref.NodesBegin(); it != rref.NodesEnd(); ++it) {  // loop over nodes
      
      i++ ;
      ScalarVariable variable1_ ((*it)->Status(key_node_Parameter), i*viscosity_oil/100.0) ;
      (*it)->Store( key_node_Parameter, variable1_);
      
      cerr <<"\n" ;
      cout <<"viscosity of oil= "<<variable1_<<"\n" ;
      
      auto xcoord = (*it)->x() ;
      auto ycoord = (*it)->y() ;
      
      // MAHYAR Question: how can I use the Status(key_node_VectorParameter) for vector?
      VectorVariable<2U> vector_variable1_ (PLAIN, PLAIN, rand()/((double64) RAND_MAX)*vx+cos(xcoord*pi/100.), rand()/((double64) RAND_MAX)*vy);
      (*it)->Store( key_node_VectorParameter, vector_variable1_);
      
      cerr <<"\n" ;
      cout <<"nodal velocity= "<<vector_variable1_<<"\n" ;
      
      TensorVariable<2U> tensor_variable1_ (PLAIN, PLAIN, PLAIN, sxx*xcoord, sxy*xcoord*ycoord, syy*ycoord);
      (*it)->Store( key_node_TensorParameter, tensor_variable1_);
      
    } // end of loop over nodes
  
    // visualization of input parameter:
  
    VTK_Interface<2U>  vtk_output;
    vtk_output.OutputDataToVTK( model, "nodal scalar variable", "nodal scalar variable", 0 );
    vtk_output.OutputDataToVTK( model, "nodal vector variable", "nodal vector variable", 0 );
    vtk_output.OutputDataToVTK( model, "nodal tensor variable", "nodal tensor variable", 0 );  // Mahyar: the output is 3x3 Matrix rather than 2x2
  
    // estimating the viscosity at the  Bary Center and IP
    cerr <<"/n" ;
    cerr <<"----- The Comparison of Old and New method of interpolations of Node value to, Bary center, IP's -----\n" ;

    for ( auto it = rref.ElementsBegin(); it != rref.ElementsEnd(); ++it ){  // loop over all elements
        size_t idx_ = (*it)->Idx()  ;
      
        cerr <<"Entering to Element #"<<idx_<<"\n" ;
        cerr <<" \n" ;
      
        int idy_ = 0 ;
        for ( auto jt : (*it)->AllNodes() ){     //loop over IP in the element
          auto variable0_ = (jt).Read(key_node_Parameter) ;
          cerr <<"  Input values of viscosity of oil at Node "<< idy_ <<", is "<<variable0_<<"\n" ;
          
          auto vector_variable0_ = (jt).Read(key_node_VectorParameter) ;
          cerr <<"  Input values of nodal vector at Node "<< idy_ <<", is "<<vector_variable0_[0]<<" ,"<<vector_variable0_[1]<<"\n" ;
          
          TensorVariable<2U> tensor_variable0_ ;
          (jt).Read(key_node_TensorParameter, tensor_variable0_ ) ;
          cerr <<"  Input values of nodal tensor at Node "<< idy_ <<", is "<<tensor_variable0_(0,0)<<" ,"<<tensor_variable0_(1,0)<<" ,"<<tensor_variable0_(1,1)<<" ,"<<tensor_variable0_(0,1)<<"\n" ;
          
          cerr<<" \n" ;
          
          idy_++ ;
        }
        cerr <<" \n" ;
      
      
        //++++++++++++++++++++++++++++++++++  Scalar
        //  at the  Bary Center with new scheme
      
        auto E1 = (*it)->AtBarycenter() ;
        double64 variable2_ = E1.Obtain(key_node_Parameter) ;
      
        cerr <<"  New Method: viscosity of oil at Bary Center() = "<<variable2_<<"\n" ;
        cerr <<" \n" ;
      
        // compare with old scheme at the  Bary Center
      
        ScalarVariable oldvariable2_;
        (*it)->PropertyValueAtBaryCenter( key_node_Parameter, oldvariable2_ );

        cerr <<"  Old Method: viscosity of oil at Bary Center() = "<<oldvariable2_<<"\n" ;
        cerr <<" \n" ;
      
      
        //++++++++++++++++++++++++++++++++++ Vector
        //  at the  Bary Center with new scheme
        auto vector_variable2_ = E1.Obtain(key_node_VectorParameter) ;
      
        cerr <<"  New Method: nodal velocity at Bary Center() = "<<vector_variable2_[0]<<", "<<vector_variable2_[1]<<"\n" ;
        cerr <<" \n" ;
      
        // compare with old scheme at the  Bary Center
      
        VectorVariable<2U> oldvector_variable2_;
        (*it)->PropertyValueAtBaryCenter( key_node_VectorParameter, oldvector_variable2_ );
      
        cerr <<"  Old Method: nodal velocity at Bary Center() = "<<oldvector_variable2_<<"\n" ;
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++ Tensor
        //  at the  Bary Center with new scheme
        TensorVariable<2U> tensor_variable2_ ;
        E1.Obtain(key_node_TensorParameter , tensor_variable2_ );
      
        cerr <<"  New Method: nodal tensor at Bary Center() = "<<tensor_variable2_(0,0)<<" ,"<<tensor_variable2_(1,0)<<" ,"<<tensor_variable2_(1,1)<<"\n" ;
        cerr <<" \n" ;
      
        // compare with old scheme at the  Bary Center
      
        TensorVariable<2U> oldtensor_variable2_;
        (*it)->PropertyValueAtBaryCenter( key_node_TensorParameter, oldtensor_variable2_ );
      
        cerr <<"  Old Method: nodal tensor at Bary Center() = "<<oldtensor_variable2_(0,0)<<" ,"<<oldtensor_variable2_(1,0)<<" ,"<<oldtensor_variable2_(1,1)<<"\n" ;
        cerr <<" \n" ;
      
      
        //++++++++++++++++++++++++++++++++++  Scalar
        // at the IP with new scheme
        idy_ = 0 ;
        for ( auto jt : (*it)->AllElementIntegrationPoints() ){     //loop over IP in the element
          
          auto variable3_ = (jt).Obtain(key_node_Parameter) ;
          cerr <<"  New Method: viscosity of oil at IP No. "<< idy_ <<", is "<<variable3_<<"\n" ;
          idy_++ ;
          
        }
        cerr <<" \n" ;

        // at the  IP compare with old scheme
        auto N_IP = (*it)->IntegrationPoints() ;
        
        for ( size_t j=0; j< N_IP; ++j ) { //loop over integeration points
          ScalarVariable oldvariable3_;
          (*it)->PropertyValueAtIntegrationPoint( key_node_Parameter, j ,oldvariable3_ );
          cerr <<"  Old Method: viscosity of oil at IP No. "<< j <<", is "<<oldvariable3_<<"\n" ;
        }
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++ Vector
        // at the IP with new scheme
        idy_ = 0 ;
        for ( auto jt : (*it)->AllElementIntegrationPoints() ){     //loop over IP in the element
          auto vector_variable3_ = (jt).Obtain(key_node_VectorParameter) ;
          cerr <<"  New Method: nodal velocity at IP No. "<< idy_ <<", is "<<vector_variable3_[0]<<", "<<vector_variable3_[1]<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
      
        // at the  IP compare with old scheme
      
        for ( size_t j=0; j< N_IP; ++j ) {
          VectorVariable<2U> oldvector_variable3_;
          (*it)->PropertyValueAtIntegrationPoint( key_node_VectorParameter, j ,oldvector_variable3_ );
          cerr <<"  Old Method: nodal velocity at IP No. "<< j <<", is "<<oldvector_variable3_<<"\n" ;
        }
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++ Tensor
        // at the IP with new scheme
    
        idy_ = 0 ;
        for ( auto jt : (*it)->AllElementIntegrationPoints() ){     //loop over IP in the element
          TensorVariable<2U> tensor_variable3_ ;
          (jt).Obtain(key_node_TensorParameter , tensor_variable3_ );
          cerr <<"  New Method: nodal tensor at IP No. "<< idy_ <<", is "<<tensor_variable3_(0,0)<<" ,"<<tensor_variable3_(1,0)<<" ,"<<tensor_variable3_(1,1)<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
      
        // at the  IP compare with old scheme
      
        for ( size_t j=0; j< N_IP; ++j ) {
          TensorVariable<2U> oldtensor_variable3_;
          (*it)->PropertyValueAtIntegrationPoint( key_node_TensorParameter, j ,oldtensor_variable3_ );
          cerr <<"  Old Method: nodal tensor at IP No. "<< j <<", is "<<oldtensor_variable3_(0,0)<<" ,"<<oldtensor_variable3_(1,0)<<" ,"<<oldtensor_variable3_(1,1)<<"\n" ;
        }
        cerr <<" \n" ;
      
      
      
        //++++++++++++++++++++++++++++++++++  Scalar
        // at the Facet IP with new scheme
        idy_ = 0 ;
        for ( auto kt : (*it)->AllFacetIntegrationPoints() ){     //loop over Facet IP in the element
          auto variable4_ = (kt).Obtain(key_node_Parameter) ;
          cerr <<"  New Method: viscosity of oil at Facet IP No. "<< idy_ <<", is "<<variable4_<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
      
        // compare with old scheme at the  Facet IP
        for ( size_t facet=0; facet< (*it)->Facets(); ++facet )
          for ( size_t ip_f=0; ip_f < (*it)->IntegrationPointsPerFacet(); ++ip_f ) {
             double64 oldvariable4_= (*it)->PropertyValueAtFacetIntegrationPoint( facet, ip_f, key_node_Parameter );
             cerr <<"  Old Method: viscosity of oil at Facet IP No. "<< facet <<", is "<<oldvariable4_<<"\n" ;
          }
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++ Vector
        // at the Facet IP with new scheme
        idy_ = 0 ;
        for ( auto kt : (*it)->AllFacetIntegrationPoints() ){     //loop over Facet IP in the element
          auto vector_variable4_ = (kt).Obtain(key_node_VectorParameter) ;
          cerr <<"  New Method: nodal velocity at Facet IP No. "<< idy_ <<", is "<<vector_variable4_[0]<<", "<<vector_variable4_[1]<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
     
        // compare with old scheme at the  Facet IP
        for ( size_t facet=0; facet< (*it)->Facets(); ++facet )
          for ( size_t ip_f=0; ip_f < (*it)->IntegrationPointsPerFacet(); ++ip_f ) {
            VectorVariable<2U> oldvector_variable4_;
            (*it)->PropertyValueAtFacetIntegrationPoint( key_node_VectorParameter, facet, ip_f,  oldvector_variable4_ );
            cerr <<"  Old Method: nodal velocity at Facet IP No. "<< facet <<", is "<<oldvector_variable4_<<"\n" ;
          }
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++  Tensor
      
        // at the Facet IP with new scheme
        idy_ = 0 ;
        for ( auto kt : (*it)->AllFacetIntegrationPoints() ){     //loop over Facet IP in the element
          TensorVariable<2U> tensor_variable4_ ;
          (kt).Obtain(key_node_TensorParameter , tensor_variable4_ );
          cerr <<"  New Method: nodal tensor at Facet IP No. "<< idy_ <<", is "<<tensor_variable4_(0,0)<<" ,"<<tensor_variable4_(1,0)<<" ,"<<tensor_variable4_(1,1)<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
      
        // compare with old scheme at the  Facet IP
        for ( size_t facet=0; facet< (*it)->Facets(); ++facet )
          for ( size_t ip_f=0; ip_f < (*it)->IntegrationPointsPerFacet(); ++ip_f ) {
            TensorVariable<2U> oldtensor_variable4_;
            (*it)->PropertyValueAtFacetIntegrationPoint( key_node_TensorParameter, facet, ip_f,  oldtensor_variable4_ );
            cerr <<"  Old Method: nodal tensor at Facet IP No. "<< facet <<", is "<<oldtensor_variable4_(0,0)<<" ,"<<oldtensor_variable4_(1,0)<<" ,"<<oldtensor_variable4_(1,1)<<"\n" ;
          }
        cerr <<" \n" ;

        //++++++++++++++++++++++++++++++++++  Scalar
        // at the  Sector IP with new scheme
        for ( auto lt : (*it)->AllSectorIntegrationPoints() ){     //loop over Sector IP in the element
          auto variable5_ = (lt).Obtain(key_node_Parameter) ;
          auto idy_ = lt.NodeIdx() ;
          cerr <<"  New Method: viscosity of oil at Sector IP "<< idy_ <<", is "<<variable5_<<"\n" ;
        }
        cerr <<" \n" ;
        
        // compare with old scheme at the Sector IP
        for ( size_t sector=0; sector< (*it)->Sectors(); ++sector )
          for ( size_t ip_s=0; ip_s < (*it)->IntegrationPointsPerSector(); ++ip_s ) {
            double64 oldvariable5_= (*it)->PropertyValueAtSectorIntegrationPoint( sector, ip_s, key_node_Parameter );
            cerr <<"  Old Method: viscosity of oil at Sector IP No. "<< sector <<", is "<<oldvariable5_<<"\n" ;
          }
        
        cerr <<" \n" ;
      
      
        //++++++++++++++++++++++++++++++++++ Vector
        // at the  Sector IP with new scheme
        for ( auto lt : (*it)->AllSectorIntegrationPoints() ){     //loop over Sector IP in the element
          auto vector_variable5_ = (lt).Obtain(key_node_VectorParameter) ;
          auto idy_ = lt.NodeIdx() ;
          cerr <<"  New Method: nodal velocity at Sector IP "<< idy_ <<", is "<<vector_variable5_[0]<<", "<<vector_variable5_[1]<<"\n" ;
        }
        cerr <<" \n" ;
      
        // compare with old scheme at the Sector IP
        for ( size_t sector=0; sector< (*it)->Sectors(); ++sector )
          for ( size_t ip_s=0; ip_s < (*it)->IntegrationPointsPerSector(); ++ip_s ) {
            VectorVariable<2U> oldvector_variable5_;
            (*it)->PropertyValueAtSectorIntegrationPoint( key_node_VectorParameter, sector, ip_s, oldvector_variable5_ );
            cerr <<"  Old Method: nodal velocity at Sector IP No. "<< sector <<", is "<<oldvector_variable5_<<"\n" ;
          }
      
        cerr <<" \n" ;
      
        //++++++++++++++++++++++++++++++++++  Tensor
        // at the Facet IP with new scheme
        idy_ = 0 ;
        for ( auto lt : (*it)->AllSectorIntegrationPoints()  ){     //loop over Sector IP in the element
          TensorVariable<2U> tensor_variable5_ ;
          (lt).Obtain(key_node_TensorParameter , tensor_variable5_ );
          cerr <<"  New Method: nodal tensor at Sector IP "<< idy_ <<", is "<<tensor_variable5_(0,0)<<" ,"<<tensor_variable5_(1,0)<<" ,"<<tensor_variable5_(1,1)<<"\n" ;
          idy_++ ;
        }
        cerr <<" \n" ;
      
        // compare with old scheme at the Sector IP
        for ( size_t sector=0; sector< (*it)->Sectors(); ++sector )
          for ( size_t ip_s=0; ip_s < (*it)->IntegrationPointsPerSector(); ++ip_s ) {
            TensorVariable<2U> oldtensor_variable5_;
            (*it)->PropertyValueAtSectorIntegrationPoint( key_node_TensorParameter, sector, ip_s, oldtensor_variable5_ );
            cerr <<"  Old Method: nodal tensor at Sector IP No. "<< sector <<", is "<<oldtensor_variable5_(0,0)<<" ,"<<oldtensor_variable5_(1,0)<<" ,"<<oldtensor_variable5_(1,1)<<"\n" ;
          }
      
        cerr <<" \n" ;
      
      
        
    } //  end of loop over all elements

} // Run()



} // csmp
