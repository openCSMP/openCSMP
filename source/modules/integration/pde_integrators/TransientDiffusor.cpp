#include "TransientDiffusor.h"
#include "Model.h"
#include "Exception.h"
#include "CSMP_highLevelUtilities.h"
#include "ModelTime.h"

using namespace std;

namespace csmp {

template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
TransientDiffusor<dim,COMPUTATION_DOMAIN>::TransientDiffusor( Model<dim>& sg,
                                                              const char* diffusivity,
                                                              const char* diffusing_variable,
                                                              const char* storage_variable,
                                                              const char* element_source_variable )
  
 :
#ifdef CSMP_WITH_SAMG_SOLVER
   settings_(),
   PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
#else
   /// add extra functionality for alternative solver if needed
   PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
#endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), element_source_variable, diffusing_variable) ),
   capacitance_lhs_( sg.Database(), storage_variable, diffusing_variable, diffusing_variable ),
   capacitance_rhs_( sg.Database(), storage_variable, diffusing_variable ),
   nodal_source_(0),
   gravity_(0),
   grad_multiplier_(1.),
   dep_var_name_(diffusing_variable)
 {
   	if ( !isoparametricElementMesh( sg ) )
		  throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
		                                     "elements are not isoparametric; use other Algorithm" ); 

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  diffusivity_key = p_ref.StorageKey(diffusivity);
    if ( diffusivity_key.place != ELEMENT  and  diffusivity_key.place != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusivity, "variable must be placed on element or constraint point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  storage_variable_key = p_ref.StorageKey(storage_variable);
    if ( (storage_variable_key.place != ELEMENT  and  storage_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         storage_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             storage_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(element_source_variable);
    if ( (source_variable_key.place != ELEMENT  and  source_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             element_source_variable, "variable must be a scalar placed on element or constraint point" ); 

    // setting the PDE operators up for transient computation
    capacitance_lhs_.MultiplyWithTimeIncrement(true);
    capacitance_lhs_.LumpedFormulation(true);
    capacitance_rhs_.MultiplyWithTimeIncrement(true);
    source_->AddAccumulateLater();
    
    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add( &capacitance_lhs_ );
    this->Add( &capacitance_rhs_ );
     
    // basic SAMG solver settings for hybrid element meshes
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( dim != 1U )  AdjustSolverSettings(); 
#endif

} // end constructor






template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
TransientDiffusor<dim,COMPUTATION_DOMAIN>::TransientDiffusor( Model<dim>& sg,
                                                                            const char* diffusivity,
                                                                            const char* diffusing_variable,
                                                                            const char* storage_variable,
                                                                            const char* element_source_variable,
                                                                            const char* point_source_variable )
                                     
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), element_source_variable, diffusing_variable) ),
   capacitance_lhs_( sg.Database(), storage_variable, diffusing_variable, diffusing_variable ),
   capacitance_rhs_( sg.Database(), storage_variable, diffusing_variable ),
   nodal_source_(new PointSource_rhsop<dim,ComputationCell>(sg.Database(), point_source_variable, diffusing_variable) ),
   gravity_(0),
   grad_multiplier_(1.),
   dep_var_name_(diffusing_variable)
 {
    if ( !isoparametricElementMesh( sg ) )
          throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):",
                                             "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  diffusivity_key = p_ref.StorageKey(diffusivity);
    if ( diffusivity_key.place != ELEMENT  and  diffusivity_key.place != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusivity, "variable must be placed on element or constraint point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  storage_variable_key = p_ref.StorageKey(storage_variable);
    if ( (storage_variable_key.place != ELEMENT  and  storage_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         storage_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             storage_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(element_source_variable);
    if ( (source_variable_key.place != ELEMENT  and  source_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             element_source_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  nsource_variable_key = p_ref.StorageKey(point_source_variable);
    if ( nsource_variable_key.place != NODE  and  nsource_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             point_source_variable, "variable must be a scalar placed on the node" ); 

    // setting the PDE operators up for transient computation
    capacitance_lhs_.MultiplyWithTimeIncrement(true);
    capacitance_lhs_.LumpedFormulation(true);
    capacitance_rhs_.MultiplyWithTimeIncrement(true);
    source_->AddAccumulateLater();
    nodal_source_->AddAccumulateLater();
    
    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add( &capacitance_lhs_ );
    this->Add( &capacitance_rhs_ );
    this->Add(  source_ );
    this->Add(  nodal_source_ );
     
    // basic SAMG solver settings for hybrid element meshes
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( dim != 1U )  AdjustSolverSettings(); 
#endif

} // end constructor (including point sources)




template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
TransientDiffusor<dim,COMPUTATION_DOMAIN>::TransientDiffusor( Model<dim>& sg,
                                                              const char* diffusivity,
                                                              const char* diffusing_variable,
                                                              const char* storage_variable,
                                                              const char* spatial_source_variable,
                                                              const char* gradient_variable,
                                                              double64 gradient_multiplier )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable) ),
   capacitance_lhs_( sg.Database(), storage_variable, diffusing_variable, diffusing_variable ),
   capacitance_rhs_( sg.Database(), storage_variable, diffusing_variable ),
   nodal_source_(0),
   gravity_(new NumIntegral_dNT_op_dV<dim,ComputationCell>(sg.Database(), gradient_variable, diffusing_variable) ),
   grad_multiplier_(gradient_multiplier),
   dep_var_name_(diffusing_variable)
 {
    if ( !isoparametricElementMesh( sg ) )
          throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):",
                                              "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  diffusivity_key = p_ref.StorageKey(diffusivity);
    if ( diffusivity_key.place != ELEMENT  and  diffusivity_key.place != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusivity, "lefthand variable must be placed on element or constraint point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  storage_variable_key = p_ref.StorageKey(storage_variable);
    if ( (storage_variable_key.place != ELEMENT  and  storage_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         storage_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             storage_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( (source_variable_key.place != ELEMENT  and  source_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  gradient_variable_key = p_ref.StorageKey(gradient_variable);
    if ( (gradient_variable_key.place != ELEMENT || gradient_variable_key.place != ELEMENT_INTEGRATION_POINT) && gradient_variable_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<dim>::(constructor):", 
                             gradient_variable, 
                            "variable must be a vector placed on element or integration point" ); 

   // setting the PDE operators up for transient computation
    capacitance_lhs_.MultiplyWithTimeIncrement(true);
    capacitance_lhs_.LumpedFormulation(true);
    capacitance_rhs_.MultiplyWithTimeIncrement(true);
    source_->AddAccumulateLater();
    
    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add( &capacitance_lhs_ );
    this->Add( &capacitance_rhs_ );
    this->Add(  gravity_ );
    gravity_->AddAccumulateLater();
     
    // basic SAMG solver settings for hybrid element meshes
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( dim != 1U )  AdjustSolverSettings(); 
#endif

} // end constructor (gravity 1)




template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
TransientDiffusor<dim,COMPUTATION_DOMAIN>::TransientDiffusor( Model<dim>& sg,
                                                                            const char* lhs_diffusivity,
                                                                            const char* rhs_diffusivity,
                                                                            const char* diffusing_variable,
                                                                            const char* storage_variable,
                                                                            const char* spatial_source_variable,
                                                                            const char* gradient_variable, 
                                                                            double64 gradient_multiplier )
                                     
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
      settings_(),
      PDE_Integrator_CRM<dim,COMPUTATION_DOMAIN>( new SAMG_Solver(&settings_) ),
   #else
      /// add extra functionality for alternative solver if needed
      PDE_Integrator<dim,COMPUTATION_DOMAIN>( new CSMP_DEFAULT_LINEAR_SOLVER() ),
   #endif
   conductance_( sg.Database(), lhs_diffusivity, diffusing_variable, diffusing_variable ),
   source_(new NumIntegral_NT_op_N_dV<dim,ComputationCell>(sg.Database(), spatial_source_variable, diffusing_variable) ),
   capacitance_lhs_( sg.Database(), storage_variable, diffusing_variable, diffusing_variable ),
   capacitance_rhs_( sg.Database(), storage_variable, diffusing_variable ),
   nodal_source_(0),
   gravity_(new NumIntegral_dNT_op_dV<dim,ComputationCell>(sg.Database(), gradient_variable, diffusing_variable) ),
   grad_multiplier_(gradient_multiplier),
   dep_var_name_(diffusing_variable)
 {
    if ( !isoparametricElementMesh( sg ) )
          throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):",
                                              "elements are not isoparametric; use other Algorithm" );

    const PropertyDatabase<dim>&   p_ref = sg.Database();
    
    csmp::Index  lhs_diffusivity_key = p_ref.StorageKey(lhs_diffusivity);
    if ( lhs_diffusivity_key.place != ELEMENT  and  lhs_diffusivity_key.place != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             lhs_diffusivity, "lefthand variable must be placed on element or constraint point" ); 

    csmp::Index  rhs_diffusivity_key = p_ref.StorageKey(rhs_diffusivity);
    if ( rhs_diffusivity_key.place != ELEMENT  and  rhs_diffusivity_key.place != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             rhs_diffusivity, "righthand variable must be placed on element or constraint point" ); 

    csmp::Index  diffusing_variable_key = p_ref.StorageKey(diffusing_variable);
    if ( diffusing_variable_key.place != NODE || diffusing_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             diffusing_variable, "variable must be a scalar placed on the node" ); 

    csmp::Index  storage_variable_key = p_ref.StorageKey(storage_variable);
    if ( (storage_variable_key.place != ELEMENT  and  storage_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         storage_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             storage_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  source_variable_key = p_ref.StorageKey(spatial_source_variable);
    if ( (source_variable_key.place != ELEMENT  and  source_variable_key.place != ELEMENT_INTEGRATION_POINT) ||
         source_variable_key.type != SCALAR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<>::(constructor):", 
                             spatial_source_variable, "variable must be a scalar placed on element or constraint point" ); 

    csmp::Index  gradient_variable_key = p_ref.StorageKey(gradient_variable);
    if ( (gradient_variable_key.place != ELEMENT || gradient_variable_key.place != ELEMENT_INTEGRATION_POINT) && gradient_variable_key.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "TransientDiffusor<dim>::(constructor):", 
                             gradient_variable, 
                            "variable must be a vector placed on element or integration point" ); 

   // setting the PDE operators up for transient computation
    capacitance_lhs_.MultiplyWithTimeIncrement(true);
    capacitance_lhs_.LumpedFormulation(true);
    capacitance_rhs_.MultiplyWithTimeIncrement(true);
    source_->AddAccumulateLater();
    
    // initializing the algorithm  
    this->Add( &conductance_ );
    this->Add(  source_ );
    this->Add( &capacitance_lhs_ );
    this->Add( &capacitance_rhs_ );
    this->Add(  gravity_ );
    gravity_->AddAccumulateLater();
     
    // basic SAMG solver settings for hybrid element meshes
#ifdef CSMP_WITH_SAMG_SOLVER
    if ( dim != 1U )  AdjustSolverSettings(); 
#endif

} // end constructor (gravity)




template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
TransientDiffusor<dim,COMPUTATION_DOMAIN>::~TransientDiffusor()
 {
    delete source_;
    delete nodal_source_;
    delete gravity_;
 } 







template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void TransientDiffusor<dim,COMPUTATION_DOMAIN>::ComputeTransientStateFullyImplicit(
                                                                          Model<dim>& model,
                                                                          double64 time_increment,
                                                                          bool verbose )
 {
    double64& model_time( ModelTime::Instance().modelTime );

    cout <<"\n\n\n\nTransientDiffusor<"<< dim;
    cout <<">: ComputeTransientStateFullyImplicit: ";
    cout <<"Computing transient '"<< dep_var_name_ <<"' at t="<< model_time;
    cout <<", after time increment: "<< time_increment <<" s."<< endl;
    this->TimeIncrement( 1. / time_increment );
    model.Apply ( *this );

 } // ComputeTransientStateFullyImplicit




template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
void TransientDiffusor<dim,COMPUTATION_DOMAIN>::AdjustSolverSettings()
 {
#ifdef CSMP_WITH_SAMG_SOLVER
    // basic solver settings for hybrid element meshes
    settings_.Set_ncgtyp(5); 
    settings_.Set_nxtyp(1); // 1=ILU, 0=Gauss Seidel
    settings_.Set_ndefault(40);
    settings_.Set_idmp(1); // suppress output
#else
    /// add extra functionality for alternative solver if needed
#endif
 } // end AdjustSolverSettings


#ifdef CSMP_WITH_SAMG_SOLVER
    template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
    SAMG_Settings& TransientDiffusor<dim,COMPUTATION_DOMAIN>::GetSolverSettings()
     { return settings_; }
#else
    /// add extra functionality for alternative solver if needed
#endif


template class TransientDiffusor<1U,Region>;
template class TransientDiffusor<2U,Region>;
template class TransientDiffusor<3U,Region>;

template class TransientDiffusor<1U,Boundary>;
template class TransientDiffusor<2U,Boundary>;
template class TransientDiffusor<3U,Boundary>;

} // end csmp
