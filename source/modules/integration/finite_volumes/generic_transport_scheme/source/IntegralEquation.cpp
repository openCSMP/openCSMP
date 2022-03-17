//
//  IntegralEquation.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 22/7/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "IntegralEquation.h"
#include "ErrorHandler.h"
#include "PropertyDatabase.h"

// ADD OTHER POTENTIAL SETS HERE
#include "VariableSet_TracerTransfer.h"

#include "DiffusionLHS.h"
#include "PoreVolumeLHS.h"
#include "PoreVolumeRHS.h"
#include "UpwindFluxLHS.h"
#include "SourceTermRHS.h"
#include "FluidCompressibilitySourceLHS.h"
//#include "SourceTermLHS.h"

using namespace std;

namespace csmp {

/**
       Default constructor - advection with sources and sinks.

   upstream weighted first-order space-time fluxes
*/
template<uint32_t dim, class VARIABLE_SET>
IntegralEquation<dim,VARIABLE_SET>::IntegralEquation( const PropertyDatabase<dim>& pref )
 : VARIABLE_SET(pref),
   Notation(*this)
 {
   const bool mult_with_dt(true);
   
   mat_operators_.insert( make_pair( Operation("PoreVolumeLHS",csmp::ADD,1,mult_with_dt), 
                                     new PoreVolumeLHS<dim>(Notation.key_SPV,Notation.key_FVPV)) );
                                     
   mat_operators_.insert( make_pair( Operation("UpwindFluxLHS",csmp::ADD,1,false), 
                                     new UpwindFluxLHS<dim>(Notation.key_ff)) );
                                     
   vec_operators_.insert( make_pair( Operation("PoreVolumeRHS",csmp::ADD,1,mult_with_dt), 
                                     new PoreVolumeRHS<dim>(Notation.key_SPV,Notation.key_FVPV,Notation.key_C1)) );

   vec_operators_.insert( make_pair( Operation("SourceTermRHS",csmp::ADD,2,false), 
                                     new SourceTermRHS<dim>( Notation.key_QV,
                                                             Notation.key_NQV,
                                                             Notation.key_C1)) );
    Out();
    
 } // end constructor








/** 
     builds equation for one of the available options
        
          Basic advection diffusion equation in finite volume form:


 */
template<uint32_t dim, class VARIABLE_SET>
IntegralEquation<dim,VARIABLE_SET>::IntegralEquation( const PropertyDatabase<dim>& pref, 
                                                      const set<csmp::ADE_TERM>& ADE_components )
 : VARIABLE_SET(pref),
   Notation(*this)
 {
    ErrorHandler& error_handler( ErrorHandler::Instance() );

    csmp::INDEX<SCALAR,SECTOR_INTEGRATION_POINT> spv_key;
    csmp::INDEX<SCALAR,NODE>                     fpv_key, adv_key;
    csmp::INDEX<SCALAR,FACET_INTEGRATION_POINT>  ff_key;
    
    for ( set<csmp::ADE_TERM>::const_iterator
          it=ADE_components.begin(); it!=ADE_components.end(); ++it )
      {
         // TODO: get rid of this super ugly switch statement
         switch( *it ) {
            case ADVECTION: {
                   const bool mult_with_dt(true);
                   // upstream weighted first-order space-time fluxes
                   mat_operators_.insert( make_pair( Operation("PoreVolumeLHS",csmp::ADD,1,mult_with_dt), 
                                                     new PoreVolumeLHS<dim>(Notation.key_SPV,Notation.key_FVPV)) );
                                                     
                   mat_operators_.insert( make_pair( Operation("UpwindFluxLHS",csmp::ADD,1,false), 
                                                     new UpwindFluxLHS<dim>(Notation.key_ff)) );
                                                     
                   vec_operators_.insert( make_pair( Operation("PoreVolumeRHS",csmp::ADD,1,mult_with_dt), 
                                                     new PoreVolumeRHS<dim>(Notation.key_SPV,Notation.key_FVPV,Notation.key_C1)) );
                }
              break;
            case FE_DIFFUSION:
                 mat_operators_.insert( make_pair( Operation("DiffusionLHS",csmp::ADD,1,false), 
                                                   new DiffusionLHS<dim>(Notation.key_D,Notation.key_C1)) );
              break;
            case DISPERSION:
                 error_handler.notice( ERROR, "IntegralEquation<dim,VARIABLE_SET>::IntegralEquation", "DISPERSION not implemented yet");
              break;
            case SOURCE_SINK:
                 vec_operators_.insert( make_pair( Operation("SourceTermRHS",csmp::ADD,2,false), 
                                                   new SourceTermRHS<dim>( Notation.key_QV,
                                                                           Notation.key_NQV,
                                                                           Notation.key_C1)) );
              break;
            case COMPRESSIBILITY:
                 mat_operators_.insert( make_pair( Operation("FluidCompressibilitySourceLHS",csmp::ADD,2,false), 
                                                   new FluidCompressibilitySourceLHS<dim>( Notation.key_SPV,
                                                                                           Notation.key_PHI,
                                                                                           Notation.key_CT,
                                                                                           Notation.key_PF0,
                                                                                           Notation.key_PF,
                                                                                           true /* interpolate_pf_to_sector_ip */ )) );
              break;
            case SORPTION: // essentially a sink term
                 error_handler.notice( ERROR, "IntegralEquation<dim,VARIABLE_SET>::IntegralEquation", "SORPTION not implemented yet");
              break;
            case DECAY: // and another sink term
                 error_handler.notice( ERROR, "IntegralEquation<dim,VARIABLE_SET>::IntegralEquation", "DECAY not implemented yet");
              break;
            default:
                error_handler.notice( ERROR, "IntegralEquation:ctor", "ADE component term not recognised." );
              return;
         }
      }

    Out();
    
 } // end constructor
 
 
 
/**
   delete all dynamic copies of  MathOperators from the heap
*/
template<uint32_t dim, class VARIABLE_SET>
IntegralEquation<dim,VARIABLE_SET>::~IntegralEquation()
 {
    for ( typename map<Operation,MatrixOperator<dim>*>::iterator
          it=mat_operators_.begin(); it!=mat_operators_.end(); ++it )
      delete (*it).second;   

    for ( typename map<Operation,VectorOperator<dim>*>::iterator
          it=vec_operators_.begin(); it!=vec_operators_.end(); ++it )
      delete (*it).second;   
          
    for ( typename map<string,MatrixOperator<dim>*>::iterator
          it=postpro_operators_.begin(); it!=postpro_operators_.end(); ++it )
      delete (*it).second;   
 }
    
    
    
/**
    operators are added in the user-specified sequence (add, subtract, multiply, specified in operators 

\f{ \mathbf{u}_{\alpha} = -\frac{k_{r\alpha}}{\mu_\alpha} \mathbb{K}\left( \nabla p_\alpha - \rho_\alpha \mathbf{g}\right)\f}

*/
template<uint32_t dim, class VARIABLE_SET>
void IntegralEquation<dim,VARIABLE_SET>::Add( const Operation& operation, MatrixOperator<dim>* mop ) // makes internal copies of these operators
 {
    ErrorHandler& error_handler( ErrorHandler::Instance() );

    // checking for incompatibilities
    if ( operation.ExecutionLevel() == 1 && 
         ( operation.OperationType() == MULTIPLY || operation.OperationType() == DIVIDE ) )
      error_handler.notice( ERROR, "IntegralEquation<dim,VARIABLE_SET>::Add", "mutiply/divide have no effect on empty matrix." );
    
    // adding lefthand operator
    pair<typename map<Operation,MatrixOperator<dim>*>::iterator,bool> it = mat_operators_.insert( make_pair(operation,mop) );
    if ( it.second == false ) {
         operation.Out();
         error_handler.notice( FATAL_ERROR, "IntegralEquation<dim,VARIABLE_SET>::Add", "unable to add operation & matrix (lefthand) operator");
      }
 }
 
 
template<uint32_t dim, class VARIABLE_SET>
void IntegralEquation<dim,VARIABLE_SET>::Add( const Operation& operation, VectorOperator<dim>* vop )
 {
    ErrorHandler& error_handler( ErrorHandler::Instance() );

    // checking for incompatibilities
    if ( operation.ExecutionLevel() == 1 && 
         ( operation.OperationType() == MULTIPLY || operation.OperationType() == DIVIDE ) )
      error_handler.notice( ERROR, "IntegralEquation<dim,VARIABLE_SET>::Add", "mutiply/divide have no effect on empty vector." );

    // adding righthand operator
    pair<typename map<Operation,VectorOperator<dim>*>::iterator,bool> it = vec_operators_.insert( make_pair(operation,vop) );
    if ( it.second == false ) {
         operation.Out();
         error_handler.notice( FATAL_ERROR, "IntegralEquation<dim,VARIABLE_SET>::Add", "unable to add operation & (righthand) vector operator");
      }
 }


/// post-processing operations involving the solution variable (these will never have to be modified)
template<uint32_t dim, class VARIABLE_SET>
void IntegralEquation<dim, VARIABLE_SET>::AddPostProcess( MatrixOperator<dim>* mop )
 {
    ErrorHandler& error_handler( ErrorHandler::Instance() );

    const string opname(typeid(mop).name());
    pair<typename map<string,MatrixOperator<dim>*>::iterator,bool> it = postpro_operators_.insert( make_pair(opname,mop) );
    if ( it.second == false ) {
         cerr <<"\n\t"<< opname;
         error_handler.notice( FATAL_ERROR, "IntegralEquation<dim,VARIABLE_SET>::Add", "unable to add operation & (righthand) vector operator");
      }
 }
 


template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::MatrixOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::LHS_OperatorsBegin() const 
  { return mat_operators_.begin(); }

template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::MatrixOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::LHS_OperatorsEnd() const 
  { return mat_operators_.end(); }
    
template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::VectorOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::RHS_OperatorsBegin() const 
  { return vec_operators_.begin(); }

template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::VectorOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::RHS_OperatorsEnd() const 
  { return vec_operators_.end(); }
    
template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::PostProcessingOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::PostProcessingOperatorsBegin() const
 { return postpro_operators_.begin(); }
 
template<uint32_t dim, class VARIABLE_SET>
typename IntegralEquation<dim,VARIABLE_SET>::PostProcessingOperatorConstIterator IntegralEquation<dim,VARIABLE_SET>::PostProcessingOperatorsEnd() const
 { return postpro_operators_.end(); }
    



template<uint32_t dim, class VARIABLE_SET>
void IntegralEquation<dim,VARIABLE_SET>::Out() const
 {
    cout <<"\nIntegralEquation<dim>::Out: ";
    cout <<"\n\tMatrix operators: \n";
    for ( typename map<Operation,MatrixOperator<dim>*>::const_iterator
          it=mat_operators_.begin(); it!=mat_operators_.end(); ++it ) {
         (*it).first.Out();
         // (*it).second->Out();   
      }

    cout <<"\n\tVector operators: \n";
    for ( typename map<Operation,VectorOperator<dim>*>::const_iterator
          it=vec_operators_.begin(); it!=vec_operators_.end(); ++it ) {
         (*it).first.Out();
         // (*it).second->Out();   
      }
          
    cout <<"\n\tPost-processing operators: \n";
    for ( typename map<string,MatrixOperator<dim>*>::const_iterator
          it=postpro_operators_.begin(); it!=postpro_operators_.end(); ++it ) {
         cout <<"\n\t"<< (*it).first;
         // (*it).second->Out();   
      }
    
 } // end Out


template class IntegralEquation<3U,variables::VariableSet_TracerTransfer>;

} // end csmp
