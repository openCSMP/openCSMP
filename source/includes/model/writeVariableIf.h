//
//  writeVariableIf.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 7/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_WRITE_VARIABLE_IF_H
#define CSMP_WRITE_VARIABLE_IF_H

#include "Index.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

namespace csmp {

/**
    Helper functions that checks complex boundary flags of a variable that shall be assigned.
    They functions only transfer values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on the Node, Element, Face etc.
    
    @author SKM 10/9/2014
*/
template<uint32_t dim, template<uint32_t> class CELL, class Var> requires CsmpVariable<dim, Var>
void writeVariableIf( CELL<dim>*, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
    std::cerr <<"\nwriteVariableIf( CELL<dim>*, const csmp::Index&, const Var&, VARIABLE_FLAG ): generic should never be called.\n";
 } // end generic specification


/// write guard for scalar variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      const csmp::Index& idx,
                      const ScalarVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(idx) != dont_overwrite )
      ptr->Store( idx, var );
 } // end version for scalars


/// write guard for scalar variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      const csmp::Index& idx,
                      const ArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(idx) != dont_overwrite )
      ptr->Store( idx, var );
 } // end version for arrays
 

/// write guard for vector variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      const csmp::Index& idx,
                      const VectorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( idx, vc );
    for ( auto i{0u}; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      const csmp::Index& idx,
                      const TensorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( idx, ts );
    for ( uint32_t i{0u}; i<dim; i++ )
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( uint32_t j{0U}; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( idx, ts );
 } // end version for tensors



/// write guard for vector variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      const csmp::Index& idx,
                      const FlaggedArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    FlaggedArrayVariable fa( var.Size() );
    ptr->Read( idx, fa );
    for ( auto i{0u}; i<var.Size(); i++ )
      // the component gets overwritten
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           fa.Flag(i) = var.Flag(i);
           fa(i)      = var[i];
        }
    ptr->Store( idx, fa );
 } // end version for vector variables



/**
    Helper functions that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on Element/Face/Interface integration points.
*/
template<uint32_t dim, template<uint32_t> class CELL, class Var> requires CsmpVariable<dim, Var>
void writeVariableIf( CELL<dim>*, uint32_t, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
    std::cerr <<"\nwriteVariableIf( CELL<dim>*, uint32_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG ): generic should never be called.\n";
 } // end generic specification

/// write guard for scalar variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const ScalarVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars
 
 
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const ArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for arrays
 
 

/// write guard for vector variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const VectorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( ip, idx, vc );
    for ( auto i{0u}; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const TensorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( ip, idx, ts );
    for ( uint32_t i{0u}; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( uint32_t j{0U}; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( ip, idx, ts );
 } // end version for tensors


/// write guard for vector variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const FlaggedArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    FlaggedArrayVariable fa( var.Size() );
    ptr->Read( ip, idx, fa );
    for ( uint32_t i{0u}; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           fa.Flag(i) = var.Flag(i);
           fa(i)      = var[i];
        }
    ptr->Store( ip, idx, fa );
 } 



/**
    Helper functions that check a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for finite volume-related integration points.
*/
template<uint32_t dim, template<uint32_t> class CELL, class Var> requires CsmpVariable<dim, Var>
void writeVariableIf( CELL<dim>*, uint32_t /* sector_or_facet */,
                      uint32_t, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
    std::cerr <<"\nwriteVariableIf( CELL<dim>*, uint32_t sector_or_facet, uint32_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG ): should never be called; only specialisation thereof.\n";
 } // end generic specification


/// write guard for scalar variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const ScalarVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars



template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const ArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for ArrayVariable
 


/// write guard for vector variables
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t sector_or_facet,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const VectorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( sector_or_facet, ip, idx, vc );
    for ( auto i{0u}; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( sector_or_facet, ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t sector_or_facet,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const TensorVariable<dim>& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( sector_or_facet, ip, idx, ts );
    for ( auto i{0u}; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( uint32_t j{0U}; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( sector_or_facet, ip, idx, ts );
    
 } // end version for tensors


template<uint32_t dim, template<uint32_t> class CELL>
void writeVariableIf( CELL<dim>* ptr,
                      uint32_t sector_or_facet,
                      uint32_t ip,
                      const csmp::Index& idx,
                      const FlaggedArrayVariable& var,
                      VARIABLE_FLAG dont_overwrite )
 {
    FlaggedArrayVariable fa( var.Size() );
    ptr->Read( sector_or_facet, ip, idx, fa );
    for ( auto i{0u}; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           fa.Flag(i) = var.Flag(i);
           fa(i)      = var[i];
        }
    ptr->Store( sector_or_facet, ip, idx, fa );
 } 




} // end csmp

#endif /* CSMP_WRITE_VARIABLE_IF_H */
