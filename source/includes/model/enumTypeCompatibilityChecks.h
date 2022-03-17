//
//  enumTypeCompatibilityChecks.h
//
//  CSMP API library
//
//  Allows the compiler to test whether the types used in the creation or access of CSMP variables are 
//  actually compatible with the types of the objects that they are applied to, for instance,
//
//  SCALAR -> ScalarVariable   or   NODE -> Node<dim>
//
//  These methods are used by static_assert(), to identify incompatibilities already during compilation.
//
//  Created by Stephan Matthai on 20/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_ENUM_TYPE_COMPATIBILITY_CHECKS_H
#define CSMP_ENUM_TYPE_COMPATIBILITY_CHECKS_H

#include "CSMP_global_enumerations.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Region;
template<uint32_t> class Boundary;
template<uint32_t> class SplitBoundary;
template<uint32_t> class Model;

/** for compile-time type compatibility checking between csmp enumerations and types
        
        use like:
        
        static_assert( TypeMatchesVariablePlacement<type, placement>::value, "type does not match variable placement");
 */
template <template<uint32_t> class T, PLACEMENT IT> struct TypeMatchesVariablePlacement {
    enum { value = false };
 };
 
// full specialisations 
template <> struct TypeMatchesVariablePlacement<Node,NODE> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<Element,ELEMENT> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<Face,FACE> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<InterFace,INTER_FACE> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<Region,REGION> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<Boundary,BOUNDARY> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<SplitBoundary,SPLIT_BOUNDARY> {
    enum { value = true };
 };
template <> struct TypeMatchesVariablePlacement<Model,MODEL> {
    enum { value = true };
 };
 

/// For variables placed on the integration points (Gauss quadrature points) of Element, Face or InterFace objects
template <template<uint32_t> class T, PLACEMENT IT> struct TypeMatchesFEM_IP_Placement {
    enum { value = false };
 };
 
// full specialisations 
template <> struct TypeMatchesFEM_IP_Placement<Element,ELEMENT_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFEM_IP_Placement<Face,FACE_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFEM_IP_Placement<InterFace,INTER_FACE_INTEGRATION_POINT> {
    enum { value = true };
 };


/// For variables placed on finite volume sector integration points
template <template<uint32_t> class T, PLACEMENT IT> struct TypeMatchesFVM_SIP_Placement {
    enum { value = false };
 };
 
// full specialisations 
template <> struct TypeMatchesFVM_SIP_Placement<Element,SECTOR_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFVM_SIP_Placement<Face,FACE_SECTOR_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFVM_SIP_Placement<InterFace,INTER_FACE_SECTOR_INTEGRATION_POINT> {
    enum { value = true };
 };


/// For variables placed on finite volume facet integration points
template <template<uint32_t> class T, PLACEMENT IT> struct TypeMatchesFVM_FIP_Placement {
    enum { value = false };
 };
 
// full specialisations 
template <> struct TypeMatchesFVM_FIP_Placement<Element,FACET_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFVM_FIP_Placement<Face,FACE_FACET_INTEGRATION_POINT> {
    enum { value = true };
 };
template <> struct TypeMatchesFVM_FIP_Placement<InterFace,INTER_FACE_FACET_INTEGRATION_POINT> {
    enum { value = true };
 };


class ScalarVariable;
class ArrayVariable;
class FlaggArrayariable;
template<uint32_t> class VectorVariable;
template<uint32_t> class TensorVariable;

/** for compile-time type compatibility checking between csmp enumerations and types
        
        use like:
        
        static_assert( TypeMatchesVariableType<type, placement>::value, "type does not match variable placement");
 */
template <class T, VARIABLE_TYPE IT> struct TypeMatchesVariableType {
    enum { value = false };
 };

template <template<uint32_t> class T, VARIABLE_TYPE IT> struct TemplateTypeMatchesVariableType {
    enum { value = false };
 };
 
// full specialisations 
template <> struct TypeMatchesVariableType<ScalarVariable,SCALAR> {
    enum { value = true };
 };
template <> struct TemplateTypeMatchesVariableType<VectorVariable,VECTOR> {
    enum { value = true };
 };
template <> struct TemplateTypeMatchesVariableType<TensorVariable,TENSOR> {
    enum { value = true };
 };
template <> struct TypeMatchesVariableType<ArrayVariable,ARRAY> {
    enum { value = true };
 };
template <> struct TypeMatchesVariableType<FlaggedArrayVariable,FLAGGEDARRAY> {
    enum { value = true };
 };

 
 
 // more generic checks
 template <class T1, class T2> struct SameType {
      enum{value = false};
   };

template <class T> struct SameType<T,T> {
    enum { value = true };
  }; 

} // end csmp

#endif /* CSMP_ENUM_TYPE_COMPATIBILITY_CHECKS_H */
