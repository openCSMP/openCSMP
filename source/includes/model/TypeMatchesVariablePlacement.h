//
//  TypeMatchesVariablePlacement.h
//  CSMP_FECFVM_Simulator
//
//  Created by Stephan Matthai on 20/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_TYPE_MATCHES_VARIABLE_PLACEMENT_H
#define CSMP_TYPE_MATCHES_VARIABLE_PLACEMENT_H

#include "CSMP_global_enumerations.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Region;
template<size_t> class Boundary;
template<size_t> class SplitBoundary;
template<size_t> class Model;

/** for compile-time type compatibility checking between csmp enumerations and types
        
        use like:
        
        static_assert( TypeMatchesVariablePlacement<type, placement>::value, "type does not match variable placement");
 */
template <template<size_t> class T, PLACEMENT IT> struct TypeMatchesVariablePlacement {
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
template <template<size_t> class T, PLACEMENT IT> struct TypeMatchesFEM_IP_Placement {
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
template <template<size_t> class T, PLACEMENT IT> struct TypeMatchesFVM_SIP_Placement {
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
template <template<size_t> class T, PLACEMENT IT> struct TypeMatchesFVM_FIP_Placement {
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
template<size_t> class VectorVariable;
template<size_t> class TensorVariable;

/** for compile-time type compatibility checking between csmp enumerations and types
        
        use like:
        
        static_assert( TypeMatchesVariableType<type, placement>::value, "type does not match variable placement");
 */
template <class T, VARIABLE_TYPE IT> struct TypeMatchesVariableType {
    enum { value = false };
 };

template <template<size_t> class T, VARIABLE_TYPE IT> struct TemplateTypeMatchesVariableType {
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

#endif /* CSMP_TYPE_MATCHES_VARIABLE_PLACEMENT_H */
