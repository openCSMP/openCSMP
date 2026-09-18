// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SPLIT_REGION_DIFFUS_HEAT_H
#define SPLIT_REGION_DIFFUS_HEAT_H


#include "Model.h"
#include "SplitBoundary.h"

namespace csmp {

template<uint32_t> class Model;

/**
   @class SplitRegionDiffusiveHeating SplitRegionDiffusiveHeating.h

   CALCULATE explicit heat sources to account for thermal conduction in split boundary, old method, not used

   @author Benoit LC, ETH Zuerich
   @section contact Contact


   @changes changes Latest Changes

   @section motivation Motivation

   @section usage Usage

   @code
        
   @endcode
   
   @section dependencies Dependencies
   
   @section issues Known issues
   
   @section testing Testing
   testing was done in 2024

*/

template<uint32_t dim>
class SplitRegionDiffusiveHeating {

public:
    SplitRegionDiffusiveHeating(Model<dim>& model, const std::string& sb_name);

    virtual ~SplitRegionDiffusiveHeating();

    void Apply();

    const std::string sb_name;
    const std::string getName() const {return sb_name;}

private:

    void ComputeNodalArea();

    Model< dim >&                  model_ref_;
    const PropertyDatabase<dim>&   prop_ref_;

    SplitBoundary< dim >& sb;

    Index
        q_key,
        t_key,
        nA_key,
        kth_key,
        thick_key;

    ScalarVariable t_in, t_out, t_mid;
    ScalarVariable q_in, q_out, q_mid;
    ScalarVariable nA;
    ScalarVariable kth;
    ScalarVariable thick;

    Node<dim>* in_node;
    Node<dim>* out_node;
    Node<dim>* mid_node;

};
} // end namespace csmp


#endif
