// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_COPY_REPLACE_VISITOR_H
#define CSMP_COPY_REPLACE_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "AP_BoolVector.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Region;



/** Replaces distributed properties with the user supplied values

    @author S.K. Matthai
    @date 2006
    
    @note updated 17/12/2016 removing flag tracking
    and OpenMP stuff, new visitor should be built for
    this. Node tracking relied on indices which
    is no longer a give.
*/
template<typename Var, uint32_t dim>
class CopyReplaceVisitor final : public Visitor<dim> {
  public:
    CopyReplaceVisitor( const PropertyDatabase<dim>&, 
                        const char* prop_a,
                        const char* prop_b );
    
    void Visit( Model<dim>* ) override final;
    void Visit( Region<dim>* ) override final;
    void Visit( Element<dim>* ) override final;
    void Visit( Face<dim>* ) override final;
    void Visit( InterFace<dim>* ) override final;
    void Visit( Node<dim>* ) override final;

  private:
    csmp::Index  prop_key_a_, prop_key_b_;
    Var          variable_;
};


} // csmp


#endif







