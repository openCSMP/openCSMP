// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRAVITY_PROJECTION_VISITOR_H
#define GRAVITY_PROJECTION_VISITOR_H

#include "Visitor.h"
#include "VectorVariable.h"

namespace csmp {

template<uint32_t> class Model;

/**
       TODO: Who? - What for?
*/
template<uint32_t dim>
class GravityProjectionVisitor final : public Visitor<dim>
  {
    public:
      GravityProjectionVisitor(  Model<dim>&,
                                 const Index& prop_idx,
                                 const Index& result_idx,
                                 VectorVariable<dim>& );

      void Visit(Element<dim>* ) override final;
      void Visit(Model<dim>* ) override final;
      void Visit(Region<dim>* ) override final;

      /// Get Result
      Index Get_PropertyIndex();
      Index Get_ResultIndex();
      void  Get_Result( Element<dim>*, VectorVariable<dim>& result );

    private:

      /// Property indices
      Index               prop_idx_;
      Index               result_idx_;

      /// Temporary data
      VectorVariable<dim> vec_;
      VectorVariable<dim> proj_;

      /// Flags
      bool                const_vec_;
      bool                debug_;
  };

} //csmp

#endif // GravityProjectionVisitor
