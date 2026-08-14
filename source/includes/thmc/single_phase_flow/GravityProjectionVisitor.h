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
class GravityProjectionVisitor : public Visitor<dim>
  {
    public:

      GravityProjectionVisitor(  Model<dim>&,
                                 const Index& prop_idx,
                                 const Index& result_idx,
                                 VectorVariable<dim>& );

      virtual ~GravityProjectionVisitor();

      virtual void Visit(Element<dim>* );
      virtual void Visit(Model<dim>* );
      virtual void Visit(Region<dim>* );

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
