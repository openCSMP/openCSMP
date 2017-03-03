#ifndef GRAVITY_PROJECTION_VISITOR_H
#define GRAVITY_PROJECTION_VISITOR_H

#include "Visitor.h"

#include <iostream>
#include "CSMP_definitions.h"

#include "Model.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"

namespace csmp{

template<size_t dim>
class GravityProjectionVisitor : public Visitor<dim>
{

  public:

    GravityProjectionVisitor(  Model<dim>& model,
                               const Index& prop_idx,
                               const Index& result_idx,
                               VectorVariable<dim> vec );


    virtual ~GravityProjectionVisitor();

    virtual void Visit(Element<dim>* element);
    virtual void Visit(Model<dim>* model);
    virtual void Visit(Region<dim>* model);

    /// Get Result
    Index Get_PropertyIndex();
    Index Get_ResultIndex();
    void  Get_Result( Element<dim>* eptr, VectorVariable<dim>& result );

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
