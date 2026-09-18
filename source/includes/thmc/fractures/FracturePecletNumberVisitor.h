// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FRACTUREPECLETNUMBERVISITOR_H
#define FRACTUREPECLETNUMBERVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"

namespace csmp{

  struct Index;
  template<uint32_t> class Model;


template<uint32_t dim>
class FracturePecletNumberVisitor : public Visitor<dim>
  {
  public:
    FracturePecletNumberVisitor( Model<dim>& model, const char* velocityTag,
                                 const char* pecletNumberTag, double nodalSourceSink);
    virtual ~FracturePecletNumberVisitor() {}

    virtual void Visit(Element<dim>* element);

  private:
    Model<dim>& model_;
    Index velKey_, pecletNumberKey_;
    const double sourceSink_;
    ScalarVariable FracNPe_;
  };
} //csmp

#endif // FRACTUREPECLETNUMBERVISITOR_H
