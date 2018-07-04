#ifndef SOURCE_VISITOR_H
#define SOURCE_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"


namespace csmp {

struct Index;
template<size_t> class Model;

template<size_t dim> 
class SourceVisitor : public Visitor<dim>
{
public:
    SourceVisitor( Model<dim>& model , std::vector<std::string> *to_initialize_keys =NULL);
    SourceVisitor(Model<dim>& model,
                   Index porosityKey,
                   Index densityDiffKey,
                   Index nfvsKey,
                   Index thicknessKey);
    ~SourceVisitor();
    virtual void Visit(Element<dim>* e);
    virtual void Visit(Model<dim>* model);

private:
    SourceVisitor();

    //! variables
    Model<dim>& model_;
    std::vector<std::string>* to_initialize_keys_;
    ScalarVariable phi, nfvs, timestep, thickness;
    csmp::Index porosityKey_, densityDiffKey_, nfvsKey_,timestepKey_,thicknessKey_;
};  

}// csmp

/**
@class SourceVisitor SourceVisitor.h

@author Thomas Driesner, ETH Zuerich
@section contact Contact
thomas.driesner@erdw.ethz.ch

@section motivation Motivation

@section usage Usage

@code
@endcode

@section dependencies Dependencies

@section issues Known issues 

@section testing Testing
*/

#endif //SOURCE_VISITOR_H

