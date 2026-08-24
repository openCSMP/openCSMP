#ifndef SOURCE_VISITOR_H
#define SOURCE_VISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"


namespace csmp {

struct Index;
template<uint32_t> class Model;


/**
@class SourceVisitor SourceVisitor.h

@author Thomas Driesner, ETH Zuerich
@c thomas.driesner@erdw.ethz.ch

*/
template<uint32_t dim>
class SourceVisitor final : public Visitor<dim>
{
public:
    SourceVisitor( Model<dim>&, std::vector<std::string>* to_initialize_keys = nullptr );
    SourceVisitor( Model<dim>&,
                   Index porosityKey,
                   Index densityDiffKey,
                   Index nfvsKey,
                   Index thicknessKey);
                   
    ~SourceVisitor() = default;
    
    void Visit( Element<dim>* ) override final;
    void Visit( Model<dim>* ) override final;

private:
    SourceVisitor();

    //! variables
    Model<dim>& model_;
    std::vector<std::string>* to_initialize_keys_;
    ScalarVariable phi, nfvs, timestep, thickness;
    csmp::Index porosityKey_, densityDiffKey_, nfvsKey_,timestepKey_,thicknessKey_;
};  

}// csmp


#endif //SOURCE_VISITOR_H

