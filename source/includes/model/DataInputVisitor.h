// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_DATA_INPUT_VISITOR_H
#define CSMP_DATA_INPUT_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class Region;
template<uint32_t> class Model;

/**

   Use this visitor to input raw data into model
   which is assumed to be ordered in input vector
   by variable components (for example vector v1x, v1y, v1z, v2x...)
   visitor will renumber elements or nodes.
   
   @author S. K. Matthai
   @date 2006

*/
template<typename Var, uint32_t dim>
class DataInputVisitor final : public Visitor<dim> {
  public:
    DataInputVisitor( Model<dim>& sg, 
                      const char* input_prop,
                      const std::vector<double>& input_data );
    
    void Visit( Model<dim>* ) override final;
    void Visit( Region<dim>* ) override final;
    void Visit( Element<dim>* ) override final;
    void Visit( Face<dim>* ) override final;
    void Visit( InterFace<dim>* ) override final;
    void Visit( Node<dim>* ) override final;
    
    void  Reset();

  private:
    const std::vector<double>&  input_data_ref_;
    csmp::Index             prop_key_;
    Var                     variable_;
    size_t                  counter_;
};

} // csmp

#endif







