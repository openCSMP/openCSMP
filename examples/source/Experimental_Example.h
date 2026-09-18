// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Experimental_Example.h
//  CSMP_API_library2014
//
//  Created by Stephan Matthai on 1/27/14.
//

#ifndef EXPERIMENTAL_EXAMPLE_H
#define EXPERIMENTAL_EXAMPLE_H

#include "Example.h"

namespace csmp {

template<uint32_t> class Model;


class  Experimental_Example : public Example {
public:
  virtual void Run();
  virtual void Specifications();
};


/*
class  Experimental_Example : public Example {
public:
    virtual void Run();
    virtual void Specifications();

protected:
    bool ImportModelAndRunChecks( const std::string& model_name, const std::string& variables_file ) const;
    bool RunChecks( const Model<3U>& model ) const;

private:
    void RunSimulation( Model<3U>& model ) const;

    template<class FLOW_FUNCTIONS>
    void Compute2PhaseFlowProperties(
        Model<3U>& mdl, FLOW_FUNCTIONS& flowfunctions,
        bool with_gravity, bool with_tensor_k
    ) const;

    void ComputeSteadyStatePressure(
        Model<3U>& mdl, bool with_gravity, bool with_tensor_k
    ) const;
};
*/


} // csmp

#endif // EXPERIMENTAL_EXAMPLE_H
