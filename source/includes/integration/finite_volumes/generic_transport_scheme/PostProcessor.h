// Copyright © 2020 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_POST_PROCESSOR_H
#define CSMP_POST_PROCESSOR_H

#include "IntegralEquation.h"

namespace csmp {

/**
    Policy of the Implicit Transport for the post-processing of the computation.
*/
template<uint32_t dim, template<uint32_t> class USER>
class PostProcessor {
  public:
    void PostProcess();

  private:
    /// shorthand for accessing the class that this is a policy of
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

} // end csmp

#endif /* CSMP_POST_PROCESSOR_H */

// template<uint32_t dim,template<uint32_t> class CELL> ModelSubDomain

