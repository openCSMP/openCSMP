// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_TWO_PHASE_MODEL_TEST_H
#define CSMP_TWO_PHASE_MODEL_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class TwoPhaseModel;

/// PL Nov 2010
class TwoPhaseModel_Test : public Test {
public:
  explicit TwoPhaseModel_Test( Model<1U>*,
                               TwoPhaseModel<1U>*,
                               const char* testName,
                               const char* krn,
                               const char* krw,
                               const char* pc,
                               size_t upToElementNumber = 0 );
  ~TwoPhaseModel_Test();
  virtual void run();

private:

   Model<1U>* model_;

   TwoPhaseModel<1U>* twoPhaseModel_;

   long upToElementNumber_;

   Index krn_key_,
         krw_key_,
         pc_key_;

   std::list<std::string> vtuOutputProps_;
};

} // csmp

#endif // TWO_PHASE_MODEL_TEST_H
