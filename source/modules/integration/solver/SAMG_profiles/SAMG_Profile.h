// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SAMG_PROFILE_H
#define SAMG_PROFILE_H

#include <string>

namespace csmp{

class SAMG_Settings;

class SAMG_Profile {
  public:
    SAMG_Profile();
    virtual ~SAMG_Profile();

    /// Solve for matrix equation
    virtual bool Solve( double modelTime ) = 0;

    std::string   FileDumpName() const { return dumpFileName_; }
    void          FileDumpName( const std::string& dumpFileName ) { dumpFileName_ = dumpFileName; }

protected:
    std::string   dumpFileName_;

};

} // csmp

#endif // SAMG_PROFILE_H

