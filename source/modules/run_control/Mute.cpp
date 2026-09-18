// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Mute.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 7/08/2016.
//

#include "Mute.h"

namespace csmp {


Mute::Mute()
    : cachebf_(nullptr), fileRedirect_(nullptr), nullOut_( new std::ofstream("/dev/null") )
  {
  }


  /// Dtor does not reinstate ostream
Mute::~Mute()
  {
  }


  /// Redirects to null stream, returns false if stream already redirected, true if redirection successful
  bool Mute::Null( std::ostream& os )
  {
    if(cachebf_)
      return false;

    cachebf_ = os.rdbuf();
    os.rdbuf( nullOut_->rdbuf() );
    return true;
  }


  /// Redirects to file, returns false if stream already redirected, true if redirection successful
  bool Mute::File( std::ostream& os, const char* fileName )
  {
    if( cachebf_ || fileRedirect_ )
      return false;

    fileRedirect_ = new std::ofstream;
    fileRedirect_->open(fileName);
    if( !fileRedirect_->is_open() )
      return false;

    cachebf_ = os.rdbuf();
    os.rdbuf( fileRedirect_->rdbuf() );
    return true;
  }


  /// Reinstates screen output, returns false if no redirection current, true if reinstation successful
bool Mute::Reinstate( std::ostream& os )
  {
    if(!cachebf_)
      return false;

    if(fileRedirect_)
      fileRedirect_->close();

    os.rdbuf(cachebf_);
    cachebf_ = nullptr;
    fileRedirect_ = nullptr;
    return true;
  }


  bool Mute::Redirected() const
  {
    return static_cast<bool>(cachebf_);
  }

}
