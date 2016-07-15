#include "Fluidproperties.h"

using namespace std;

namespace csmp
{
  Fluidproperties::Fluidproperties()
  {
    InitToZero();
  }

  Fluidproperties::~Fluidproperties()
  {
  }

  Fluidproperties::Fluidproperties(const Fluidproperties& fp)
  {
    *this = fp;
  }

  Fluidproperties& Fluidproperties::operator=(const Fluidproperties& fp)
  {
    if (&fp != this) 
      {
        t     = fp.t; 
        p     = fp.p;
        x     = fp.x;
        wt    = fp.wt;
        smf   = fp.smf;
        rho   = fp.rho;
        h     = fp.h;
        cp    = fp.cp;
        beta  = fp.beta;
        s     = fp.s;
        mf    = fp.mf;
        mu    = fp.mu;
        state = fp.state;
      }
    return *this;
  }

}//csmp
