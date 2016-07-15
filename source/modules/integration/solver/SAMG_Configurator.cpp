#include "SAMG_Configurator.h"
#include <stdexcept>

namespace csmp {
  SAMG_Settings* SAMG_Configurator::GetSettingsByName( const std::string& name ) {
    if (name == "default") {
      return new SAMG_Settings();
    } else if (name == "two_unknowns") {
      SAMG_Settings* settings = new SAMG_Settings();
      settings->Set_napproach(2);
      settings->Set_nxtyp(2);
      settings->Set_ncgtyp(5);
      settings->Set_nred(1);
      return settings;
    } else {
      throw std::domain_error("Unknown name: " + name);
    }
      
  }
} // end namespace csp
