#ifndef AP_TEST_UTILITIES_H
#define AP_TEST_UTILITIES_H

// utilities of A. Paluszny created for testing the FV machinery
#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Element;
template<size_t,template<uint32_t> class> class FiniteVolumeTraits;

void rhinoOutput( const Model<3U>& );

void printElementStencil( const Element<3U>&, const std::string& );

void printElement( const Element<3U>&, const std::string& );

///creates a layer, sets a color, and sets it as the current layer
void createLayer ( std::stringstream & ss, const std::string & layer, size_t id_element );

void printFiniteVolumes( const Model<3U>& );

void writeSurfaceFacet( std::stringstream & ss, const FiniteVolumeTraits<3,Element>& efvt, 
                        const Element<3U>& e, size_t iFacet );
 
void writeNormal ( std::stringstream & ss, const FiniteVolumeTraits<3,Element>& efvt, 
                   const Element<3U>& e, size_t iFacet, bool bInverted );
  
} // end csmp

#endif 
