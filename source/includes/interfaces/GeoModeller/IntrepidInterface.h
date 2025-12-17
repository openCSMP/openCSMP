#ifndef INTREPID_INTERFACE_H
#define INTREPID_INTERFACE_H

#include <iostream>
#include <cmath>
#include <cstdint>

namespace csmp {

template<uint32_t> class VSet;
class ModelTopology;

/// interface to import GeoModeller / CGAL tetra/triangle meshes; superseded by VTK connection
class IntrepidInterface {
	public:
 
    /// reads Intrepid meshfile into a VSet, identifying its regions in model topology
  	void Read( const char* filename, VSet<3U>&, ModelTopology& ) const;
 
		/// import model from Intrepid file format (see below for an example)
// useful ?!		Model<3U>* Read( const char* filename, const char* element_var) const;

  private:
    /// reverts node order of tetrahedral elements that do not comply with right-hand rule
    void RepairElementOrientations( VSet<3U>& ) const;
  
    /// writes regions file allowing user to select subregions that shall be part of computational model
    void OutputRegionsFile() const;
};


} // end csmp

#endif
