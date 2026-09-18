// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  UG4_UGX_FileExport.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 3/12/2023.
//

#ifndef UG4_UGX_FILE_EXPORT_H
#define UG4_UGX_FILE_EXPORT_H

#include "CSMP_definitions.h"

namespace csmp {

/// RGBA (0-1) colours to distinguish 'subsets' in PROMESH
std::tuple<float,float,float,float,std::string> goodColour( int id );

template<uint32_t> class Model;
template<uint32_t> class Region;
template<uint32_t> class Boundary;
template<uint32_t> class MeshManager;

// TODO: write functions that output edge nodes as pointer vectors
// TODO: make sure that methods work for quadratic element meshes



/**

Export interface for writing text files in the .ugx format using XML for the formatting.
Explanations of the UG4 output tags.

Drawing on the information provided in 'file_io_ugx*' in conjunction with the GridWriterUGX class.
Some enumerations in grid/grid_base_objects.h

To assign properties to the grid, attachment traits are used.
These are: vertex, edge, face and volume.

@attention Design choices had to be made to crerate the additional edge and face data
that UG4 needs for during the output:

Importantly, CSMP models store edges and faces only if these are needed as lower-dimensional entities
for computations (well and fracture modelling, assignment of boundary condtions etc.). In UG4 they are always
needed and their numbering must be unique.

The choice made for the implementation of this CSMP-UG4 interface is that those edges and faces that
are  present in the CSMP model  are later identified in the edges and faces vectors created for UG.
What is output is the "model domain" plus the faces (and interfaces) touching it.
Their 'idx' is mapped that of the edge and faces created for UG.

Region and Boundary output enlists all the nodes, edges and faces that belong to these,
although UG4 only wants them to be listed once. Here, the feature in Promesh which eliminates
excess edges and faces is relied upon.

@attention properties are no associated with regions but communicated as node, edge, face and volume attachments.
Since the CSMP model only contains property values for those elements, faces, interfaces and nodes
present in there, no-data values are printed for all the other ones needed by UG.

@section Example

One of the UGX files from the examples directorty ofg UG4.
@Code

<?xml version="1.0" encoding="utf-8"?>
<grid name="defGrid">
	<vertices coords="3">0 0 0 0.5 ... </vertices>
  
	<edges>1 8  ...</edges>
	<triangles>1 8 232  ...</quadrilaterals>
	<tetrahedrons>90 36 226 232 226 ...</tetrahedrons>
	<hexahedrons>10 92 196 98 6 84 193 102 ...</hexahedrons>
	<pyramids>80 90 8 1 232 ...</pyramids>
  
  <volume_attachment name="PORO" type="double" passOn="1" global="1">0.34971 0.337518 0.326172  ...</volume_attachment>
	<volume_attachment name="PERM" type="double" passOn="1" global="1">1.91727 2.05433 2.14788 ...</volume_attachment>

	<subset_handler name="defSH">
		<subset name="Flex" color="0.588235 0.588235 1 1" state="262144">
			<vertices>12 13 14 15 16 17 ...</vertices>
			<edges>8 9 10 11 12 13 14 15 ...</edges>
			<faces>338 339 340 341 ...</faces>
			<volumes>80 81 82 83 ...</volumes>
		</subset>
    
		<subset name="Fixed" color="1 0 0 1" state="393216">
			<vertices>36 226 90 232 ...</vertices>
			<edges>536 73 582 539 501 ...</edges>
			<faces>84 10 80 85 86 81 11 ...</faces>
			<volumes>0 1 2 3 4 5 6 7 8 9 10 ...</volumes>
		</subset>
    
		<subset name="Force" color="0 1 0 1" state="393216">
			<vertices>63 182 224 144 145 225</vertices>
			<edges>256 456 454 181 182 459 461</edges>
			<faces>517 523</faces>
		</subset>

	</subset_handler>
	<selector name="defSel"/>
	<projection_handler name="defPH" subset_handler="0">
		<default type="default">0 0</default>
	</projection_handler>
</grid>

@endcode

@section XML tags used

The hierarchical XML tags above are getting processed here using the functionality in XML.
Comments in XML files are inserted using <!-- text -->
The tags suported here were defined by Andreas Vogel, and this is the correspondence list to CSMP types:

grid -> mesh
vertices -> nodes
edges, triangles, quadrilaterals -> line elements, triangular elements, quadrilateral elements etc.
subset = region or boundary
subset_handler -> class without specific purpose yet. Only a single subset handler needs to be present.
Its name is "defSH" (default subset handler)

projection_handler  for the transformation of coordinates to another systems

vertex, edge, face, volume -> element type independent identifier of elements of a specific topologic dimensionality

@attention the .ugx format requires lists for all the subentities of the actual finite elements in a mesh.
These are listed using the implicit integer numbers of the nodes ( vertices). It follows, that a 3D models needs:

<edges> unique pairs of nodes that define the edges of all elements in the model

<triangles> the vertex ids (0..n-1) that define the faces of simplex elements
<quadrilaterals> the faces of hexadral, prism, and pyramid elements

Note that these numbers have to be repeated for each region and boundary of the mesh, using the same numbering implicit to their order.


@section Storing the properties of a model

Variables are stored as "attachments."


@section Variable Types

common_attachments.h defines the types of variables that can be associated with nodes, edges, faces and elements. These are spelled out also in`
global_attachments.h -> attachment_info_traits -> distinguishing  double, float, int, bool, vector1, vector2, vector3, vector4.

Vectors and matrices are typedef'd  in  /common/math/ug_math_types.h

matrix22, matrix33, matrix44

All these types are referred to via strings in the output .ugx file.

@attention UGX file only stores straight-sided elements with corner nodes only. Any refinement occurs inside of UG.
This may include obtaining body-fitted representations of curves etc.

@author SKM
#date 3/12/2-23

*/
template<uint32_t dim>
class UG4_UGX_FileExport {
  public:
    explicit UG4_UGX_FileExport( const Model<dim>& );

    /// filename without extension which will be appended
    bool Write_UGX_FileASCII( const Model<dim>&, const std::string& file_name );

  private:
    UG4_UGX_FileExport() = delete;
    UG4_UGX_FileExport( const UG4_UGX_FileExport& ) = delete;
    
    /// Loops over the unique regions in the model to make a vector of  the unique edges (node-iD pairs) in the model (corner nodes only)
    size_t CollectEdges( const Model<dim>& );
    
    /// creates two sorted vectors for triangular and quadrilateral elements, respectively
    size_t CollectFaces( const Model<dim>& );
    
    /// creates sorted vectors for tetra, hexa, prism, and pyramid elements
    size_t CollectVolumes( const Model<dim>& );
    
    // SUBSET mapping
    /// Finds the ID numbers of the edges in the global edge vector inside of the region
    size_t CollectEdgesInRegion( const Region<dim>&, std::vector<size_t>& region_edges ) const;
    size_t CollectFacesInRegion( const Region<dim>&, std::vector<size_t>& region_faces ) const;
    size_t CollectVolumesInRegion( const Region<dim>&, std::vector<size_t>& region_volumes ) const;

    size_t CollectEdgesInBoundary( const Boundary<dim>&, std::vector<size_t>& region_edges ) const;
    size_t CollectFacesInBoundary( const Boundary<dim>&, std::vector<size_t>& region_faces ) const;
    
    // helpers
    bool HasTriangles() const      { return !tria_faces_.empty(); }
    bool HasQuadrilaterals() const { return !quad_faces_.empty(); }
    bool HasTetrahedra() const     { return !tetra_volumes_.empty(); }
    bool HasHexahedra() const      { return !hexa_volumes_.empty(); }
    bool HasPrisms() const         { return !prism_volumes_.empty(); }
    bool HasPyramids() const       { return !pyra_volumes_.empty(); }
    
    /// returns CSMP line element/face/interface  idx if it exists or UNSPECIFIED if not
    unsigned long EquivalentEdgeInCSMP( size_t ug_idx ) const;

    /// returns CSMP line element/face/interface  idx if it exists or UNSPECIFIED if not
    unsigned long EquivalentFaceInCSMP( size_t ug_idx ) const;
    
    /// using the node numbering in the region 'Model' the property values at the vertices are written out
    void WriteVertexVariableValue( std::ofstream&, const csmp::Index&, const Region<dim>&, size_t ug_cell_idx, bool& print_whitespace ) const;

    void WriteEdgeVariableValue( std::ofstream&, const csmp::Index&, const Region<dim>&, size_t ug_cell_idx, bool& print_whitespace ) const;
    ///
    void WriteFaceVariableValue( std::ofstream&, const csmp::Index&, const MeshManager<dim>&, size_t ug_cell_idx, bool& print_whitespace ) const;
    ///
    void WriteVolumeVariableValue( std::ofstream&, const csmp::Index&, const MeshManager<dim>&, size_t ug_cell_idx, bool& print_whitespace ) const;

    // translation of CSMP variable types to UG variable types
    std::string UG_VariableType( const csmp::Index& ) const;
        
    // extra datasets needed for UGX output; note that the set of global node ids gets sorted!
    // UG also does not require that nodes, edges, faces and volumes are listed in any particular numbering sense.
    // Thus, sorting allows to use them as search keys.
    // Note also that the UGX format only supports elements with linear shape functions, higher-order elements are created inside of UG
    // global node ids, followed by element ID of edge
    std::vector<std::pair<size_t,size_t>> edges_; ///< global edge end-node IDs for all edges in the model as needed by UG (n...total)
    //     csmp-idx, edge-#
    std::map<size_t,size_t>  csmp_edges_, csmp_edges_inverted_; ///< mapping between CSMP element/face IDs and edge number in UG

    std::vector<std::pair<std::set<size_t>,std::array<size_t,3>>>  tria_faces_; ///< triangle faces numbered continuously from 0..faces-1, second original node order
    std::vector<std::pair<std::set<size_t>,std::array<size_t,4>>>  quad_faces_;
                                                     
    std::vector<std::pair<std::set<size_t>,size_t>>  tetra_volumes_,     ///< volumes = tets + hexa + prism + pyra, numbered continuously
                                                     hexa_volumes_,
                                                     prism_volumes_,
                                                     pyra_volumes_;      ///< volumetric elements present only in 3D models
    //     csmp-idx,ug-id
    std::map<size_t,size_t> csmp_faces_, csmp_faces_inverted_;           ///< mapping between CSMP element/face IDs and face-number in UG
  
    std::vector<bool> edges_in_csmp_, faces_in_csmp_; ///<  vectors that can be queried before a search of a CSMP object is carried out to see whether it exists at all
  
    // printing of region=subset ids (0..n-1) which requires searching the edges, faces, and volumes vectors
    bool with_region_edge_output_   = true;
    bool with_region_face_output_   = true;
    bool with_region_volume_output_ = true;

    bool with_boundary_edge_output_   = true;
    bool with_boundary_face_output_   = true;
};


} // end csmp

#endif /* UG4_UGX_FILE_EXPORT_H */
