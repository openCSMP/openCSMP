//
//  quadraturePointOutputToVTK.cpp
//  CSMP_API_examples
//
//  Created by Stephan Matthai on 9/06/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "CSMP_definitions.h"
#include "QuadraturePointOutputToVTK.h"

using namespace std;

namespace csmp {

/**
    Extrapolates scalar integration point properties to unique (new) nodes of quadratic tetrahedra and triangles.
    
    The legacy VTK ASCII text file format is used.

    Quadratic triangular elements only.
*/


	QuadraturePointOutputToVTK::QuadraturePointOutputToVTK()
	{
	}

	QuadraturePointOutputToVTK::~QuadraturePointOutputToVTK()
	{
	}


	void QuadraturePointOutputToVTK::OutputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK(const Model<3U>& model, const char* region, const char* variable_name)
 {
     const size_t dim(3U); // REMOVE WHEN TEMPLATIZING
     const Region<dim>& model_domain(model.Region(region));
     const csmp::Index key(model.Database().StorageKey(variable_name));
     assert( key.type == SCALAR );
   
     if ( key.place != ELEMENT_INTEGRATION_POINT ) {
          ErrorHandler::Instance().notice( CSMP_ERROR, "outputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK:",
                                          "output variable must be placed on the element integration point; nothing was done." );
          return;
       }
     if ( !(*model_domain.ElementsBegin())->FE()->UsesLocalCoordinates() ) {
          ErrorHandler::Instance().notice( CSMP_ERROR, "outputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK:",
                                          "the mesh must consist of numerically integrated elements; nothing was done." );
          return;
       }

//const string file_name(string(model.Name()) + "_" + region + "_" + variable_name);
     const string file_name( variable_name );
     string var_name(variable_name);
     replaceWhiteSpaceBy( var_name, '_' );

     const string  outfile(file_name + ".vtk");

     // 0. opening data output file in ascii format
     // -------------------------------------------
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs ) {
          ErrorHandler::Instance().notice( CSMP_ERROR, "outputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK:",
                                           outfile, "output file could not be opened; nothing was done." );
          return;
       }

     // 1. writing the VTK file header
     // ------------------------------
     ofs <<"# vtk DataFile Version 2.0"<< endl;
     ofs <<"Finite-element dataset (CSMP): variable: "<< var_name << endl;
     ofs <<"ASCII"<< endl << endl;
   
   
     // 2. VTK point data (each element is assigned unique nodes, some of which coincide in position)
     // ---------------------------------------------------------------------------------------------
     ofs <<"DATASET UNSTRUCTURED_GRID"<< endl;
     size_t nodes(0u), cells(0U);
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
          nodes += (*it)->Nodes();
          cells++;
       }
     // NOTE: in VTK coordinates always are stored in single precision
     ofs <<"POINTS " << nodes <<" float"<< endl;

     DenseMatrix<DM_MIN> COORD;
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
          (*it)->NodeCoordinateMatrix( COORD );
          for ( size_t i=0; i<(*it)->Nodes(); ++i ) {
              for ( size_t j=0; j<dim; j++ ) ofs << COORD(i,j) <<" ";
              if ( dim == 2 ) ofs << 0.;
              ofs << endl;
            }
          ofs << endl;
       }
   
     // 3. writing CELLS (cell-size and member nodes (point))
     // -----------------------------------------------------
     ofs <<"CELLS "<< cells <<" "<< nodes + cells << endl;
     size_t node(0U);
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
            // nodes per cell
            ofs << (*it)->Nodes() <<" ";
            // member nodes (running node index)
            for ( size_t i=0; i<(*it)->Nodes(); ++i ) ofs << node++ <<" ";
            ofs << endl;
        }

     // 4. writing CELL_TYPES
     // ---------------------
     ofs <<"CELL_TYPES "<< cells << endl;
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it )
       if ( (*it)->FE_Type() == ISOPARAMETRIC_QUADRATIC_TRIANGLE ||
            (*it)->FE_Type() == ISOPARAMETRIC_QUADRATIC_TETRAHEDRON )
         ofs << parseElementType( (*it)->FE_Type() ) << endl;
     ofs << endl;
   
     // 5. writing POINT_DATA point-type data values
     // --------------------------------------------
     // Unfortunately the data can only be output as nodal variables
     ofs <<"POINT_DATA "<< nodes << endl;
     ofs.setf( ios::scientific );
     ofs <<"SCALARS "<< var_name <<" double"<< endl;
     ofs <<"LOOKUP_TABLE default" << endl; // table must always be created
   
     // extrapolating the integration point data to the node points and writing them to file
     vector<double64> IVAR, NVAR;
     for ( auto it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
           // collecting the variable values from the integration points
           IVAR.resize( (*it)->IntegrationPoints() );
           for ( size_t i=0U; i<(*it)->IntegrationPoints(); ++i )
             IVAR[i] = (*it)->Read( i, key );
           // extrapolation
           const size_t n_node_variables(1U);
           (*it)->FE()->ExtrapolateIntegrationPointVariableToNodes( n_node_variables, IVAR, NVAR );
           for ( size_t i=0U; i<NVAR.size(); ++i ) ofs << NVAR[i] <<" ";
           ofs << endl;
       }

 } // end outputQuadraturePointPropertiesAsDiscontinuousNodeVariablesToVTK

} // end csmp
