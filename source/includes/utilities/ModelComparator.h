// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_MODEL_COMPARATOR_H
#define CSMP_MODEL_COMPARATOR_H

#include "Region.h"
#include "Point.h"
#include "PropertyAtPointVisitor.h"

namespace csmp {
  
  struct Index;
  class ScalarVariable;
  template<uint32_t> class Model;
  template<uint32_t> class VSet;


  /**
      @brief Interface to compare two instances of CSMP models.
    
      Comparisons assume identical meshes with the same numbering of nodes.

      Compares Models, Vsets. For usage example, please refer to ModelComparator_Test
      for the time being.
      
      @attention since VSet and VData have overloaded operator==() functions these provide an alternative way to compare different polygonal datasets.

      @todo (3) Different error norms
      @todo (3) Different properties of same model
      @todo (3) Cor diff var types
      @todo (3) Cefactor(outsource) functions once diff var types implemented
      @todo (3) Cefactor static, no sense anymore - DONE(or not so much???)
      @todo (3) Cefaults for var file param
      @todo (3) Compare list of properties

     @author Philip Lang
     @date 2011
  */
  template<uint32_t dim>
  class ModelComparator
  {
  public:

    /////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for identical meshes with identical numbering of nodes
    /////////////////////////////////////////////////////////////////////////////////////////

    /// do the properties (named parameters) with the same type and placement have the same values in both models
    double CompareModels( const Model<dim>& model1, const Model<dim>& model2,
                          const char* property1, const char* property2, const char* region = "Model" ) const;

    /// using the polgonal mesh and data container VSet, this method compares the values of the specified properties among them; properties must have same placement and type
    double CompareVSets( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                         const char* variablesFile1, const char* variablesFile2, const char* region = "Model" ) const;
    
    /// reads input VSets from file
    double CompareVSets( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                         const char* variablesFile1, const char* variablesFile2, const char* region = "Model" ) const;


    ////////////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for identical meshes but does not depend on node numbering
    ////////////////////////////////////////////////////////////////////////////////////////////////

    double CompareModelsRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                         const char* property1, const char* property2, const char* region = "Model" );

    double CompareVSetsRenumberedNodes( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                                        const char* variablesFile1, const char* variablesFile2, const char* region = "Model" );
    
    double CompareVSetsRenumberedNodes( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                                        const char* variablesFile1, const char* variablesFile2, const char* region = "Model" );


    /////////////////////////////////////////////////////////////////////////////////////////
    /// Comparison is done for meshes that can be in general different from each other     
    /////////////////////////////////////////////////////////////////////////////////////////

    double CompareModelsAtPoints( Model<dim>& model1, Model<dim>& model2,
                                  const char* property1, const char* property2, const char* region = "Model" );

    double CompareVSetsAtPoints( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                                 const char* variablesFile1, const char* variablesFile2, const char* region = "Model" );
    
    double CompareVSetsAtPoints( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                                 const char* variablesFile1, const char* variablesFile2, const char* region = "Model" );

  private:

    double CompareRegionScalarVariable( const Model<dim>& model1, const Model<dim>& model2,
                                        const Index model1Key, const Index model2Key, const char* region = "Model" ) const;

    double CompareScalarDequesL2( const std::deque<ScalarVariable>& deque1,
                                  const std::deque<ScalarVariable>& deque2 ) const;

    void ReadRegionElementScalarVariable( const Region<dim>& region,
                                          const Index propKey,
                                          std::deque<ScalarVariable>& scalarDeque ) const;

    void ReadRegionNodalScalarVariable( const Region<dim>& region,
                                        const Index propKey,
                                        std::deque<ScalarVariable>& scalarDeque ) const;

    double CompareRegionScalarVariableRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                                       const char* property1, const char* property2,
                                                       const char* region = "Model" );

    void ReadRegionNodalScalarVariableAndNodeCoordinates( const Region<dim>& region,
                                                          const Index propKey,
                                                          std::map<Point<dim>, ScalarVariable >& points_and_values );

    double CompareRegionScalarVariableAtPoints( Model<dim>& model1, Model<dim>& model2,
                                                const char* property1, const char* property2, const char* region = "Model" );

    void ReadRegionNodalScalarVariableAndNodeCoordinates( const Region<dim>& region,
                                                          const Index propKey,
                                                          std::vector<ScalarVariable> & values,
                                                          std::map<size_t, std::vector<double> >& points);
  };

} // csmp

#endif
